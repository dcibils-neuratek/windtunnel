/* ------------------------------------------------------------------
   lbm.js — 3D Lattice Boltzmann fluid solver (D3Q19)

   This is a real, incompressible Navier-Stokes solver. The lattice
   Boltzmann equation

       f_q(x + c_q dt, t + dt) = f_q(x,t) - (1/tau) [ f_q - f_q^eq ]

   recovers the incompressible Navier-Stokes equations in the low-Mach
   limit (Chapman-Enskog expansion), with kinematic viscosity

       nu = c_s^2 (tau - 1/2),   c_s^2 = 1/3   (lattice units)

   Features:
     - D3Q19 velocity set, BGK single-relaxation-time collision
     - Smagorinsky sub-grid model (LES) so we stay stable at high Re
     - Halfway bounce-back no-slip walls on the solid body
     - Momentum-exchange force evaluation -> real drag / lift / side force

   Layout note: populations are stored structure-of-arrays,
   f[q * n + cell], so the streaming gather for a fixed direction walks
   memory sequentially. The hot kernel below is fully unrolled over the
   19 directions - ugly, but roughly 3x faster than the loop version.
------------------------------------------------------------------ */
(function (global) {
  'use strict';

  const Q = 19;
  //                          0  1  2  3  4  5  6   7  8  9 10  11 12 13 14  15 16 17 18
  const EX = new Int32Array([ 0, 1,-1, 0, 0, 0, 0,  1,-1, 1,-1,  1,-1, 1,-1,  0, 0, 0, 0]);
  const EY = new Int32Array([ 0, 0, 0, 1,-1, 0, 0,  1,-1,-1, 1,  0, 0, 0, 0,  1,-1, 1,-1]);
  const EZ = new Int32Array([ 0, 0, 0, 0, 0, 1,-1,  0, 0, 0, 0,  1,-1,-1, 1,  1,-1,-1, 1]);
  const OPP = new Int32Array([0, 2, 1, 4, 3, 6, 5,  8, 7,10, 9, 12,11,14,13, 16,15,18,17]);
  const W = new Float64Array(Q);
  W[0] = 1 / 3;
  for (let q = 1; q <= 6; q++) W[q] = 1 / 18;
  for (let q = 7; q < Q; q++) W[q] = 1 / 36;

  const FLUID = 0, SOLID_OBJ = 1, SOLID_WALL = 2;
  const W0 = 1 / 3, WF = 1 / 18, WE = 1 / 36;
  const SQ18 = 18 * Math.SQRT2;

  class LBM {
    constructor(nx, ny, nz) {
      this.nx = nx; this.ny = ny; this.nz = nz;
      this.n = nx * ny * nz;
      const n = this.n;

      this.f = new Float32Array(n * Q);
      this.g = new Float32Array(n * Q);
      this.solid = new Uint8Array(n);
      this.near = new Uint8Array(n);

      this.ux = new Float32Array(n);
      this.uy = new Float32Array(n);
      this.uz = new Float32Array(n);
      this.rho = new Float32Array(n);

      // surface velocity of solid cells — a rolling road and spinning wheels
      // are no-slip walls that happen to be moving
      this.uwx = new Float32Array(n);
      this.uwy = new Float32Array(n);
      this.uwz = new Float32Array(n);
      this.movingWalls = false;

      this.off = new Int32Array(Q);   // cell offset of c_q
      this.src = new Int32Array(Q);   // read base:  f[src[q] + i]
      this.dst = new Int32Array(Q);   // write base: g[dst[q] + i]
      for (let q = 0; q < Q; q++) {
        this.off[q] = EX[q] + nx * (EY[q] + ny * EZ[q]);
        this.src[q] = q * n - this.off[q];
        this.dst[q] = q * n;
      }

      this.u0 = 0.06;
      this.tau = 0.5035;
      this.les = true;
      this.csm = 0.14;

      this.force = new Float64Array(3);
      this.frontalArea = 1;
      this.solidCells = 0;
      this.ground = false;
      this.steps = 0;

      this._fq = new Float64Array(Q);
      this._wv = new Float64Array(3);
      this._col = new Uint8Array(ny * nz);   // frontal (x) projection
      this._plan = new Uint8Array(nx * nz);  // planform (y) projection
      this._low = new Int16Array(nx * nz);   // lowest body cell per column
      this.rideHeight = 0;

      this.reset();
    }

    idx(x, y, z) { return x + this.nx * (y + this.ny * z); }
    get nu() { return (this.tau - 0.5) / 3; }

    _writeEq(arr, i, r, vx, vy, vz) {
      const n = this.n;
      const usq = 1.5 * (vx * vx + vy * vy + vz * vz);
      for (let q = 0; q < Q; q++) {
        const cu = 3 * (EX[q] * vx + EY[q] * vy + EZ[q] * vz);
        arr[q * n + i] = W[q] * r * (1 + cu + 0.5 * cu * cu - usq);
      }
    }

    reset() {
      const { n, f, g, solid, ux, uy, uz, rho } = this;
      for (let i = 0; i < n; i++) {
        const v = solid[i] ? 0 : this.u0;
        this._writeEq(f, i, 1, v, 0, 0);
        this._writeEq(g, i, 1, v, 0, 0);
        ux[i] = v; uy[i] = 0; uz[i] = 0; rho[i] = 1;
      }
      this.steps = 0;
      this.force[0] = this.force[1] = this.force[2] = 0;
    }

    _seedCell(i) {
      this._writeEq(this.f, i, 1, this.u0, 0, 0);
      this._writeEq(this.g, i, 1, this.u0, 0, 0);
      this.ux[i] = this.u0; this.uy[i] = 0; this.uz[i] = 0; this.rho[i] = 1;
    }

    /**
     * Rasterise the body into the lattice.
     * @param inside  fn(x,y,z)->bool in lattice coordinates
     * @param box     {x0,y0,z0,x1,y1,z1}
     * @param ground  add a no-slip road at y = 0
     * @param wallVel optional fn(x,y,z,out)->bool giving the surface velocity
     *                of a body cell (rotating wheels); out is [vx,vy,vz]
     * @param roadVel optional streamwise speed of the road surface
     */
    voxelize(inside, box, ground, wallVel, roadVel) {
      const { nx, ny, nz, n, solid } = this;
      const uwx = this.uwx, uwy = this.uwy, uwz = this.uwz;
      this.ground = !!ground;
      this.movingWalls = false;

      for (let i = 0; i < n; i++) {
        if (solid[i] === SOLID_OBJ) {
          solid[i] = FLUID; this._seedCell(i);
          uwx[i] = uwy[i] = uwz[i] = 0;
        }
      }

      const gy = ground ? SOLID_WALL : FLUID;
      const rv = roadVel || 0;
      if (rv) this.movingWalls = true;
      for (let z = 0; z < nz; z++) {
        for (let x = 0; x < nx; x++) {
          const i = x + nx * ny * z;
          if (solid[i] !== gy) { solid[i] = gy; if (gy === FLUID) this._seedCell(i); }
          uwx[i] = gy === SOLID_WALL ? rv : 0;
          uwy[i] = uwz[i] = 0;
        }
      }

      const x0 = Math.max(1, Math.floor(box.x0)), x1 = Math.min(nx - 2, Math.ceil(box.x1));
      const y0 = Math.max(1, Math.floor(box.y0)), y1 = Math.min(ny - 2, Math.ceil(box.y1));
      const z0 = Math.max(1, Math.floor(box.z0)), z1 = Math.min(nz - 2, Math.ceil(box.z1));

      const col = this._col; col.fill(0);
      const plan = this._plan; plan.fill(0);
      const low = this._low; low.fill(32767);
      const wv = this._wv;
      let cells = 0;
      for (let z = z0; z <= z1; z++) {
        for (let y = y0; y <= y1; y++) {
          let i = x0 + nx * (y + ny * z);
          for (let x = x0; x <= x1; x++, i++) {
            if (inside(x + 0.5, y + 0.5, z + 0.5)) {
              solid[i] = SOLID_OBJ;
              this._writeEq(this.f, i, 1, 0, 0, 0);
              this.ux[i] = this.uy[i] = this.uz[i] = 0; this.rho[i] = 1;
              uwx[i] = uwy[i] = uwz[i] = 0;
              if (wallVel && wallVel(x + 0.5, y + 0.5, z + 0.5, wv)) {
                uwx[i] = wv[0]; uwy[i] = wv[1]; uwz[i] = wv[2];
                this.movingWalls = true;
              }
              col[y + ny * z] = 1;
              const pk = x + nx * z;
              plan[pk] = 1;
              if (y < low[pk]) low[pk] = y;
              cells++;
            }
          }
        }
      }
      let a = 0;
      for (let k = 0; k < col.length; k++) a += col[k];
      this.frontalArea = Math.max(1, a);
      let p = 0;
      for (let k = 0; k < plan.length; k++) p += plan[k];
      this.planformArea = Math.max(1, p);
      this.solidCells = cells;

      // Underbody gap, in cells. Columns that sit on the road (wheels) give
      // zero, so the median of the non-zero gaps is the ride height. A gap of
      // only one or two cells cannot support a shear layer, which makes a
      // moving road wildly over-predict.
      const gaps = [];
      for (let k = 0; k < plan.length; k++) {
        if (plan[k] && low[k] < 32767) { const g = low[k] - 1; if (g > 0) gaps.push(g); }
      }
      gaps.sort((a, b) => a - b);
      this.rideHeight = gaps.length ? gaps[gaps.length >> 1] : 0;
      this._rebuildNear();
    }

    _rebuildNear() {
      const { nx, ny, nz, solid, near, off } = this;
      near.fill(0);
      for (let z = 1; z < nz - 1; z++) {
        for (let y = 1; y < ny - 1; y++) {
          let i = 1 + nx * (y + ny * z);
          for (let x = 1; x < nx - 1; x++, i++) {
            if (solid[i]) continue;
            for (let q = 1; q < Q; q++) {
              if (solid[i - off[q]]) { near[i] = 1; break; }
            }
          }
        }
      }
    }

    /** Advance the flow by one lattice time step. */
    step() {
      const nx = this.nx, ny = this.ny, nz = this.nz, nn = this.n;
      const F = this.f, G = this.g;
      const solid = this.solid, near = this.near, off = this.off;
      const src = this.src, dst = this.dst;
      const ux = this.ux, uy = this.uy, uz = this.uz, rho = this.rho;
      const uwx = this.uwx, uwy = this.uwy, uwz = this.uwz, mov = this.movingWalls;
      const scratch = this._fq;
      const tau0 = this.tau, les = this.les, cs2 = this.csm * this.csm;
      const omega0 = 1 / tau0;

      const s0 = src[0], s1 = src[1], s2 = src[2], s3 = src[3], s4 = src[4],
        s5 = src[5], s6 = src[6], s7 = src[7], s8 = src[8], s9 = src[9],
        s10 = src[10], s11 = src[11], s12 = src[12], s13 = src[13], s14 = src[14],
        s15 = src[15], s16 = src[16], s17 = src[17], s18 = src[18];
      const d0 = dst[0], d1 = dst[1], d2 = dst[2], d3 = dst[3], d4 = dst[4],
        d5 = dst[5], d6 = dst[6], d7 = dst[7], d8 = dst[8], d9 = dst[9],
        d10 = dst[10], d11 = dst[11], d12 = dst[12], d13 = dst[13], d14 = dst[14],
        d15 = dst[15], d16 = dst[16], d17 = dst[17], d18 = dst[18];

      let Fx = 0, Fy = 0, Fz = 0;

      for (let z = 1; z < nz - 1; z++) {
        for (let y = 1; y < ny - 1; y++) {
          let i = 1 + nx * (y + ny * z);
          for (let x = 1; x < nx - 1; x++, i++) {
            if (solid[i]) continue;

            let f0, f1, f2, f3, f4, f5, f6, f7, f8, f9,
              f10, f11, f12, f13, f14, f15, f16, f17, f18;

            if (near[i]) {
              // --- streaming with halfway bounce-back + momentum exchange ---
              for (let q = 0; q < Q; q++) {
                const j = i - off[q];
                const s = solid[j];
                if (s) {
                  const fd = F[OPP[q] * nn + i];      // population that hit the wall
                  // Moving no-slip wall (rolling road, spinning wheel): the
                  // reflected population is biased by the wall's tangential
                  // velocity, f_q = f_qbar + 2*w_q*rho*(c_q . u_w)/cs^2.
                  let corr = 0;
                  if (mov) {
                    const cu = EX[q] * uwx[j] + EY[q] * uwy[j] + EZ[q] * uwz[j];
                    if (cu !== 0) corr = 6 * W[q] * rho[i] * cu;
                  }
                  scratch[q] = fd + corr;
                  if (s === SOLID_OBJ) {
                    // Momentum exchange c_qbar*(f_in + f_out), measured
                    // relative to the rest equilibrium w_q*rho0. Over a closed
                    // surface that reference sums to zero, so it changes
                    // nothing there (and kills a nasty float32 cancellation) -
                    // but for a body sitting on the road the surface is NOT
                    // closed, and without it the unbalanced ambient pressure
                    // swamps the aerodynamic load.
                    const d = 2 * (fd - W[q]) + corr;
                    Fx -= d * EX[q]; Fy -= d * EY[q]; Fz -= d * EZ[q];
                  }
                } else {
                  scratch[q] = F[src[q] + i];
                }
              }
              f0 = scratch[0]; f1 = scratch[1]; f2 = scratch[2]; f3 = scratch[3];
              f4 = scratch[4]; f5 = scratch[5]; f6 = scratch[6]; f7 = scratch[7];
              f8 = scratch[8]; f9 = scratch[9]; f10 = scratch[10]; f11 = scratch[11];
              f12 = scratch[12]; f13 = scratch[13]; f14 = scratch[14]; f15 = scratch[15];
              f16 = scratch[16]; f17 = scratch[17]; f18 = scratch[18];
            } else {
              f0 = F[s0 + i]; f1 = F[s1 + i]; f2 = F[s2 + i]; f3 = F[s3 + i];
              f4 = F[s4 + i]; f5 = F[s5 + i]; f6 = F[s6 + i]; f7 = F[s7 + i];
              f8 = F[s8 + i]; f9 = F[s9 + i]; f10 = F[s10 + i]; f11 = F[s11 + i];
              f12 = F[s12 + i]; f13 = F[s13 + i]; f14 = F[s14 + i]; f15 = F[s15 + i];
              f16 = F[s16 + i]; f17 = F[s17 + i]; f18 = F[s18 + i];
            }

            // --- moments ---
            const r = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8 + f9 +
              f10 + f11 + f12 + f13 + f14 + f15 + f16 + f17 + f18;
            const inv = 1 / r;
            const vx = (f1 - f2 + f7 - f8 + f9 - f10 + f11 - f12 + f13 - f14) * inv;
            const vy = (f3 - f4 + f7 - f8 - f9 + f10 + f15 - f16 + f17 - f18) * inv;
            const vz = (f5 - f6 + f11 - f12 - f13 + f14 + f15 - f16 - f17 + f18) * inv;
            rho[i] = r; ux[i] = vx; uy[i] = vy; uz[i] = vz;

            // --- equilibrium (O(u^2) truncation) ---
            const A = 1 - 1.5 * (vx * vx + vy * vy + vz * vz);
            const t0 = r * W0, t1 = r * WF, t2 = r * WE;
            const ax = vx, ay = vy, az = vz;
            const bxy = vx + vy, bxY = vx - vy;
            const bxz = vx + vz, bxZ = vx - vz;
            const byz = vy + vz, byZ = vy - vz;

            const e0 = t0 * A;
            const hx = 4.5 * ax * ax, hy = 4.5 * ay * ay, hz = 4.5 * az * az;
            const e1 = t1 * (A + 3 * ax + hx), e2 = t1 * (A - 3 * ax + hx);
            const e3 = t1 * (A + 3 * ay + hy), e4 = t1 * (A - 3 * ay + hy);
            const e5 = t1 * (A + 3 * az + hz), e6 = t1 * (A - 3 * az + hz);
            const h7 = 4.5 * bxy * bxy, h9 = 4.5 * bxY * bxY;
            const h11 = 4.5 * bxz * bxz, h13 = 4.5 * bxZ * bxZ;
            const h15 = 4.5 * byz * byz, h17 = 4.5 * byZ * byZ;
            const e7 = t2 * (A + 3 * bxy + h7), e8 = t2 * (A - 3 * bxy + h7);
            const e9 = t2 * (A + 3 * bxY + h9), e10 = t2 * (A - 3 * bxY + h9);
            const e11 = t2 * (A + 3 * bxz + h11), e12 = t2 * (A - 3 * bxz + h11);
            const e13 = t2 * (A + 3 * bxZ + h13), e14 = t2 * (A - 3 * bxZ + h13);
            const e15 = t2 * (A + 3 * byz + h15), e16 = t2 * (A - 3 * byz + h15);
            const e17 = t2 * (A + 3 * byZ + h17), e18 = t2 * (A - 3 * byZ + h17);

            // --- non-equilibrium parts ---
            const n1 = f1 - e1, n2 = f2 - e2, n3 = f3 - e3, n4 = f4 - e4;
            const n5 = f5 - e5, n6 = f6 - e6, n7 = f7 - e7, n8 = f8 - e8;
            const n9 = f9 - e9, n10 = f10 - e10, n11 = f11 - e11, n12 = f12 - e12;
            const n13 = f13 - e13, n14 = f14 - e14, n15 = f15 - e15, n16 = f16 - e16;
            const n17 = f17 - e17, n18 = f18 - e18;

            // --- relaxation rate (+ Smagorinsky eddy viscosity) ---
            let om = omega0;
            if (les) {
              const pxx = n1 + n2 + n7 + n8 + n9 + n10 + n11 + n12 + n13 + n14;
              const pyy = n3 + n4 + n7 + n8 + n9 + n10 + n15 + n16 + n17 + n18;
              const pzz = n5 + n6 + n11 + n12 + n13 + n14 + n15 + n16 + n17 + n18;
              const pxy = n7 + n8 - n9 - n10;
              const pxz = n11 + n12 - n13 - n14;
              const pyz = n15 + n16 - n17 - n18;
              const mag = Math.sqrt(2 * (pxx * pxx + pyy * pyy + pzz * pzz +
                2 * (pxy * pxy + pxz * pxz + pyz * pyz)));
              om = 2 / (tau0 + Math.sqrt(tau0 * tau0 + SQ18 * cs2 * mag * inv));
            }

            // --- collision ---
            G[d0 + i] = f0 - om * (f0 - e0);
            G[d1 + i] = f1 - om * n1; G[d2 + i] = f2 - om * n2;
            G[d3 + i] = f3 - om * n3; G[d4 + i] = f4 - om * n4;
            G[d5 + i] = f5 - om * n5; G[d6 + i] = f6 - om * n6;
            G[d7 + i] = f7 - om * n7; G[d8 + i] = f8 - om * n8;
            G[d9 + i] = f9 - om * n9; G[d10 + i] = f10 - om * n10;
            G[d11 + i] = f11 - om * n11; G[d12 + i] = f12 - om * n12;
            G[d13 + i] = f13 - om * n13; G[d14 + i] = f14 - om * n14;
            G[d15 + i] = f15 - om * n15; G[d16 + i] = f16 - om * n16;
            G[d17 + i] = f17 - om * n17; G[d18 + i] = f18 - om * n18;
          }
        }
      }

      this._boundaries(G);
      this.f = G; this.g = F;
      this.force[0] = Fx; this.force[1] = Fy; this.force[2] = Fz;
      this.steps++;
    }

    /** Velocity inlet, zero-gradient outlet, free-stream side walls. */
    _boundaries(g) {
      const nx = this.nx, ny = this.ny, nz = this.nz, n = this.n, u0 = this.u0;
      const ux = this.ux, uy = this.uy, uz = this.uz, rho = this.rho, solid = this.solid;

      const setEq = (i) => {
        this._writeEq(g, i, 1, u0, 0, 0);
        ux[i] = u0; uy[i] = 0; uz[i] = 0; rho[i] = 1;
      };

      for (let z = 0; z < nz; z++) {
        for (let y = 0; y < ny; y++) setEq(nx * (y + ny * z));              // inlet
      }
      for (let z = 0; z < nz; z++) {
        for (let x = 0; x < nx; x++) {
          setEq(x + nx * ((ny - 1) + ny * z));                              // ceiling
          const i0 = x + nx * ny * z;
          if (!solid[i0]) setEq(i0);                                        // floor
        }
      }
      for (let y = 0; y < ny; y++) {
        for (let x = 0; x < nx; x++) {
          setEq(x + nx * y);                                                // z = 0
          setEq(x + nx * (y + ny * (nz - 1)));                              // z = nz-1
        }
      }
      const xo = nx - 1;                                                    // outlet
      for (let z = 0; z < nz; z++) {
        for (let y = 0; y < ny; y++) {
          const i = xo + nx * (y + ny * z), j = i - 1;
          for (let q = 0; q < Q; q++) g[q * n + i] = g[q * n + j];
          ux[i] = ux[j]; uy[i] = uy[j]; uz[i] = uz[j]; rho[i] = rho[j];
        }
      }
    }

    /** Trilinear velocity sample at lattice position (x,y,z). */
    sample(x, y, z, out) {
      const nx = this.nx, ny = this.ny, nz = this.nz;
      if (x < 0.5 || y < 0.5 || z < 0.5 || x > nx - 1.5 || y > ny - 1.5 || z > nz - 1.5) {
        out[0] = this.u0; out[1] = 0; out[2] = 0; return out;
      }
      const xf = x - 0.5, yf = y - 0.5, zf = z - 0.5;
      const x0 = xf | 0, y0 = yf | 0, z0 = zf | 0;
      const tx = xf - x0, ty = yf - y0, tz = zf - z0;
      const ux = this.ux, uy = this.uy, uz = this.uz;
      let vx = 0, vy = 0, vz = 0;
      for (let k = 0; k < 2; k++) {
        const wz = k ? tz : 1 - tz;
        for (let j = 0; j < 2; j++) {
          const wy = (j ? ty : 1 - ty) * wz;
          const base = nx * ((y0 + j) + ny * (z0 + k));
          const ia = x0 + base, ib = ia + 1;
          const wa = (1 - tx) * wy, wb = tx * wy;
          vx += ux[ia] * wa + ux[ib] * wb;
          vy += uy[ia] * wa + uy[ib] * wb;
          vz += uz[ia] * wa + uz[ib] * wb;
        }
      }
      out[0] = vx; out[1] = vy; out[2] = vz;
      return out;
    }

    isSolid(x, y, z) {
      const xi = x | 0, yi = y | 0, zi = z | 0;
      if (xi < 0 || yi < 0 || zi < 0 || xi >= this.nx || yi >= this.ny || zi >= this.nz) return false;
      return this.solid[xi + this.nx * (yi + this.ny * zi)] !== 0;
    }
  }

  LBM.Q = Q; LBM.EX = EX; LBM.EY = EY; LBM.EZ = EZ; LBM.W = W;
  LBM.FLUID = FLUID; LBM.SOLID_OBJ = SOLID_OBJ; LBM.SOLID_WALL = SOLID_WALL;
  global.LBM = LBM;
})(window);
