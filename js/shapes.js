/* ------------------------------------------------------------------
   shapes.js — test bodies

   Every body is a small list of analytic primitives defined in a
   normalised space where the body length spans x in [-0.5, 0.5],
   y is up (y = 0 is the ground contact plane for wheeled bodies) and
   z is spanwise.

   The SAME primitive list drives both the rendered geometry and the
   lattice voxelisation, so what you see is exactly what the solver
   feels — no mesh/collision mismatch.
------------------------------------------------------------------ */
(function (global) {
  'use strict';

  // ---- NACA 4-digit aerofoil ----------------------------------------
  function nacaThickness(xc, t) {
    return 5 * t * (0.2969 * Math.sqrt(xc) - 0.1260 * xc - 0.3516 * xc * xc +
      0.2843 * xc * xc * xc - 0.1036 * xc * xc * xc * xc);
  }
  function nacaCamber(xc, m, p) {
    if (m === 0) return 0;
    return xc < p
      ? (m / (p * p)) * (2 * p * xc - xc * xc)
      : (m / ((1 - p) * (1 - p))) * ((1 - 2 * p) + 2 * p * xc - xc * xc);
  }

  // ---- 3x3 rotation from Euler XYZ (matches THREE.Euler 'XYZ') ------
  function eulerMat(rx, ry, rz) {
    const a = Math.cos(rx), b = Math.sin(rx);
    const c = Math.cos(ry), d = Math.sin(ry);
    const e = Math.cos(rz), f = Math.sin(rz);
    const ae = a * e, af = a * f, be = b * e, bf = b * f;
    return [
      c * e, -c * f, d,
      af + be * d, ae - bf * d, -b * c,
      bf - ae * d, be + af * d, a * c
    ];
  }

  // ---- primitive point test (point already in primitive space) ------
  function insidePrim(p, q) {
    const s = p.s;
    switch (p.t) {
      case 'box':
        return Math.abs(q[0]) <= s[0] && Math.abs(q[1]) <= s[1] && Math.abs(q[2]) <= s[2];
      case 'ell': {
        const a = q[0] / s[0], b = q[1] / s[1], c = q[2] / s[2];
        return a * a + b * b + c * c <= 1;
      }
      case 'cyl':                                  // axis = local +y
        return Math.abs(q[1]) <= s[1] && (q[0] * q[0] + q[2] * q[2]) <= s[0] * s[0];
      case 'cone': {                               // apex at +y, base at -y
        if (Math.abs(q[1]) > s[1]) return false;
        const r = s[0] * (s[1] - q[1]) / (2 * s[1]);
        return (q[0] * q[0] + q[2] * q[2]) <= r * r;
      }
      case 'foil': {                               // chord s[0] along x, span s[2] along z
        const ch = s[0], span = s[2];
        if (Math.abs(q[2]) > span * 0.5) return false;
        const xc = (q[0] + ch * 0.5) / ch;
        if (xc < 0 || xc > 1) return false;
        const t = p.naca ? p.naca[2] : 0.12;
        const m = p.naca ? p.naca[0] : 0;
        const pp = p.naca ? p.naca[1] : 0.4;
        const yt = nacaThickness(xc, t) * ch;
        const yc = nacaCamber(xc, m, pp) * ch;
        return Math.abs(q[1] - yc) <= yt;
      }
    }
    return false;
  }

  // ---- body catalogue ------------------------------------------------
  // s = half-extents (box), radii (ell), [radius, halfHeight] (cyl/cone),
  //     [chord, -, span] (foil)
  const BODIES = {
    sedan: {
      name: 'Sedan', ground: true, lengthM: 4.6, refCd: '0.30 – 0.35',
      parts: [
        { t: 'ell', c: [0.00, 0.155, 0], s: [0.500, 0.090, 0.185] },
        { t: 'ell', c: [0.05, 0.235, 0], s: [0.250, 0.080, 0.160] },
        { t: 'box', c: [0.00, 0.120, 0], s: [0.430, 0.055, 0.180] },
        { t: 'cyl', wheel: true, c: [-0.290, 0.075, 0.175], s: [0.075, 0.035], r: [Math.PI / 2, 0, 0] },
        { t: 'cyl', wheel: true, c: [-0.290, 0.075, -0.175], s: [0.075, 0.035], r: [Math.PI / 2, 0, 0] },
        { t: 'cyl', wheel: true, c: [0.280, 0.075, 0.175], s: [0.075, 0.035], r: [Math.PI / 2, 0, 0] },
        { t: 'cyl', wheel: true, c: [0.280, 0.075, -0.175], s: [0.075, 0.035], r: [Math.PI / 2, 0, 0] }
      ]
    },
    sports: {
      name: 'Sports car', ground: true, lengthM: 4.3, refCd: '0.29 – 0.36',
      parts: [
        { t: 'ell', c: [0.00, 0.115, 0], s: [0.500, 0.070, 0.200] },
        { t: 'ell', c: [0.07, 0.180, 0], s: [0.190, 0.055, 0.150] },
        { t: 'box', c: [0.00, 0.100, 0], s: [0.450, 0.045, 0.195] },
        { t: 'box', c: [0.450, 0.245, 0], s: [0.045, 0.009, 0.170] },
        { t: 'box', c: [0.430, 0.205, 0.165], s: [0.020, 0.045, 0.008] },
        { t: 'box', c: [0.430, 0.205, -0.165], s: [0.020, 0.045, 0.008] },
        { t: 'cyl', wheel: true, c: [-0.300, 0.070, 0.195], s: [0.070, 0.045], r: [Math.PI / 2, 0, 0] },
        { t: 'cyl', wheel: true, c: [-0.300, 0.070, -0.195], s: [0.070, 0.045], r: [Math.PI / 2, 0, 0] },
        { t: 'cyl', wheel: true, c: [0.290, 0.070, 0.195], s: [0.070, 0.050], r: [Math.PI / 2, 0, 0] },
        { t: 'cyl', wheel: true, c: [0.290, 0.070, -0.195], s: [0.070, 0.050], r: [Math.PI / 2, 0, 0] }
      ]
    },
    truck: {
      name: 'Box truck', ground: true, lengthM: 7.0, refCd: '0.60 – 0.90',
      fill: 0.92,
      parts: [
        { t: 'box', c: [0.100, 0.330, 0], s: [0.380, 0.220, 0.200] },
        { t: 'box', c: [-0.340, 0.190, 0], s: [0.160, 0.110, 0.190] },
        { t: 'cyl', wheel: true, c: [-0.320, 0.080, 0.190], s: [0.080, 0.040], r: [Math.PI / 2, 0, 0] },
        { t: 'cyl', wheel: true, c: [-0.320, 0.080, -0.190], s: [0.080, 0.040], r: [Math.PI / 2, 0, 0] },
        { t: 'cyl', wheel: true, c: [0.300, 0.080, 0.190], s: [0.080, 0.040], r: [Math.PI / 2, 0, 0] },
        { t: 'cyl', wheel: true, c: [0.300, 0.080, -0.190], s: [0.080, 0.040], r: [Math.PI / 2, 0, 0] }
      ]
    },
    airplane: {
      name: 'Airliner', ground: false, lengthM: 38, refCd: '0.025 – 0.045',
      ref: 'planform', fill: 0.30,
      parts: [
        { t: 'ell', c: [0.00, 0, 0], s: [0.500, 0.055, 0.055] },
        { t: 'foil', c: [0.040, -0.020, 0], s: [0.290, 0, 0.900], naca: [0.03, 0.4, 0.15] },
        { t: 'foil', c: [0.430, 0.020, 0], s: [0.130, 0, 0.340], naca: [0.0, 0.4, 0.12] },
        { t: 'foil', c: [0.430, 0.110, 0], s: [0.160, 0, 0.200], naca: [0.0, 0.4, 0.12], r: [Math.PI / 2, 0, 0] },
        { t: 'ell', c: [0.010, -0.058, 0.230], s: [0.090, 0.035, 0.035] },
        { t: 'ell', c: [0.010, -0.058, -0.230], s: [0.090, 0.035, 0.035] }
      ]
    },
    wing: {
      name: 'Wing (NACA 4418)', ground: false, lengthM: 2.0, refCd: 'stalls near α ≈ 15°',
      ref: 'planform', fill: 0.42,
      parts: [
        { t: 'foil', c: [0, 0, 0], s: [0.950, 0, 0.900], naca: [0.04, 0.4, 0.18] }
      ]
    },
    sphere: {
      name: 'Sphere', ground: false, lengthM: 1.0, refCd: '0.47 (sub-critical)',
      fill: 0.79,
      parts: [{ t: 'ell', c: [0, 0, 0], s: [0.28, 0.28, 0.28] }]
    },
    cube: {
      name: 'Cube', ground: false, lengthM: 1.0, refCd: '1.05',
      fill: 1.0,
      parts: [{ t: 'box', c: [0, 0, 0], s: [0.22, 0.22, 0.22] }]
    },
    cylinder: {
      name: 'Cylinder (spanwise)', ground: false, lengthM: 1.0, refCd: '1.17',
      fill: 1.0,
      parts: [{ t: 'cyl', c: [0, 0, 0], s: [0.25, 0.45], r: [Math.PI / 2, 0, 0] }]
    },
    teardrop: {
      name: 'Teardrop', ground: false, lengthM: 2.0, refCd: '0.04',
      fill: 0.79,
      parts: [
        { t: 'ell', c: [-0.220, 0, 0], s: [0.280, 0.130, 0.130] },
        { t: 'cone', c: [0.130, 0, 0], s: [0.130, 0.350], r: [0, 0, -Math.PI / 2] }
      ]
    }
  };

  // ---- uploaded models ------------------------------------------------
  /**
   * Guess which model axis is "up" and which is "forward".
   *
   * Bounding-box extents are a poor signal — an airliner's span and length
   * are within a few percent, so ranking by size is a coin flip. Silhouette
   * area is far better: a vehicle is *shaped* to present the least area along
   * its direction of travel and the most when viewed from above.
   *
   * For a closed surface the projected silhouette along d is exactly
   * ½·Σ|Aᵢ·d| over triangle area vectors — every ray enters once and exits
   * once. That is one cheap pass with no rasterisation, and it degrades
   * gracefully on open meshes.
   *
   * Falls back to extent ranking when the three projections are too close to
   * call (a sphere or cube has no meaningful "forward").
   */
  function guessAxes(size, P) {
    const byExtent = () => {
      const o = [0, 1, 2].sort((a, b) => size[b] - size[a]);
      return { forward: o[0], up: o[2] };
    };
    if (!P || P.length < 9) return byExtent();

    const proj = [0, 0, 0];
    for (let i = 0; i < P.length; i += 9) {
      const ux = P[i + 3] - P[i], uy = P[i + 4] - P[i + 1], uz = P[i + 5] - P[i + 2];
      const vx = P[i + 6] - P[i], vy = P[i + 7] - P[i + 1], vz = P[i + 8] - P[i + 2];
      proj[0] += Math.abs(uy * vz - uz * vy);       // 2 * area vector components
      proj[1] += Math.abs(uz * vx - ux * vz);
      proj[2] += Math.abs(ux * vy - uy * vx);
    }
    const lo = Math.min(proj[0], proj[1], proj[2]);
    const hi = Math.max(proj[0], proj[1], proj[2]);
    if (!(lo > 0) || hi / lo < 1.15) return byExtent();   // too symmetric to tell

    const forward = proj.indexOf(lo);                     // least frontal area
    let up = proj.indexOf(hi);                            // biggest planform
    if (up === forward) up = [0, 1, 2].find(a => a !== forward);
    return { forward, up };
  }

  /**
   * Rotate/scale a triangle soup into body space: nose at -x, up at +y,
   * length normalised to 1 so it matches the built-in bodies.
   * @param P        Float32Array (9 per triangle), modified in place
   * @param forward  model axis index the nose points along
   * @param up       model axis index that points up
   * @param flip     nose points along -forward instead of +forward
   * @param ground   rest the model on y = 0 instead of centring it
   */
  function orient(P, forward, up, flip, ground) {
    // forward and up must be different axes, or f x u collapses and the
    // model flattens into a plane
    if (forward === up) forward = (up + 1) % 3;
    const f = [0, 0, 0], u = [0, 0, 0];
    f[forward] = flip ? -1 : 1;
    u[up] = 1;
    // side = f x u ; body rows are (-f, u, -side)
    const s = [
      f[1] * u[2] - f[2] * u[1],
      f[2] * u[0] - f[0] * u[2],
      f[0] * u[1] - f[1] * u[0]
    ];
    const R = [-f[0], -f[1], -f[2], u[0], u[1], u[2], -s[0], -s[1], -s[2]];

    let x0 = Infinity, y0 = Infinity, z0 = Infinity;
    let x1 = -Infinity, y1 = -Infinity, z1 = -Infinity;
    for (let i = 0; i < P.length; i += 3) {
      const a = P[i], b = P[i + 1], c = P[i + 2];
      const x = R[0] * a + R[1] * b + R[2] * c;
      const y = R[3] * a + R[4] * b + R[5] * c;
      const z = R[6] * a + R[7] * b + R[8] * c;
      P[i] = x; P[i + 1] = y; P[i + 2] = z;
      if (x < x0) x0 = x; if (x > x1) x1 = x;
      if (y < y0) y0 = y; if (y > y1) y1 = y;
      if (z < z0) z0 = z; if (z > z1) z1 = z;
    }
    const k = 1 / Math.max(1e-9, x1 - x0);          // length spans 1.0
    const cx = (x0 + x1) * 0.5, cz = (z0 + z1) * 0.5;
    const cy = ground ? y0 : (y0 + y1) * 0.5;       // wheels on the road
    for (let i = 0; i < P.length; i += 3) {
      P[i] = (P[i] - cx) * k;
      P[i + 1] = (P[i + 1] - cy) * k;
      P[i + 2] = (P[i + 2] - cz) * k;
    }
    return k;
  }

  /** Register an uploaded model as a selectable test body. */
  function addModel(key, name, P, opts) {
    const o = opts || {};
    orient(P, o.forward, o.up, o.flip, o.ground);
    const vox = Voxel.build(P, o.res || 96);

    // Measure the true frontal area from the voxels rather than assuming a
    // bounding-box fill factor. A thin-winged aircraft fills maybe 15% of its
    // box and a brick fills 100%; guessing wrong sizes the model badly and
    // throws the tunnel blockage off.
    const nx = vox.nx, ny = vox.ny, nz = vox.nz, g = vox.grid;
    let cols = 0;
    for (let z = 0; z < nz; z++) {
      for (let y = 0; y < ny; y++) {
        const base = nx * (y + ny * z);
        for (let x = 0; x < nx; x++) {
          if (g[base + x]) { cols++; break; }
        }
      }
    }
    const frontal = cols * vox.cell * vox.cell;
    const bbox = (2 * vox.ext[1]) * (2 * vox.ext[2]);
    const fill = bbox > 0 ? Math.max(0.05, Math.min(1, frontal / bbox)) : 0.75;

    BODIES[key] = {
      name, ground: !!o.ground, lengthM: o.lengthM || 4.5,
      refCd: 'unknown — measure it', custom: true,
      fill, parts: [],
      vox, tris: P
    };
    return BODIES[key];
  }

  // ---- a positioned, oriented, scaled body ---------------------------
  class Body {
    constructor(key) { this.setType(key); this.pos = [0, 0, 0]; this.scale = 30; this.pitch = 0; this.yaw = 0; }

    setType(key) {
      this.key = key;
      this.def = BODIES[key];
      if (this.def.custom) {                    // voxel-backed uploaded model
        const v = this.def.vox;
        // ground models were dropped onto y = 0, so their vertical budget is
        // the full height; free-flying ones are centred, so it is the half
        this.ext = [
          Math.max(v.ext[0], 1e-3),
          Math.max(this.def.ground ? v.ext[1] * 2 : v.ext[1], 1e-3),
          Math.max(v.ext[2], 1e-3)
        ];
        this.radius = Math.max(this.ext[0], this.ext[1], this.ext[2]);
        this.parts = [];
        return;
      }
      this.parts = this.def.parts.map(p => {
        const q = Object.assign({}, p);
        q.r = p.r || [0, 0, 0];
        q.m = eulerMat(q.r[0], q.r[1], q.r[2]);   // primitive rotation
        return q;
      });
      // axis-aligned half-extents in body space (used for a tight AABB)
      const e = [0, 0, 0];
      for (const p of this.parts) {
        const s = p.s;
        let h;
        switch (p.t) {
          case 'box': case 'ell': h = [s[0], s[1], s[2]]; break;
          case 'cyl': case 'cone': h = [s[0], s[1], s[0]]; break;
          case 'foil': h = [s[0] * 0.5, s[0] * 0.22, s[2] * 0.5]; break;
          default: h = [0, 0, 0];
        }
        const m = p.m;                       // primitive -> body rotation
        for (let i = 0; i < 3; i++) {
          const w = Math.abs(m[i * 3]) * h[0] + Math.abs(m[i * 3 + 1]) * h[1] + Math.abs(m[i * 3 + 2]) * h[2];
          e[i] = Math.max(e[i], Math.abs(p.c[i]) + w);
        }
      }
      this.ext = e;
      this.radius = Math.max(e[0], e[1], e[2]);
    }

    /**
     * Pick a size (lattice cells per unit length). Real tunnels keep the
     * model under ~10% of the test-section area because blockage inflates
     * Cd badly, so size the body by its estimated frontal area rather than
     * by bounding box alone — that lets a thin wing span much wider than a
     * blunt body without cheating on blockage.
     * Ground vehicles measure their height from y = 0, not the centre.
     */
    fit(nx, ny, nz) {
      const e = this.ext;
      const hy = this.def.ground ? e[1] : 2 * e[1];   // vertical size
      const wz = 2 * e[2];                            // spanwise size
      const fill = this.def.fill || 0.7;              // frontal area / bbox area
      const blockage = Math.sqrt(0.10 * ny * nz / (hy * wz * fill));
      this.scale = Math.min(
        blockage,
        nx * 0.34 / (2 * e[0]),   // leave room to develop a wake
        ny * 0.55 / hy,           // clear of the ceiling
        nz * 0.75 / wz);          // clear of the side walls
      return this.scale;
    }

    /** Recompute the body -> world rotation (yaw about y, then pitch about z). */
    _updateMat() {
      const cy = Math.cos(this.yaw), sy = Math.sin(this.yaw);
      const cp = Math.cos(this.pitch), sp = Math.sin(this.pitch);
      // R = Ry(yaw) * Rz(pitch)
      this.M = [
        cy * cp, -cy * sp, sy,
        sp, cp, 0,
        -sy * cp, sy * sp, cy
      ];
    }

    /** Point test in lattice coordinates. */
    makeTest() {
      this._updateMat();
      const M = this.M, P = this.pos, inv = 1 / this.scale;

      if (this.def.custom) {
        // Uploaded model: the mesh was voxelised once into body space, so
        // this is a transform plus one array read — the same cost as the
        // analytic bodies, and independent of triangle count.
        const v = this.def.vox;
        const g = v.grid, gnx = v.nx, gny = v.ny, gnz = v.nz;
        const gx0 = v.min[0], gy0 = v.min[1], gz0 = v.min[2], gi = v.inv;
        return function (X, Y, Z) {
          const dx = X - P[0], dy = Y - P[1], dz = Z - P[2];
          const lx = (M[0] * dx + M[3] * dy + M[6] * dz) * inv;
          const ly = (M[1] * dx + M[4] * dy + M[7] * dz) * inv;
          const lz = (M[2] * dx + M[5] * dy + M[8] * dz) * inv;
          const a = ((lx - gx0) * gi) | 0;
          if (a < 0 || a >= gnx) return false;
          const b = ((ly - gy0) * gi) | 0;
          if (b < 0 || b >= gny) return false;
          const c = ((lz - gz0) * gi) | 0;
          if (c < 0 || c >= gnz) return false;
          return g[a + gnx * (b + gny * c)] !== 0;
        };
      }

      const parts = this.parts;
      const l = [0, 0, 0], qv = [0, 0, 0];
      return function (X, Y, Z) {
        const dx = X - P[0], dy = Y - P[1], dz = Z - P[2];
        // world -> body : multiply by M^T then scale down
        l[0] = (M[0] * dx + M[3] * dy + M[6] * dz) * inv;
        l[1] = (M[1] * dx + M[4] * dy + M[7] * dz) * inv;
        l[2] = (M[2] * dx + M[5] * dy + M[8] * dz) * inv;
        for (let k = 0; k < parts.length; k++) {
          const p = parts[k], m = p.m, c = p.c;
          const ax = l[0] - c[0], ay = l[1] - c[1], az = l[2] - c[2];
          qv[0] = m[0] * ax + m[3] * ay + m[6] * az;
          qv[1] = m[1] * ax + m[4] * ay + m[7] * az;
          qv[2] = m[2] * ax + m[5] * ay + m[8] * az;
          if (insidePrim(p, qv)) return true;
        }
        return false;
      };
    }

    /** Does this body have wheels we can spin? */
    get hasWheels() { return this.parts.some(p => p.wheel); }

    /**
     * Surface velocity of the rotating wheels, in lattice units.
     *
     * In the tunnel frame the car stands still and the road sweeps past at
     * +u0. A wheel rolling without slip must match the road at its contact
     * patch, so omega = u0/R about the axle, which puts the contact point at
     * +u0 and the crown at -u0 — the crown meets the oncoming air at twice
     * the free-stream speed, which is exactly why a rolling road changes the
     * answer.
     *
     * @returns fn(X,Y,Z,out) -> true when the point is inside a wheel
     */
    makeWheelVel(u0) {
      this._updateMat();
      const M = this.M, P = this.pos, s = this.scale, inv = 1 / s;
      const wheels = [];
      for (const p of this.parts) {
        if (!p.wheel) continue;
        const m = p.m;
        // the cylinder's own axis is its local +y; take it to body space,
        // then to world. Normalise the sense so it always points along +z.
        let ab = [m[1], m[4], m[7]];
        if (ab[2] < 0) ab = [-ab[0], -ab[1], -ab[2]];
        const aw = [0, 0, 0], cw = [0, 0, 0];
        for (let i = 0; i < 3; i++) {
          aw[i] = M[i * 3] * ab[0] + M[i * 3 + 1] * ab[1] + M[i * 3 + 2] * ab[2];
          cw[i] = P[i] + s * (M[i * 3] * p.c[0] + M[i * 3 + 1] * p.c[1] + M[i * 3 + 2] * p.c[2]);
        }
        const R = p.s[0] * s;
        wheels.push({ p, m, c: p.c, axis: aw, centre: cw, omega: R > 0 ? u0 / R : 0 });
      }
      if (!wheels.length) return null;

      const l = [0, 0, 0], qv = [0, 0, 0];
      return function (X, Y, Z, out) {
        const dx = X - P[0], dy = Y - P[1], dz = Z - P[2];
        l[0] = (M[0] * dx + M[3] * dy + M[6] * dz) * inv;
        l[1] = (M[1] * dx + M[4] * dy + M[7] * dz) * inv;
        l[2] = (M[2] * dx + M[5] * dy + M[8] * dz) * inv;
        for (let k = 0; k < wheels.length; k++) {
          const w = wheels[k], m = w.m, c = w.c;
          const ax = l[0] - c[0], ay = l[1] - c[1], az = l[2] - c[2];
          qv[0] = m[0] * ax + m[3] * ay + m[6] * az;
          qv[1] = m[1] * ax + m[4] * ay + m[7] * az;
          qv[2] = m[2] * ax + m[5] * ay + m[8] * az;
          if (!insidePrim(w.p, qv)) continue;
          // u = omega * (axis x r), r measured from the wheel centre
          const rx = X - w.centre[0], ry = Y - w.centre[1], rz = Z - w.centre[2];
          const a = w.axis, o = w.omega;
          out[0] = o * (a[1] * rz - a[2] * ry);
          out[1] = o * (a[2] * rx - a[0] * rz);
          out[2] = o * (a[0] * ry - a[1] * rx);
          return true;
        }
        return false;
      };
    }

    /** Tight axis-aligned bounding box in lattice coordinates. */
    box() {
      this._updateMat();
      const M = this.M, e = this.ext, s = this.scale, P = this.pos;
      const w = [0, 0, 0];
      for (let i = 0; i < 3; i++) {
        w[i] = s * (Math.abs(M[i * 3]) * e[0] + Math.abs(M[i * 3 + 1]) * e[1] +
          Math.abs(M[i * 3 + 2]) * e[2]) + 2;
      }
      return {
        x0: P[0] - w[0], x1: P[0] + w[0],
        y0: P[1] - w[1], y1: P[1] + w[1],
        z0: P[2] - w[2], z1: P[2] + w[2]
      };
    }
  }

  // ---- three.js geometry for a primitive ------------------------------
  function primGeometry(p) {
    const s = p.s;
    switch (p.t) {
      case 'box': return new THREE.BoxGeometry(2 * s[0], 2 * s[1], 2 * s[2]);
      case 'ell': {
        const g = new THREE.SphereGeometry(1, 28, 18);
        g.scale(s[0], s[1], s[2]);
        return g;
      }
      case 'cyl': return new THREE.CylinderGeometry(s[0], s[0], 2 * s[1], 24);
      case 'cone': return new THREE.CylinderGeometry(0.0001, s[0], 2 * s[1], 24);
      case 'foil': {
        const ch = s[0], span = s[2];
        const m = p.naca ? p.naca[0] : 0, pp = p.naca ? p.naca[1] : 0.4, t = p.naca ? p.naca[2] : 0.12;
        const N = 40, pts = [];
        for (let i = 0; i <= N; i++) {                      // upper surface, LE -> TE
          const xc = 0.5 * (1 - Math.cos(Math.PI * i / N));
          pts.push(new THREE.Vector2((xc - 0.5) * ch, (nacaCamber(xc, m, pp) + nacaThickness(xc, t)) * ch));
        }
        for (let i = N - 1; i > 0; i--) {                   // lower surface, TE -> LE
          const xc = 0.5 * (1 - Math.cos(Math.PI * i / N));
          pts.push(new THREE.Vector2((xc - 0.5) * ch, (nacaCamber(xc, m, pp) - nacaThickness(xc, t)) * ch));
        }
        const shape = new THREE.Shape(pts);
        const g = new THREE.ExtrudeGeometry(shape, { depth: span, bevelEnabled: false, curveSegments: 4 });
        g.translate(0, 0, -span * 0.5);
        return g;
      }
    }
    return new THREE.BoxGeometry(1, 1, 1);
  }

  function buildMesh(body, material) {
    const group = new THREE.Group();
    if (body.def.custom) {
      const g = new THREE.BufferGeometry();
      g.setAttribute('position', new THREE.BufferAttribute(body.def.tris, 3));
      g.computeVertexNormals();
      group.add(new THREE.Mesh(g, material));
      return group;
    }
    for (const p of body.parts) {
      const mesh = new THREE.Mesh(primGeometry(p), material);
      mesh.position.set(p.c[0], p.c[1], p.c[2]);
      mesh.rotation.set(p.r[0], p.r[1], p.r[2], 'XYZ');
      group.add(mesh);
    }
    return group;
  }

  global.Shapes = {
    BODIES, Body, buildMesh, nacaThickness, nacaCamber,
    addModel, orient, guessAxes
  };
})(window);
