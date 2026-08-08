/* ------------------------------------------------------------------
   voxel.js — turn a triangle soup into a solid occupancy grid

   The solver needs to ask "is this point inside the body?" for every
   lattice cell, every time the body moves. Testing points against
   triangles directly is hopeless (a downloaded car is ~200k triangles),
   so we voxelise ONCE into a body-space grid and the per-cell test
   becomes an array lookup. Rotation and translation then cost nothing.

   Method: ray parity along all three axes, majority vote, OR'd with a
   rasterised surface shell.

   Parity alone is exact for a clean watertight mesh but degenerates on
   the self-intersecting, hole-ridden meshes people actually download.
   Voting across three independent axes recovers most of that, and the
   fraction of cells where the axes disagree is a genuine watertightness
   metric — we report it so the UI can warn instead of silently
   simulating a wrong shape.

   The result is deliberately CONSERVATIVE: every cell the surface passes
   through is marked solid, unconditionally. That costs about half a cell
   of outward bias, and buys the guarantee that features thinner than one
   cell — spoilers, mirrors, wings — still block flow instead of vanishing.
   For CFD that trade is the right way round.
------------------------------------------------------------------ */
(function (global) {
  'use strict';

  const BUCKETS = 48;          // 2D acceleration grid for the ray casts
  const MAX_SAMPLES = 4096;    // per-triangle cap for surface rasterisation

  function bounds(P) {
    let x0 = Infinity, y0 = Infinity, z0 = Infinity;
    let x1 = -Infinity, y1 = -Infinity, z1 = -Infinity;
    for (let i = 0; i < P.length; i += 3) {
      const x = P[i], y = P[i + 1], z = P[i + 2];
      if (x < x0) x0 = x; if (x > x1) x1 = x;
      if (y < y0) y0 = y; if (y > y1) y1 = y;
      if (z < z0) z0 = z; if (z > z1) z1 = z;
    }
    return [x0, y0, z0, x1, y1, z1];
  }

  /**
   * @param P    Float32Array, 9 floats per triangle, already in body space
   * @param res  target cells along the longest axis
   */
  function build(P, res) {
    res = res || 96;
    const ntri = (P.length / 9) | 0;
    if (!ntri) throw new Error('Mesh has no triangles.');

    const b = bounds(P);
    const sx = b[3] - b[0], sy = b[4] - b[1], sz = b[5] - b[2];
    const span = Math.max(sx, sy, sz);
    if (!(span > 0)) throw new Error('Mesh is degenerate (zero size).');
    const cell = span / res;
    const PAD = 2;

    const nx = Math.max(3, Math.ceil(sx / cell) + 2 * PAD);
    const ny = Math.max(3, Math.ceil(sy / cell) + 2 * PAD);
    const nz = Math.max(3, Math.ceil(sz / cell) + 2 * PAD);
    // centre the mesh bbox inside the padded grid
    const min = [
      b[0] - (nx * cell - sx) * 0.5,
      b[1] - (ny * cell - sy) * 0.5,
      b[2] - (nz * cell - sz) * 0.5
    ];
    const dim = [nx, ny, nz];
    const n = nx * ny * nz;
    if (n > 40e6) throw new Error('Voxel grid too large; lower the resolution.');

    const votes = new Uint8Array(n);       // 0..3 — how many axes say "inside"
    const shell = new Uint8Array(n);

    rasteriseSurface(P, ntri, shell, min, cell, dim);
    for (let axis = 0; axis < 3; axis++) parityFill(P, ntri, votes, min, cell, dim, axis);

    const grid = new Uint8Array(n);
    let solid = 0, disagree = 0;
    for (let i = 0; i < n; i++) {
      const v = votes[i];
      const inside = v >= 2 || shell[i] === 1;
      if (inside) solid++;
      if (v === 1 || v === 2) disagree++;
      grid[i] = inside ? 1 : 0;
    }

    return {
      grid, nx, ny, nz, cell,
      min,
      inv: 1 / cell,
      ext: [sx * 0.5, sy * 0.5, sz * 0.5],
      centre: [(b[0] + b[3]) * 0.5, (b[1] + b[4]) * 0.5, (b[2] + b[5]) * 0.5],
      bbox: b,
      stats: {
        tris: ntri, solid,
        // how unreliable the parity vote was — a good watertightness proxy
        openness: solid ? disagree / solid : 0
      }
    };
  }

  /** Mark every cell the triangle surface passes through. */
  function rasteriseSurface(P, ntri, shell, min, cell, dim) {
    const [nx, ny, nz] = dim;
    const inv = 1 / cell;
    for (let t = 0; t < ntri; t++) {
      const o = t * 9;
      const ax = P[o], ay = P[o + 1], az = P[o + 2];
      const bx = P[o + 3], by = P[o + 4], bz = P[o + 5];
      const cx = P[o + 6], cy = P[o + 7], cz = P[o + 8];
      // sample density from the two longest edges
      const e1 = Math.max(Math.abs(bx - ax), Math.abs(by - ay), Math.abs(bz - az));
      const e2 = Math.max(Math.abs(cx - ax), Math.abs(cy - ay), Math.abs(cz - az));
      let n1 = Math.ceil(e1 * inv * 2) + 1, n2 = Math.ceil(e2 * inv * 2) + 1;
      if (n1 * n2 > MAX_SAMPLES) {
        const k = Math.sqrt(MAX_SAMPLES / (n1 * n2));
        n1 = Math.max(2, (n1 * k) | 0); n2 = Math.max(2, (n2 * k) | 0);
      }
      for (let i = 0; i <= n1; i++) {
        const u = i / n1;
        const jmax = Math.ceil((1 - u) * n2);
        for (let j = 0; j <= jmax; j++) {
          const v = Math.min(j / n2, 1 - u);
          const x = ax + (bx - ax) * u + (cx - ax) * v;
          const y = ay + (by - ay) * u + (cy - ay) * v;
          const z = az + (bz - az) * u + (cz - az) * v;
          const gx = ((x - min[0]) * inv) | 0;
          const gy = ((y - min[1]) * inv) | 0;
          const gz = ((z - min[2]) * inv) | 0;
          if (gx < 0 || gy < 0 || gz < 0 || gx >= nx || gy >= ny || gz >= nz) continue;
          shell[gx + nx * (gy + ny * gz)] = 1;
        }
      }
    }
  }

  /**
   * Cast one ray per cell column along `axis`, sort the crossings and fill
   * between alternate pairs. Triangles are bucketed in the perpendicular
   * plane so each ray only tests a handful of them.
   */
  function parityFill(P, ntri, votes, min, cell, dim, axis) {
    const A = axis, B = (axis + 1) % 3, C = (axis + 2) % 3;
    const nA = dim[A], nB = dim[B], nC = dim[C];
    const minA = min[A], minB = min[B], minC = min[C];
    const spanB = nB * cell, spanC = nC * cell;

    // ---- bucket triangles by their (B,C) footprint (CSR) ----
    const nbB = Math.min(BUCKETS, nB), nbC = Math.min(BUCKETS, nC);
    const kB = nbB / spanB, kC = nbC / spanC;
    const nb = nbB * nbC;
    const count = new Int32Array(nb + 1);
    const lo = new Int32Array(ntri * 2), hi = new Int32Array(ntri * 2);

    for (let t = 0; t < ntri; t++) {
      const o = t * 9;
      const b0 = P[o + B], b1 = P[o + 3 + B], b2 = P[o + 6 + B];
      const c0 = P[o + C], c1 = P[o + 3 + C], c2 = P[o + 6 + C];
      let bl = Math.floor((Math.min(b0, b1, b2) - minB) * kB);
      let bh = Math.floor((Math.max(b0, b1, b2) - minB) * kB);
      let cl = Math.floor((Math.min(c0, c1, c2) - minC) * kC);
      let ch = Math.floor((Math.max(c0, c1, c2) - minC) * kC);
      bl = bl < 0 ? 0 : bl; bh = bh >= nbB ? nbB - 1 : bh;
      cl = cl < 0 ? 0 : cl; ch = ch >= nbC ? nbC - 1 : ch;
      lo[t * 2] = bl; hi[t * 2] = bh; lo[t * 2 + 1] = cl; hi[t * 2 + 1] = ch;
      if (bl > bh || cl > ch) continue;
      for (let cc = cl; cc <= ch; cc++) {
        for (let bb = bl; bb <= bh; bb++) count[bb + nbB * cc + 1]++;
      }
    }
    for (let i = 0; i < nb; i++) count[i + 1] += count[i];
    const items = new Int32Array(count[nb]);
    const cursor = count.slice(0, nb);
    for (let t = 0; t < ntri; t++) {
      const bl = lo[t * 2], bh = hi[t * 2], cl = lo[t * 2 + 1], ch = hi[t * 2 + 1];
      if (bl > bh || cl > ch) continue;
      for (let cc = cl; cc <= ch; cc++) {
        for (let bb = bl; bb <= bh; bb++) items[cursor[bb + nbB * cc]++] = t;
      }
    }

    // ---- one ray per column ----
    const strideA = A === 0 ? 1 : (A === 1 ? dim[0] : dim[0] * dim[1]);
    const strideB = B === 0 ? 1 : (B === 1 ? dim[0] : dim[0] * dim[1]);
    const strideC = C === 0 ? 1 : (C === 1 ? dim[0] : dim[0] * dim[1]);
    const hits = new Float64Array(512);

    for (let ic = 0; ic < nC; ic++) {
      // irrational-ish offset keeps rays off shared triangle edges
      const pc = minC + (ic + 0.5003719) * cell;
      const bc = Math.min(nbC - 1, Math.max(0, Math.floor((pc - minC) * kC)));
      for (let ib = 0; ib < nB; ib++) {
        const pb = minB + (ib + 0.4996281) * cell;
        const bb = Math.min(nbB - 1, Math.max(0, Math.floor((pb - minB) * kB)));
        const s = count[bb + nbB * bc], e = count[bb + nbB * bc + 1];
        if (s === e) continue;

        let nh = 0;
        for (let k = s; k < e && nh < hits.length; k++) {
          const o = items[k] * 9;
          const b0 = P[o + B], c0 = P[o + C];
          const b1 = P[o + 3 + B], c1 = P[o + 3 + C];
          const b2 = P[o + 6 + B], c2 = P[o + 6 + C];
          // barycentric test in the (B,C) projection
          const d = (c1 - c2) * (b0 - b2) + (b2 - b1) * (c0 - c2);
          if (d === 0) continue;
          const id = 1 / d;
          const l0 = ((c1 - c2) * (pb - b2) + (b2 - b1) * (pc - c2)) * id;
          if (l0 < 0 || l0 > 1) continue;
          const l1 = ((c2 - c0) * (pb - b2) + (b0 - b2) * (pc - c2)) * id;
          if (l1 < 0 || l1 > 1) continue;
          const l2 = 1 - l0 - l1;
          if (l2 < 0 || l2 > 1) continue;
          hits[nh++] = l0 * P[o + A] + l1 * P[o + 3 + A] + l2 * P[o + 6 + A];
        }
        if (nh < 2) continue;

        const h = hits.subarray(0, nh);
        h.sort();                                  // typed arrays sort numerically
        const colBase = ib * strideB + ic * strideC;
        for (let k = 0; k + 1 < nh; k += 2) {
          let i0 = Math.ceil((h[k] - minA) / cell - 0.5);
          let i1 = Math.floor((h[k + 1] - minA) / cell - 0.5);
          if (i0 < 0) i0 = 0;
          if (i1 >= nA) i1 = nA - 1;
          for (let ia = i0; ia <= i1; ia++) votes[colBase + ia * strideA]++;
        }
      }
    }
  }

  global.Voxel = { build, bounds };
})(window);
