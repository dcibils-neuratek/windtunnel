/* ------------------------------------------------------------------
   lbmWorker.js — one slab of the lattice, on one core

   Each worker rebuilds an LBM object over the SAME SharedArrayBuffers the
   main thread allocated, then owns a z-slab of the domain. Because the
   solver is double buffered, every worker reads f and writes a disjoint
   part of g, so there is no locking and no halo exchange — only a barrier
   before the buffers are swapped.

   The kernel itself is not duplicated here: this runs the very same
   stepSlab() that the single-threaded path runs, so the two can never
   drift apart.
------------------------------------------------------------------ */
'use strict';
importScripts('lbm.js');

let sim = null, id = 0, nWorkers = 1, z0 = 1, z1 = 1, bz0 = 0, bz1 = 1;
let ctrl = null, forces = null;

// control slots shared with the pool (Int32Array)
const C_ARRIVE = 1, C_GEN = 2;

/** Sense-reversing barrier. */
function barrier() {
  const gen = Atomics.load(ctrl, C_GEN);
  if (Atomics.add(ctrl, C_ARRIVE, 1) === nWorkers - 1) {
    Atomics.store(ctrl, C_ARRIVE, 0);
    Atomics.add(ctrl, C_GEN, 1);
    Atomics.notify(ctrl, C_GEN);
  } else {
    while (Atomics.load(ctrl, C_GEN) === gen) Atomics.wait(ctrl, C_GEN, gen);
  }
}

self.onmessage = (e) => {
  const m = e.data;

  if (m.cmd === 'init') {
    id = m.id; nWorkers = m.nWorkers;
    sim = new LBM(m.nx, m.ny, m.nz, { buffers: m.buffers });
    ctrl = new Int32Array(m.buffers.ctrl);
    forces = new Float64Array(m.forces);
    // split the interior [1, nz-1) as evenly as the slabs allow
    const lo = 1, hi = m.nz - 1, span = hi - lo;
    z0 = lo + Math.floor(span * id / nWorkers);
    z1 = lo + Math.floor(span * (id + 1) / nWorkers);
    bz0 = id === 0 ? 0 : z0;
    bz1 = id === nWorkers - 1 ? m.nz : z1;
    self.postMessage({ ok: 'init', id, z0, z1 });
    return;
  }

  if (m.cmd === 'params') {             // tau / LES / moving walls changed
    sim.tau = m.tau; sim.les = m.les; sim.u0 = m.u0;
    sim.movingWalls = m.movingWalls;
    return;
  }

  if (m.cmd === 'step') {
    sim.parity = m.parity;
    for (let s = 0; s < m.n; s++) {
      // Interior slab, then this slab's share of the domain boundaries. The
      // boundary pass writes g and reads only f, so it needs no barrier of
      // its own — which keeps the whole step down to a single barrier and
      // leaves no serial section for one thread to run alone.
      sim.stepSlab(z0, z1);
      // Boundary ranges must tile [0,nz) exactly, so the end workers extend
      // over the spanwise end planes; otherwise the outlet cells at z=0 and
      // z=nz-1 would keep their inlet value instead of the zero-gradient copy.
      sim._boundaries(sim.g, bz0, bz1, id === 0, id === nWorkers - 1);

      const b = (s * nWorkers + id) * 3;          // pool totals these up
      forces[b] = sim._sf[0]; forces[b + 1] = sim._sf[1]; forces[b + 2] = sim._sf[2];

      barrier();                                  // all of g is written
      sim.parity ^= 1;                            // every worker agrees
    }
    self.postMessage({ ok: 'step', id });
    return;
  }
};
