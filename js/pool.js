/* ------------------------------------------------------------------
   pool.js — run the solver across every available core

   Domain decomposition along z. The lattice lives in SharedArrayBuffers so
   the workers operate on it in place: no copying, no halo messages. A batch
   of steps costs exactly two postMessages per worker, everything in between
   is Atomics.

   SharedArrayBuffer is gated behind cross-origin isolation, so this is only
   available when the page is served with COOP/COEP headers (serve.py does).
   Opened straight off disk as a file:// URL, it cannot run and the app stays
   single-core.
------------------------------------------------------------------ */
(function (global) {
  'use strict';

  const C_PARITY = 0, C_ARRIVE = 1, C_GEN = 2;

  /** Why multi-core is unavailable, or null when it is available. */
  function unavailableReason() {
    if (typeof SharedArrayBuffer === 'undefined') {
      return location.protocol === 'file:'
        ? 'needs a server — open via serve.py, not from disk'
        : 'SharedArrayBuffer unavailable in this browser';
    }
    if (!global.crossOriginIsolated) {
      return location.protocol === 'file:'
        ? 'needs a server — open via serve.py, not from disk'
        : 'page is not cross-origin isolated (COOP/COEP headers missing)';
    }
    if (typeof Worker === 'undefined') return 'Web Workers unavailable';
    return null;
  }

  /**
   * Default worker count.
   *
   * Deliberately about half the logical cores rather than all-but-one.
   * hardwareConcurrency counts SMT siblings, and this kernel is memory
   * bound, so pairs of threads on one physical core mostly fight over the
   * same load/store ports. Oversubscribing measurably *lost* to single core
   * in testing. The Benchmark button tunes it properly per machine.
   */
  function suggestedWorkers() {
    const hc = global.navigator && navigator.hardwareConcurrency || 4;
    return Math.max(2, Math.min(8, Math.floor(hc / 2)));
  }

  class Pool {
    constructor(sim, nWorkers) {
      this.sim = sim;
      this.n = nWorkers;
      this.workers = [];
      this.busy = false;
      this.maxBatch = 16;
      // per-slab force for every (step, worker); the pool totals them after
      // the batch, so no worker has to act as a serial leader
      this.forces = new Float64Array(
        new SharedArrayBuffer(this.maxBatch * nWorkers * 3 * 8));
      this.stepForces = [];
    }

    async start() {
      const s = this.sim;
      const buffers = s.buffers();
      Atomics.store(new Int32Array(buffers.ctrl), C_ARRIVE, 0);
      const ready = [];
      for (let i = 0; i < this.n; i++) {
        const w = new Worker('js/lbmWorker.js');
        this.workers.push(w);
        ready.push(new Promise((res, rej) => {
          w.onmessage = (e) => { if (e.data.ok === 'init') res(e.data); };
          w.onerror = (e) => rej(new Error(e.message || 'worker failed'));
        }));
        w.postMessage({
          cmd: 'init', id: i, nWorkers: this.n,
          nx: s.nx, ny: s.ny, nz: s.nz,
          buffers, forces: this.forces.buffer
        });
      }
      this.slabs = await Promise.all(ready);
      return this;
    }

    /** Push solver settings the workers keep their own copy of. */
    syncParams() {
      const s = this.sim;
      for (const w of this.workers) {
        w.postMessage({
          cmd: 'params', tau: s.tau, les: s.les, u0: s.u0,
          movingWalls: s.movingWalls
        });
      }
    }

    /** Advance n steps across all cores. Resolves when the batch is done. */
    step(n) {
      n = Math.max(1, Math.min(this.maxBatch, n | 0));
      this.busy = true;
      return new Promise((resolve, reject) => {
        let left = this.n;
        const parity = this.sim.parity;
        for (const w of this.workers) {
          w.onmessage = () => { if (--left === 0) finish(); };
          w.onerror = (e) => { this.busy = false; reject(e); };
          w.postMessage({ cmd: 'step', n, parity });
        }
        const finish = () => {
          // sum each step's slab contributions into a per-step total
          this.stepForces.length = 0;
          for (let s = 0; s < n; s++) {
            let fx = 0, fy = 0, fz = 0;
            for (let w = 0; w < this.n; w++) {
              const o = (s * this.n + w) * 3;
              fx += this.forces[o]; fy += this.forces[o + 1]; fz += this.forces[o + 2];
            }
            this.stepForces.push([fx, fy, fz]);
          }
          this.sim.parity ^= (n & 1);        // keep the page's view in step
          this.sim.steps += n;
          this.busy = false;
          resolve(this.stepForces);
        };
      });
    }

    dispose() {
      for (const w of this.workers) w.terminate();
      this.workers.length = 0;
    }
  }

  global.LBMPool = { Pool, unavailableReason, suggestedWorkers };
})(window);
