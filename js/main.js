/* ------------------------------------------------------------------
   main.js — scene, visualisation, UI and interaction
------------------------------------------------------------------ */
(function () {
  'use strict';

  const RES = {
    fast: [72, 34, 34],
    balanced: [96, 44, 44],
    detail: [124, 56, 56]
  };
  const RHO_AIR = 1.225;      // kg/m^3
  const NU_AIR = 1.48e-5;     // m^2/s

  const S = {                 // UI state
    res: 'balanced',
    compute: 'single',        // 'single' | 'multi' (worker pool)
    workers: 0,               // 0 = auto; the Benchmark button tunes it
    shape: 'sedan',
    windMs: 30,
    lengthM: 4.6,
    aoa: 0,
    yaw: 0,
    ground: true,
    les: true,
    reKnob: 0.62,             // 0..1 -> relaxation time
    rollingRoad: false,       // moving ground + spinning wheels; opt-in,
                              // because it needs a resolved ride height
    running: true,
    smoke: true,
    ribbons: true,            // smoke-wand streaklines
    nParticles: 9000,
    flowSpeed: 14,
    slice: true,
    sliceAxis: 'xy',
    sliceField: 'speed',
    streamlines: false,
    showBody: true,
    units: localStorage.getItem('wt-units') === 'imperial' ? 'imperial' : 'metric',
    theme: localStorage.getItem('wt-theme') === 'irix' ? 'irix' : 'modern'
  };

  // scene colours per look-and-feel; the CSS half lives in style.css
  const THEME = {
    modern: { bg: 0x070a10, tunnel: 0x2b3d55, road: 0x151b24, line: 0x8fe6ff },
    irix: { bg: 0x191934, tunnel: 0x7d76b4, road: 0x2b2b52, line: 0xa9ecff }
  };

  let sim, body;
  let scene, camera, renderer, orbit, raycaster;
  let bodyGroup, tunnelBox, roadMesh, sliceMesh, sliceTex, sliceData;
  let points, pPos, pCol, pAge, pGeo;
  let lineGeo, lineMesh, linePos;
  let needsVoxelize = true, dragging = null;
  let pool = null, poolReady = null, batchInFlight = false;
  const tmpV = [0, 0, 0];
  const forceEMA = [0, 0, 0];
  let fCount = 0, settle = 0;

  // ---- measurement ------------------------------------------------------
  // Coefficient history, so a run can be judged and exported rather than just
  // glanced at. Ring buffer of post-settle samples, one per lattice step.
  const HIST = 6000;
  const hCd = new Float32Array(HIST), hCl = new Float32Array(HIST), hCy = new Float32Array(HIST);
  const scratch = new Float64Array(HIST);
  let hN = 0;                 // total samples taken (may exceed HIST)

  /** Forces are meaningless until the wake has developed; restart averaging. */
  function resetAverages(steps) {
    forceEMA[0] = forceEMA[1] = forceEMA[2] = 0;
    fCount = 0; hN = 0; havePrev = false;
    settle = steps === undefined ? 150 : steps;
  }

  /** Reference area in cells — planform for wings, frontal for everything else. */
  function refCells() {
    return body.def.ref === 'planform' ? sim.planformArea : sim.frontalArea;
  }

  const prevF = new Float64Array(3);
  let havePrev = false;

  /**
   * Record one sample of the force coefficients.
   *
   * The momentum-exchange force carries a strong odd/even parity oscillation —
   * measured lag-1 autocorrelation of -0.98, lag-2 of +0.98 — which is a
   * well-known lattice Boltzmann artefact: the mean is correct but the value
   * alternates every timestep. Averaging consecutive pairs cancels it exactly,
   * leaving the mean untouched while cutting the scatter ~16x. Without this
   * the trace is an unreadable band and the error bar measures the artefact
   * rather than the flow.
   */
  /** Fold one completed timestep's force into the running statistics. */
  function accumulate() {
    if (settle > 0) { settle--; return; }
    const a = 1 / Math.min(++fCount, 500);
    for (let k = 0; k < 3; k++) forceEMA[k] += (sim.force[k] - forceEMA[k]) * a;
    pushSample();
  }

  function pushSample() {
    if (!havePrev) { prevF.set(sim.force); havePrev = true; return; }
    const q = 0.5 * sim.u0 * sim.u0 * refCells();
    const i = hN % HIST;
    hCd[i] = 0.5 * (sim.force[0] + prevF[0]) / q;
    hCl[i] = 0.5 * (sim.force[1] + prevF[1]) / q;
    hCy[i] = 0.5 * (sim.force[2] + prevF[2]) / q;
    prevF.set(sim.force);
    hN++;
  }

  /** Copy the ring into chronological order; returns the sample count. */
  function ordered(src, out) {
    const n = Math.min(hN, HIST);
    const s = hN <= HIST ? 0 : hN % HIST;
    for (let i = 0; i < n; i++) out[i] = src[(s + i) % HIST];
    return n;
  }

  /**
   * Mean and standard error by the blocking method.
   *
   * Consecutive timesteps in a shedding wake are heavily correlated, so the
   * textbook sigma/sqrt(N) understates the true error by a large factor.
   * Averaging into blocks longer than the correlation time decorrelates them,
   * and the scatter *between block means* is an honest error estimate.
   */
  function blockStats(x, n) {
    const B = 16;
    const bl = Math.floor(n / B);
    if (bl < 4) {                                  // too short to block
      let m = 0;
      for (let i = 0; i < n; i++) m += x[i];
      m /= n;
      let v = 0;
      for (let i = 0; i < n; i++) { const d = x[i] - m; v += d * d; }
      return { mean: m, sem: Math.sqrt(v / Math.max(1, n - 1) / n), blocked: false };
    }
    const used = B * bl;
    const bm = new Float64Array(B);
    let mean = 0;
    for (let b = 0; b < B; b++) {
      let s = 0;
      for (let i = b * bl; i < (b + 1) * bl; i++) s += x[i];
      bm[b] = s / bl; mean += bm[b];
    }
    mean /= B;
    let v = 0;
    for (let b = 0; b < B; b++) { const d = bm[b] - mean; v += d * d; }
    return { mean, sem: Math.sqrt(v / (B - 1) / B), blocked: true, used };
  }

  /**
   * Full statistics for the current run, or null while there is too little
   * data. `state` is what a technician actually needs to know: is this number
   * finished moving?
   */
  function runStats() {
    const n = ordered(hCd, scratch);
    if (n < 60) return null;
    const cd = blockStats(scratch, n);

    // still drifting? compare the two halves against the error bar
    let m1 = 0, m2 = 0;
    const h = n >> 1;
    for (let i = 0; i < h; i++) m1 += scratch[i];
    for (let i = h; i < n; i++) m2 += scratch[i];
    m1 /= h; m2 /= (n - h);
    const drift = Math.abs(m1 - m2);

    const nl = ordered(hCl, scratch); const cl = blockStats(scratch, nl);
    const ny = ordered(hCy, scratch); const cy = blockStats(scratch, ny);

    const rel = 1.96 * cd.sem / Math.max(1e-9, Math.abs(cd.mean));
    const drifting = drift > 2.5 * cd.sem;
    const state = n < 400 ? 'averaging' : drifting ? 'drifting' : rel < 0.03 ? 'converged' : 'noisy';
    return { n, cd, cl, cy, drift, rel, state };
  }
  let stepsPerFrame = 1, msPerStep = 8, fpsEMA = 60, lastT = performance.now();

  // ---- colour maps ---------------------------------------------------
  const CMAP = [
    [0.00, 0.02, 0.05, 0.18],
    [0.22, 0.06, 0.42, 0.82],
    [0.46, 0.10, 0.82, 0.72],
    [0.70, 0.92, 0.86, 0.24],
    [1.00, 1.00, 0.24, 0.16]
  ];
  function cmap(t, out) {
    t = t < 0 ? 0 : t > 1 ? 1 : t;
    for (let i = 1; i < CMAP.length; i++) {
      if (t <= CMAP[i][0]) {
        const a = CMAP[i - 1], b = CMAP[i];
        const k = (t - a[0]) / (b[0] - a[0]);
        out[0] = a[1] + (b[1] - a[1]) * k;
        out[1] = a[2] + (b[2] - a[2]) * k;
        out[2] = a[3] + (b[3] - a[3]) * k;
        return out;
      }
    }
    out[0] = 1; out[1] = 0.24; out[2] = 0.16; return out;
  }
  // Slice scaling — shared by the renderer and the legend so the colour key
  // can never drift from what is actually drawn.
  const K_SPEED = 0.62;        // t = |u|/u0 * K   -> bar tops out at 1/K
  const K_CP = 1.1;            // diverging(Cp * K)
  const K_VORT = 1 / 0.28;     // diverging(omega/u0 * K)

  function diverging(t, out) {         // -1 .. +1  (blue - dark - red)
    const a = Math.min(1, Math.abs(t));
    if (t < 0) { out[0] = 0.15 * a; out[1] = 0.45 * a + 0.04; out[2] = 1.0 * a + 0.06; }
    else { out[0] = 1.0 * a + 0.06; out[1] = 0.35 * a + 0.04; out[2] = 0.12 * a; }
    return out;
  }

  // ---- minimal orbit controls ----------------------------------------
  class Orbit {
    constructor(cam, dom) {
      this.cam = cam; this.dom = dom; this.enabled = true;
      this.target = new THREE.Vector3();
      this.theta = 0.9; this.phi = 1.15; this.radius = 100;
      this.min = 8; this.max = 900;
      this._drag = 0; this._px = 0; this._py = 0;
      dom.addEventListener('pointerdown', e => {
        if (!this.enabled) return;
        this._drag = (e.button === 0 && !e.shiftKey) ? 1 : 2;
        this._px = e.clientX; this._py = e.clientY;
        dom.setPointerCapture(e.pointerId);
      });
      dom.addEventListener('pointermove', e => {
        if (!this._drag || !this.enabled) return;
        const dx = e.clientX - this._px, dy = e.clientY - this._py;
        this._px = e.clientX; this._py = e.clientY;
        if (this._drag === 1) {
          this.theta -= dx * 0.006;
          this.phi = Math.max(0.05, Math.min(Math.PI - 0.05, this.phi - dy * 0.006));
        } else {
          const k = this.radius * 0.0016;
          const right = new THREE.Vector3().setFromMatrixColumn(this.cam.matrix, 0);
          const up = new THREE.Vector3().setFromMatrixColumn(this.cam.matrix, 1);
          this.target.addScaledVector(right, -dx * k).addScaledVector(up, dy * k);
        }
        this.update();
      });
      const end = () => { this._drag = 0; };
      dom.addEventListener('pointerup', end);
      dom.addEventListener('pointercancel', end);
      dom.addEventListener('wheel', e => {
        e.preventDefault();
        this.radius = Math.max(this.min, Math.min(this.max, this.radius * Math.pow(1.0016, e.deltaY)));
        this.update();
      }, { passive: false });
    }
    update() {
      const sp = Math.sin(this.phi), r = this.radius;
      this.cam.position.set(
        this.target.x + r * sp * Math.sin(this.theta),
        this.target.y + r * Math.cos(this.phi),
        this.target.z + r * sp * Math.cos(this.theta));
      this.cam.lookAt(this.target);
      this.cam.updateMatrixWorld();
    }
  }

  // ---- scene ----------------------------------------------------------
  function initScene() {
    const canvas = document.getElementById('view');
    renderer = new THREE.WebGLRenderer({ canvas, antialias: true });
    renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
    renderer.outputEncoding = THREE.sRGBEncoding;

    scene = new THREE.Scene();
    scene.background = new THREE.Color(0x070a10);
    scene.fog = new THREE.FogExp2(0x070a10, 0.0022);

    camera = new THREE.PerspectiveCamera(48, 1, 0.5, 4000);
    orbit = new Orbit(camera, canvas);
    raycaster = new THREE.Raycaster();

    scene.add(new THREE.HemisphereLight(0x9fc4ff, 0x101418, 0.85));
    const key = new THREE.DirectionalLight(0xffffff, 1.15);
    key.position.set(-1, 2, 1.5); scene.add(key);
    const rim = new THREE.DirectionalLight(0x66aaff, 0.5);
    rim.position.set(2, 0.6, -1.4); scene.add(rim);

    window.addEventListener('resize', resize);
    canvas.addEventListener('pointerdown', onPointerDown);
    window.addEventListener('pointermove', onPointerMove);
    window.addEventListener('pointerup', () => { dragging = null; orbit.enabled = true; });
    resize();
  }

  function resize() {
    const w = window.innerWidth, h = window.innerHeight;
    renderer.setSize(w, h, false);
    camera.aspect = w / h; camera.updateProjectionMatrix();
  }

  // ---- build everything for the current resolution ---------------------
  function buildWorld() {
    const [nx, ny, nz] = RES[S.res];
    if (pool) { pool.dispose(); pool = null; }
    // allocate shared whenever the platform allows it, so the multi-core
    // toggle never needs a rebuild
    sim = new LBM(nx, ny, nz, { shared: !LBMPool.unavailableReason() });
    poolReady = null;
    applySolverParams();

    if (tunnelBox) scene.remove(tunnelBox);
    const eg = new THREE.EdgesGeometry(new THREE.BoxGeometry(nx, ny, nz));
    tunnelBox = new THREE.LineSegments(eg, new THREE.LineBasicMaterial({ color: 0x2b3d55 }));
    tunnelBox.position.set(nx / 2, ny / 2, nz / 2);
    scene.add(tunnelBox);

    if (roadMesh) scene.remove(roadMesh);
    roadMesh = new THREE.Mesh(
      new THREE.PlaneGeometry(nx, nz),
      new THREE.MeshStandardMaterial({ color: 0x151b24, roughness: 0.95, metalness: 0 }));
    roadMesh.rotation.x = -Math.PI / 2;
    roadMesh.position.set(nx / 2, 1, nz / 2);
    scene.add(roadMesh);

    buildBody(true);            // before the particles — seeding needs the body
    buildParticles(nx, ny, nz);
    buildRibbons();
    buildSlice();
    buildStreamlines();

    orbit.target.set(nx * 0.46, ny * 0.34, nz / 2);
    orbit.radius = nx * 1.35;
    orbit.update();

    points.visible = S.smoke;
    ribbons.visible = S.ribbons;
    sliceMesh.visible = S.slice;
    lineMesh.visible = S.streamlines;
    applyThemeColors();         // the meshes above were just recreated
  }

  function buildParticles(nx, ny, nz) {
    if (points) { scene.remove(points); points.geometry.dispose(); }
    const P = S.nParticles;
    pPos = new Float32Array(P * 3);
    pCol = new Float32Array(P * 3);
    pAge = new Float32Array(P);
    pGeo = new THREE.BufferGeometry();
    pGeo.setAttribute('position', new THREE.BufferAttribute(pPos, 3));
    pGeo.setAttribute('color', new THREE.BufferAttribute(pCol, 3));

    // soft round sprite
    const c = document.createElement('canvas'); c.width = c.height = 64;
    const g = c.getContext('2d');
    const grd = g.createRadialGradient(32, 32, 0, 32, 32, 32);
    grd.addColorStop(0, 'rgba(255,255,255,1)');
    grd.addColorStop(0.35, 'rgba(255,255,255,0.55)');
    grd.addColorStop(1, 'rgba(255,255,255,0)');
    g.fillStyle = grd; g.fillRect(0, 0, 64, 64);
    const tex = new THREE.CanvasTexture(c);

    points = new THREE.Points(pGeo, new THREE.PointsMaterial({
      size: Math.max(0.7, nx / 90), map: tex, vertexColors: true,
      transparent: true, depthWrite: false, blending: THREE.AdditiveBlending, opacity: 0.9
    }));
    scene.add(points);
    for (let i = 0; i < P; i++) respawn(i, true);
  }

  function buildSlice() {
    if (sliceMesh) { scene.remove(sliceMesh); sliceMesh.geometry.dispose(); sliceTex.dispose(); }
    const { nx, ny, nz } = sim;
    const w = nx, h = S.sliceAxis === 'xy' ? ny : nz;
    sliceData = new Uint8Array(w * h * 4);
    sliceTex = new THREE.DataTexture(sliceData, w, h, THREE.RGBAFormat);
    sliceTex.minFilter = sliceTex.magFilter = THREE.LinearFilter;
    sliceTex.needsUpdate = true;
    const geo = new THREE.PlaneGeometry(nx, h);
    sliceMesh = new THREE.Mesh(geo, new THREE.MeshBasicMaterial({
      map: sliceTex, transparent: true, opacity: 0.72, depthWrite: false, side: THREE.DoubleSide
    }));
    if (S.sliceAxis === 'xz') sliceMesh.rotation.x = -Math.PI / 2;
    scene.add(sliceMesh);
    positionSlice();
  }

  function positionSlice() {
    if (!sliceMesh || !body) return;
    const { nx, ny, nz } = sim;
    if (S.sliceAxis === 'xy') sliceMesh.position.set(nx / 2, ny / 2, Math.round(body.pos[2]));
    else sliceMesh.position.set(nx / 2, Math.round(body.pos[1]), nz / 2);
  }

  function buildStreamlines() {
    if (lineMesh) { scene.remove(lineMesh); lineMesh.geometry.dispose(); }
    const NL = 220, LEN = 90;
    linePos = new Float32Array(NL * LEN * 2 * 3);
    lineGeo = new THREE.BufferGeometry();
    lineGeo.setAttribute('position', new THREE.BufferAttribute(linePos, 3));
    lineMesh = new THREE.LineSegments(lineGeo, new THREE.LineBasicMaterial({
      color: 0x8fe6ff, transparent: true, opacity: 0.38, depthWrite: false
    }));
    lineMesh.frustumCulled = false;
    lineMesh._nl = NL; lineMesh._len = LEN;
    scene.add(lineMesh);
  }

  function buildBody(recentre) {
    const { nx, ny, nz } = sim;
    const keep = body ? { p: body.pos.slice(), pitch: body.pitch, yaw: body.yaw } : null;
    body = new Shapes.Body(S.shape);
    body.fit(nx, ny, nz);
    const def = body.def;
    S.ground = def.ground;
    S.lengthM = def.lengthM;

    if (recentre || !keep) {
      body.pos = [nx * 0.36, def.ground ? 1 : ny * 0.5, nz * 0.5];
      body.pitch = 0; body.yaw = 0;
    } else {
      body.pos = keep.p; body.pitch = keep.pitch; body.yaw = keep.yaw;
      if (def.ground) body.pos[1] = 1;
    }
    // the nose points at -x, so a nose-up angle of attack is a NEGATIVE
    // rotation about +z
    body.pitch = -S.aoa * Math.PI / 180;
    body.yaw = S.yaw * Math.PI / 180;

    if (bodyGroup) { scene.remove(bodyGroup); }
    const mat = new THREE.MeshStandardMaterial({
      color: 0xdbe4ef, metalness: 0.45, roughness: 0.32, envMapIntensity: 1
    });
    bodyGroup = Shapes.buildMesh(body, mat);
    bodyGroup.rotation.order = 'YZX';
    scene.add(bodyGroup);

    roadMesh.visible = def.ground;
    document.getElementById('refcd').textContent = def.refCd;
    document.getElementById('lengthM').value = def.lengthM;
    document.getElementById('lengthMv').textContent = fmtLength(def.lengthM);
    syncBodyTransform();
  }

  function syncBodyTransform() {
    bodyGroup.position.set(body.pos[0], body.pos[1], body.pos[2]);
    bodyGroup.rotation.set(0, body.yaw, body.pitch, 'YZX');
    bodyGroup.scale.setScalar(body.scale);
    bodyGroup.visible = S.showBody;
    needsVoxelize = true;
    resetAverages(90);
    updateRoadNote();
  }

  /** Re-tint the 3D scene for the current theme (also after a rebuild). */
  function applyThemeColors() {
    const c = THEME[S.theme];
    if (!scene) return;
    scene.background = new THREE.Color(c.bg);
    if (scene.fog) scene.fog.color = new THREE.Color(c.bg);
    if (tunnelBox) tunnelBox.material.color.setHex(c.tunnel);
    if (roadMesh) roadMesh.material.color.setHex(c.road);
    if (lineMesh) lineMesh.material.color.setHex(c.line);
  }

  function applyTheme(t) {
    S.theme = t;
    localStorage.setItem('wt-theme', t);
    document.body.classList.toggle('irix', t === 'irix');
    document.getElementById('tModern').classList.toggle('on', t === 'modern');
    document.getElementById('tIrix').classList.toggle('on', t === 'irix');
    applyThemeColors();
  }

  /** Explain what the tunnel floor is currently doing. */
  function updateRoadNote() {
    const el = document.getElementById('roadNote');
    const btn = document.getElementById('tRoad');
    if (!el || !btn) return;
    btn.classList.toggle('off', !S.ground);
    if (!S.ground) {
      el.textContent = 'No ground plane for this body.';
      el.classList.remove('warn');
      return;
    }
    if (!S.rollingRoad) {
      el.textContent = 'Fixed floor: a boundary layer grows along it and the wheels ' +
        'are stationary — the classic static-floor tunnel.';
      el.classList.remove('warn');
      return;
    }
    const gap = sim.rideHeight;
    let t = body.hasWheels
      ? 'Road sweeps past at wind speed and the wheels roll with it: contact patch ' +
      'at +U, crown at −U, so the crown meets the air at twice the wind speed. '
      : 'Road sweeps past at wind speed. No identifiable wheels on this body, so ' +
      'nothing is rotating. ';
    // At an unresolved ride height a moving road blows up: measured -2% on a
    // 5.8-cell gap but +242% on a 1.8-cell one.
    const bad = gap > 0 && gap < 4;
    t += bad
      ? `⚠ Ride height is only ${gap} cells — too coarse for the underbody shear layer, ` +
      'so this result is not trustworthy. Raise the body or use a finer grid.'
      : `Ride height ${gap} cells.`;
    el.textContent = t;
    el.classList.toggle('warn', bad);
  }

  /**
   * Switch between the single-threaded solver and the worker pool.
   *
   * The lattice already lives in SharedArrayBuffers, so switching costs
   * nothing but spawning or terminating workers — the flow state carries
   * straight over.
   */
  function setCompute(mode) {
    const why = LBMPool.unavailableReason();
    if (mode === 'multi' && why) mode = 'single';
    S.compute = mode;
    const $ = id => document.getElementById(id);
    $('cSingle').classList.toggle('on', mode === 'single');
    $('cMulti').classList.toggle('on', mode === 'multi');
    $('cMulti').classList.toggle('off', !!why);

    if (mode === 'single') {
      if (pool) { pool.dispose(); pool = null; }
      poolReady = null;
      batchInFlight = false;
      computeNote(why ? 'Single core. Multi-core ' + why + '.' : 'Single core.');
      return;
    }
    if (pool || poolReady) return;
    const n = S.workers || LBMPool.suggestedWorkers();
    computeNote('Starting ' + n + ' worker threads…');
    poolReady = new LBMPool.Pool(sim, n).start().then(p => {
      if (S.compute !== 'multi') { p.dispose(); return; }
      pool = p; pool.syncParams();
      poolReady = null;
      computeNote(n + ' worker threads on ' +
        (navigator.hardwareConcurrency || '?') + ' logical cores, one z-slab each.');
    }).catch(err => {
      poolReady = null;
      console.error(err);
      setCompute('single');
      computeNote('Worker pool failed to start — staying single core.');
    });
  }
  function computeNote(t) {
    const el = document.getElementById('computeNote');
    if (el) el.textContent = t;
  }

  /**
   * Time the solver on this machine and pick the best worker count.
   *
   * The optimum is very machine-specific — it depends on physical vs logical
   * cores and on memory bandwidth, and this kernel is bandwidth bound — so
   * measuring beats guessing. Runs on a scratch lattice so the live flow is
   * untouched.
   */
  async function benchmark() {
    const btn = document.getElementById('cBench');
    if (btn.disabled) return;
    btn.disabled = true;
    const wasRunning = S.running;
    S.running = false;
    try {
      const [nx, ny, nz] = RES[S.res];
      const build = (shared) => {
        const s = new LBM(nx, ny, nz, { shared });
        s.tau = sim.tau; s.les = sim.les; s.u0 = sim.u0;
        const rolling = S.rollingRoad && S.ground;
        s.voxelize(body.makeTest(), body.box(), S.ground,
          rolling ? body.makeWheelVel(s.u0) : null, rolling ? s.u0 : 0);
        return s;
      };
      computeNote('Benchmarking 1 core…');
      await new Promise(r => setTimeout(r, 30));
      const A = build(false);
      for (let i = 0; i < 12; i++) A.step();                  // warm the JIT
      let best1 = Infinity;
      for (let r = 0; r < 3; r++) {
        const t = performance.now();
        for (let i = 0; i < 20; i++) A.step();
        best1 = Math.min(best1, (performance.now() - t) / 20);
      }
      const mlups = ms => (A.n / ms / 1000).toFixed(1);
      let line = `1 core ${best1.toFixed(1)} ms/step (${mlups(best1)} MLUPS)`;

      const why = LBMPool.unavailableReason();
      if (!why) {
        const hc = navigator.hardwareConcurrency || 4;
        const tries = [...new Set([2, Math.max(2, hc >> 1), Math.max(2, hc - 1)])]
          .filter(v => v >= 2 && v <= 16).sort((a, b) => a - b);
        let bestN = 0, bestMs = Infinity;
        for (const nw of tries) {
          computeNote(`Benchmarking ${nw} workers…`);
          await new Promise(r => setTimeout(r, 30));
          const B = build(true);
          const p = await new LBMPool.Pool(B, nw).start();
          p.syncParams();
          await p.step(8);                                    // warm
          let ms = Infinity;
          for (let r = 0; r < 3; r++) {
            const t = performance.now();
            await p.step(8);
            ms = Math.min(ms, (performance.now() - t) / 8);
          }
          p.dispose();
          if (ms < bestMs) { bestMs = ms; bestN = nw; }
        }
        line += ` · best ${bestN} workers ${bestMs.toFixed(1)} ms/step ` +
          `(${mlups(bestMs)} MLUPS, ${(best1 / bestMs).toFixed(2)}×)`;
        if (bestMs < best1) {
          S.workers = bestN;
          if (S.compute === 'multi') { setCompute('single'); setCompute('multi'); }
        } else {
          line += ' — single core wins here, staying on it';
          setCompute('single');
        }
      }
      computeNote(line);
    } catch (e) {
      computeNote('Benchmark failed: ' + e.message);
    } finally {
      S.running = wasRunning;
      btn.disabled = false;
    }
  }

  function applySolverParams() {
    // reKnob 0..1 maps tau from viscous (stable, low Re) to nearly inviscid
    const tau = 0.5 + Math.pow(10, -1.1 - 2.3 * S.reKnob);
    sim.tau = tau;
    sim.les = S.les;
    sim.u0 = 0.06;
    if (pool) pool.syncParams();          // workers hold their own copy
  }

  // ---- particles -------------------------------------------------------
  function respawn(i, anywhere) {
    const { nx, ny, nz } = sim;
    // seed a rake that just covers the body, so the smoke actually hits it
    // instead of carpeting the whole test section
    const ey = body.ext[1] * body.scale, ez = body.ext[2] * body.scale;
    const cy = S.ground ? body.pos[1] + ey * 0.5 : body.pos[1];
    const cz = body.pos[2];
    const sy = Math.min(ny * 0.45, Math.max(3, ey * (S.ground ? 1.1 : 1.7)));
    const sz = Math.min(nz * 0.45, Math.max(3, ez * 2.0));
    const x = anywhere ? 1 + Math.random() * (nx - 4) : 1.5 + Math.random() * 3;
    const y = Math.max(1.2, Math.min(ny - 1.5, cy + (Math.random() * 2 - 1) * sy));
    const z = Math.max(1.2, Math.min(nz - 1.5, cz + (Math.random() * 2 - 1) * sz));
    pPos[i * 3] = x; pPos[i * 3 + 1] = y; pPos[i * 3 + 2] = z;
    pAge[i] = Math.random() * 40;
  }

  function updateParticles(dt) {
    const { nx, ny, nz } = sim;
    const P = S.nParticles, k = S.flowSpeed * dt;
    const col = [0, 0, 0], u0 = sim.u0;
    for (let i = 0; i < P; i++) {
      const j = i * 3;
      let x = pPos[j], y = pPos[j + 1], z = pPos[j + 2];
      sim.sample(x, y, z, tmpV);
      x += tmpV[0] * k; y += tmpV[1] * k; z += tmpV[2] * k;
      pAge[i] += dt;
      if (x > nx - 2 || x < 0.5 || y < 0.8 || y > ny - 1.2 || z < 0.8 || z > nz - 1.2 ||
        pAge[i] > 380 || sim.isSolid(x, y, z)) {
        respawn(i, false);
        x = pPos[j]; y = pPos[j + 1]; z = pPos[j + 2];
      }
      pPos[j] = x; pPos[j + 1] = y; pPos[j + 2] = z;
      const sp = Math.hypot(tmpV[0], tmpV[1], tmpV[2]) / u0;
      cmap(sp * K_SPEED, col);
      const fade = Math.min(1, pAge[i] / 12) * (1 - Math.max(0, (pAge[i] - 300) / 80));
      pCol[j] = col[0] * fade; pCol[j + 1] = col[1] * fade; pCol[j + 2] = col[2] * fade;
    }
    pGeo.attributes.position.needsUpdate = true;
    pGeo.attributes.color.needsUpdate = true;
  }

  // ---- smoke ribbons (streaklines) -------------------------------------
  // A real tunnel's smoke wand emits continuously from a fixed nozzle, and
  // what you see is the *streakline*: every particle ever released from that
  // point, joined up. That is not the same as a streamline once the flow is
  // unsteady — streaklines are the things that visibly roll up into vortices
  // behind a bluff body, which is exactly what makes tunnel footage readable.
  //
  // Each wand keeps a polyline: advect every point, then push a fresh point
  // in at the nozzle. Drawn as camera-facing quads because WebGL will not
  // give us line thickness.
  const RIB_LEN = 88;
  let ribbons, rbGeo, rbPos, rbCol, rbSeeds, rbPts, rbSpd, rbN = 0;

  function buildRibbons() {
    if (ribbons) { scene.remove(ribbons); ribbons.geometry.dispose(); }
    const ny = 7, nz = 5;
    rbN = ny * nz;
    rbSeeds = new Float32Array(rbN * 3);
    rbPts = new Float32Array(rbN * RIB_LEN * 3);
    rbSpd = new Float32Array(rbN * RIB_LEN);   // cached at advection time
    const verts = rbN * RIB_LEN * 2;
    rbPos = new Float32Array(verts * 3);
    rbCol = new Float32Array(verts * 3);
    const idx = [];
    for (let w = 0; w < rbN; w++) {
      const base = w * RIB_LEN * 2;
      for (let i = 0; i < RIB_LEN - 1; i++) {
        const a = base + i * 2;
        idx.push(a, a + 1, a + 2, a + 1, a + 3, a + 2);
      }
    }
    rbGeo = new THREE.BufferGeometry();
    rbGeo.setAttribute('position', new THREE.BufferAttribute(rbPos, 3));
    rbGeo.setAttribute('color', new THREE.BufferAttribute(rbCol, 3));
    rbGeo.setIndex(idx);
    ribbons = new THREE.Mesh(rbGeo, new THREE.MeshBasicMaterial({
      vertexColors: true, transparent: true, opacity: 0.62,
      side: THREE.DoubleSide, depthWrite: false, blending: THREE.AdditiveBlending
    }));
    ribbons.frustumCulled = false;
    scene.add(ribbons);
    seedRibbons();
  }

  /** Park the wands on a rake upstream, sized to the body. */
  function seedRibbons() {
    if (!rbN || !body) return;
    const { ny, nz } = sim;
    const ey = body.ext[1] * body.scale, ez = body.ext[2] * body.scale;
    const cy = S.ground ? body.pos[1] + ey * 0.5 : body.pos[1];
    const sy = Math.min(ny * 0.42, Math.max(3, ey * (S.ground ? 1.15 : 1.5)));
    const sz = Math.min(nz * 0.42, Math.max(3, ez * 1.5));
    const NY = 7, NZ = 5;
    for (let a = 0; a < NY; a++) {
      for (let b = 0; b < NZ; b++) {
        const w = a * NZ + b;
        const y = clamp(cy + (a / (NY - 1) - 0.5) * 2 * sy, 1.4, ny - 1.6);
        const z = clamp(cz(b, NZ, sz), 1.4, nz - 1.6);
        rbSeeds[w * 3] = 2.0; rbSeeds[w * 3 + 1] = y; rbSeeds[w * 3 + 2] = z;
        for (let i = 0; i < RIB_LEN; i++) {      // start collapsed at the nozzle
          const o = (w * RIB_LEN + i) * 3;
          rbPts[o] = 2.0; rbPts[o + 1] = y; rbPts[o + 2] = z;
        }
      }
    }
    function cz(b, NZ, sz) { return body.pos[2] + (b / (NZ - 1) - 0.5) * 2 * sz; }
  }

  function updateRibbons(dt) {
    const k = S.flowSpeed * dt;
    for (let w = 0; w < rbN; w++) {
      const b = w * RIB_LEN * 3;
      const sb = w * RIB_LEN;
      for (let i = 0; i < RIB_LEN; i++) {
        const o = b + i * 3;
        const x = rbPts[o], y = rbPts[o + 1], z = rbPts[o + 2];
        sim.sample(x, y, z, tmpV);
        rbSpd[sb + i] = Math.hypot(tmpV[0], tmpV[1], tmpV[2]);
        const nx2 = x + tmpV[0] * k, ny2 = y + tmpV[1] * k, nz2 = z + tmpV[2] * k;
        // never let smoke wander inside the body
        if (!sim.isSolid(nx2, ny2, nz2)) {
          rbPts[o] = nx2; rbPts[o + 1] = ny2; rbPts[o + 2] = nz2;
        }
      }
      // slide the trail along and emit a fresh puff at the nozzle
      rbPts.copyWithin(b + 3, b, b + (RIB_LEN - 1) * 3);
      rbSpd.copyWithin(sb + 1, sb, sb + RIB_LEN - 1);
      rbPts[b] = rbSeeds[w * 3];
      rbPts[b + 1] = rbSeeds[w * 3 + 1];
      rbPts[b + 2] = rbSeeds[w * 3 + 2];
      rbSpd[sb] = rbSpd[sb + 1];
    }
    ribbonGeometry();
  }

  /** Expand each polyline into camera-facing quads. */
  function ribbonGeometry() {
    const cx = camera.position.x, cy = camera.position.y, cz = camera.position.z;
    const half = Math.max(0.14, sim.nx / 320);
    const col = [0, 0, 0], u0 = sim.u0;
    for (let w = 0; w < rbN; w++) {
      const b = w * RIB_LEN * 3;
      let vo = w * RIB_LEN * 2 * 3;
      for (let i = 0; i < RIB_LEN; i++) {
        const o = b + i * 3;
        const x = rbPts[o], y = rbPts[o + 1], z = rbPts[o + 2];
        const pa = b + Math.max(0, i - 1) * 3, pb = b + Math.min(RIB_LEN - 1, i + 1) * 3;
        let tx = rbPts[pa] - rbPts[pb],
          ty = rbPts[pa + 1] - rbPts[pb + 1],
          tz = rbPts[pa + 2] - rbPts[pb + 2];
        let tl = Math.hypot(tx, ty, tz);
        if (tl < 1e-6) { tx = 1; ty = 0; tz = 0; tl = 1; }
        tx /= tl; ty /= tl; tz /= tl;
        // side = tangent x view, so the ribbon always faces the camera
        const vx = cx - x, vy = cy - y, vz = cz - z;
        let sx = ty * vz - tz * vy, sy = tz * vx - tx * vz, sz2 = tx * vy - ty * vx;
        const sl = Math.hypot(sx, sy, sz2) || 1;
        sx /= sl; sy /= sl; sz2 /= sl;
        // taper: emerges from the nozzle, dissolves at the tail
        const t = i / (RIB_LEN - 1);
        const grow = Math.min(1, i / 4);
        const fade = 1 - t * t;
        const hw = half * grow * (0.35 + 0.65 * fade);
        rbPos[vo] = x + sx * hw; rbPos[vo + 1] = y + sy * hw; rbPos[vo + 2] = z + sz2 * hw;
        rbPos[vo + 3] = x - sx * hw; rbPos[vo + 4] = y - sy * hw; rbPos[vo + 5] = z - sz2 * hw;
        cmap(rbSpd[w * RIB_LEN + i] / u0 * K_SPEED, col);
        // mostly white like real smoke, tinted a little by speed. Kept dim
        // because additive blending piles up wherever ribbons overlap, and
        // the model has to stay readable through them.
        const g = 0.34 * fade * grow;
        const r0 = (0.55 + 0.45 * col[0]) * g,
          g0 = (0.55 + 0.45 * col[1]) * g,
          b0 = (0.55 + 0.45 * col[2]) * g;
        rbCol[vo] = r0; rbCol[vo + 1] = g0; rbCol[vo + 2] = b0;
        rbCol[vo + 3] = r0; rbCol[vo + 4] = g0; rbCol[vo + 5] = b0;
        vo += 6;
      }
    }
    rbGeo.attributes.position.needsUpdate = true;
    rbGeo.attributes.color.needsUpdate = true;
  }

  // ---- slice texture ---------------------------------------------------
  function updateSlice() {
    const { nx, ny, nz, ux, uy, uz, rho, solid } = sim;
    const u0 = sim.u0, col = [0, 0, 0];
    const xy = S.sliceAxis === 'xy';
    const h = xy ? ny : nz;
    const k = xy ? Math.max(1, Math.min(nz - 2, Math.round(body.pos[2])))
      : Math.max(1, Math.min(ny - 2, Math.round(body.pos[1])));
    const at = xy ? (x, y) => x + nx * (y + ny * k) : (x, z) => x + nx * (k + ny * z);

    for (let b = 0; b < h; b++) {
      for (let a = 0; a < nx; a++) {
        const i = at(a, b);
        const o = (a + b * nx) * 4;
        if (solid[i]) { sliceData[o] = 20; sliceData[o + 1] = 24; sliceData[o + 2] = 30; sliceData[o + 3] = 255; continue; }
        let t;
        if (S.sliceField === 'speed') {
          t = Math.hypot(ux[i], uy[i], uz[i]) / u0 * K_SPEED;
          cmap(t, col);
        } else if (S.sliceField === 'pressure') {
          // BGK lattice Boltzmann carries a well-known odd/even checkerboard in
          // the density field. It is a lattice artefact, not pressure, so average
          // the in-plane neighbours for display — the solver is untouched.
          const ip = a < nx - 1 ? at(a + 1, b) : i, im = a > 0 ? at(a - 1, b) : i;
          const jp = b < h - 1 ? at(a, b + 1) : i, jm = b > 0 ? at(a, b - 1) : i;
          const r = (rho[i] * 2 + rho[ip] + rho[im] + rho[jp] + rho[jm]) / 6;
          const cp = ((r - 1) / 3) / (0.5 * u0 * u0);
          diverging(cp * K_CP, col);
        } else {                                   // vorticity
          let w;
          if (xy) {
            const ip = a < nx - 1 ? at(a + 1, b) : i, im = a > 0 ? at(a - 1, b) : i;
            const jp = b < h - 1 ? at(a, b + 1) : i, jm = b > 0 ? at(a, b - 1) : i;
            w = (uy[ip] - uy[im]) * 0.5 - (ux[jp] - ux[jm]) * 0.5;
          } else {
            const ip = a < nx - 1 ? at(a + 1, b) : i, im = a > 0 ? at(a - 1, b) : i;
            const jp = b < h - 1 ? at(a, b + 1) : i, jm = b > 0 ? at(a, b - 1) : i;
            w = (ux[jp] - ux[jm]) * 0.5 - (uz[ip] - uz[im]) * 0.5;
          }
          diverging(w / u0 * K_VORT, col);
        }
        sliceData[o] = col[0] * 255; sliceData[o + 1] = col[1] * 255;
        sliceData[o + 2] = col[2] * 255; sliceData[o + 3] = 235;
      }
    }
    sliceTex.needsUpdate = true;
  }

  /**
   * Colour key for the slice. Painted with the very same cmap()/diverging()
   * the slice uses, so it cannot go stale if a palette is retuned.
   */
  function drawLegend() {
    const c = document.getElementById('legendBar');
    const g = c.getContext('2d');
    const w = c.width, h = c.height;
    const img = g.createImageData(w, h);
    const col = [0, 0, 0];
    const speed = S.sliceField === 'speed';
    for (let x = 0; x < w; x++) {
      const p = x / (w - 1);
      if (speed) cmap(p, col); else diverging(p * 2 - 1, col);
      const r = col[0] * 255, gr = col[1] * 255, b = col[2] * 255;
      for (let y = 0; y < h; y++) {
        const o = (x + y * w) * 4;
        img.data[o] = r; img.data[o + 1] = gr; img.data[o + 2] = b; img.data[o + 3] = 255;
      }
    }
    if (speed) {                       // tick marking undisturbed wind speed
      const xt = Math.round(K_SPEED * (w - 1));
      for (let y = 0; y < h; y++) {
        const o = (xt + y * w) * 4;
        img.data[o] = 255; img.data[o + 1] = 255; img.data[o + 2] = 255;
      }
    }
    g.putImageData(img, 0, 0);

    const L = {
      speed: ['still', '▲ wind speed', (1 / K_SPEED).toFixed(1) + '× wind'],
      pressure: ['suction −' + (1 / K_CP).toFixed(1), 'ambient', '+' + (1 / K_CP).toFixed(1) + ' Cₚ'],
      vorticity: ['↻ clockwise', 'no spin', 'anticlockwise ↺']
    }[S.sliceField];
    document.getElementById('legMin').textContent = L[0];
    document.getElementById('legMid').textContent = L[1];
    document.getElementById('legMax').textContent = L[2];
  }

  // ---- streamlines -----------------------------------------------------
  function updateStreamlines() {
    const { nx, ny, nz } = sim;
    const NL = lineMesh._nl, LEN = lineMesh._len;
    // same rake the smoke uses, so the two visualisations agree
    const ey = body.ext[1] * body.scale, ez = body.ext[2] * body.scale;
    const cy = S.ground ? body.pos[1] + ey * 0.5 : body.pos[1];
    const sy = Math.min(ny * 0.45, Math.max(3, ey * (S.ground ? 1.1 : 1.7)));
    const sz = Math.min(nz * 0.45, Math.max(3, ez * 2.0));
    let p = 0;
    for (let l = 0; l < NL; l++) {
      const gy = Math.floor(l / 15), gz = l % 15;
      let x = 2.0;
      let y = cy + (gy / 14 - 0.5) * 2 * sy;
      let z = body.pos[2] + (gz / 14 - 0.5) * 2 * sz;
      y = Math.max(1.2, Math.min(ny - 1.5, y));
      z = Math.max(1.2, Math.min(nz - 1.5, z));
      const h = 0.9;
      for (let s = 0; s < LEN; s++) {
        const x0 = x, y0 = y, z0 = z;
        sim.sample(x, y, z, tmpV);
        const sp = Math.hypot(tmpV[0], tmpV[1], tmpV[2]);
        if (sp < 1e-6) { for (; s < LEN; s++) { linePos[p++] = x0; linePos[p++] = y0; linePos[p++] = z0; linePos[p++] = x0; linePos[p++] = y0; linePos[p++] = z0; } break; }
        const step = h / sp;
        const mx = x + tmpV[0] * step * 0.5, my = y + tmpV[1] * step * 0.5, mz = z + tmpV[2] * step * 0.5;
        sim.sample(mx, my, mz, tmpV);
        const sp2 = Math.max(1e-6, Math.hypot(tmpV[0], tmpV[1], tmpV[2]));
        x += tmpV[0] / sp2 * h; y += tmpV[1] / sp2 * h; z += tmpV[2] / sp2 * h;
        if (x > nx - 2 || x < 1 || y < 1 || y > ny - 2 || z < 1 || z > nz - 2) { x = x0; y = y0; z = z0; }
        linePos[p++] = x0; linePos[p++] = y0; linePos[p++] = z0;
        linePos[p++] = x; linePos[p++] = y; linePos[p++] = z;
      }
    }
    lineGeo.attributes.position.needsUpdate = true;
  }

  // ---- interaction -----------------------------------------------------
  function ndc(e) {
    return new THREE.Vector2((e.clientX / window.innerWidth) * 2 - 1,
      -(e.clientY / window.innerHeight) * 2 + 1);
  }
  function onPointerDown(e) {
    if (e.button !== 0) return;
    raycaster.setFromCamera(ndc(e), camera);
    const hits = raycaster.intersectObjects(bodyGroup.children, false);
    if (!hits.length) return;
    orbit.enabled = false;
    const horizontal = e.shiftKey;
    const normal = horizontal ? new THREE.Vector3(0, 1, 0) : new THREE.Vector3(0, 0, 1);
    const plane = new THREE.Plane().setFromNormalAndCoplanarPoint(
      normal, new THREE.Vector3(body.pos[0], body.pos[1], body.pos[2]));
    const hit = new THREE.Vector3();
    raycaster.ray.intersectPlane(plane, hit);
    dragging = {
      plane, horizontal,
      off: [hit.x - body.pos[0], hit.y - body.pos[1], hit.z - body.pos[2]]
    };
  }
  function onPointerMove(e) {
    if (!dragging) return;
    raycaster.setFromCamera(ndc(e), camera);
    const hit = new THREE.Vector3();
    if (!raycaster.ray.intersectPlane(dragging.plane, hit)) return;
    const { nx, ny, nz } = sim;
    const r = body.radius * body.scale + 2;
    body.pos[0] = clamp(hit.x - dragging.off[0], r + 2, nx - r - 2);
    if (dragging.horizontal) {
      body.pos[2] = clamp(hit.z - dragging.off[2], r + 2, nz - r - 2);
    } else if (!S.ground) {
      body.pos[1] = clamp(hit.y - dragging.off[1], r + 2, ny - r - 2);
    }
    syncBodyTransform();
    positionSlice();
  }
  const clamp = (v, a, b) => v < a ? a : v > b ? b : v;

  // ---- readouts ---------------------------------------------------------
  let statTick = 0;
  function updateStats() {
    // wings and aircraft are quoted on planform (wing) area, everything else
    // on projected frontal area — the usual conventions in each field
    const planform = body.def.ref === 'planform';
    const A_cells = refCells();
    const u0 = sim.u0;
    const st = runStats();
    const ready = !!st;
    const q_lat = 0.5 * u0 * u0 * A_cells;
    const Cd = ready ? st.cd.mean : forceEMA[0] / q_lat;
    const Cl = ready ? st.cl.mean : forceEMA[1] / q_lat;
    const Cs = ready ? st.cy.mean : forceEMA[2] / q_lat;

    const mPerCell = S.lengthM / body.scale;
    const A_real = A_cells * mPerCell * mPerCell;
    const U = S.windMs;
    const q = 0.5 * RHO_AIR * U * U;
    const drag = Cd * q * A_real;
    const lift = Cl * q * A_real;
    const power = drag * U;

    const reReal = U * S.lengthM / NU_AIR;
    const reSim = u0 * body.scale / sim.nu;

    const set = (id, v) => { document.getElementById(id).textContent = v; };
    const dots = settle > 0 ? 'settling…' : 'averaging…';

    // coefficient + 95% confidence half-width; both dimensionless, so the
    // unit switch does not touch them
    const setCoef = (id, mean, sem) => {
      const el = document.getElementById(id);
      el.textContent = '';
      if (!ready) { el.textContent = id === 'cd' ? dots : ''; return; }
      const v = document.createElement('span');
      v.textContent = mean.toFixed(3);
      el.appendChild(v);
      const c = document.createElement('span');
      c.className = 'ci';
      c.textContent = ' ±' + (1.96 * sem).toFixed(3);
      el.appendChild(c);
    };
    setCoef('cd', Cd, ready ? st.cd.sem : 0);
    setCoef('cl', Cl, ready ? st.cl.sem : 0);
    setCoef('cs', Cs, ready ? st.cy.sem : 0);

    const badge = document.getElementById('convState');
    badge.textContent = ready ? st.state + ' · n=' + st.n : (settle > 0 ? 'settling' : 'starting');
    badge.className = 'cstate ' + (ready ? st.state : 'averaging');
    drawTrace(st);
    set('area', fmtArea(A_real));
    document.getElementById('arealbl').textContent =
      planform ? 'Planform (ref.) area' : 'Frontal (ref.) area';
    set('q', fmtPressure(q));
    set('drag', ready ? fmtForce(drag) : '');
    // spell out which way the vertical force actually points — a bare minus
    // sign is easy to misread, and downforce is usually the number you want
    const down = lift < 0;
    set('liftlbl', ready ? (down ? 'Downforce' : 'Lift') : 'Lift / downforce');
    set('lift', ready ? (down ? '↓ ' : '↑ ') + fmtForce(Math.abs(lift)) : '');
    document.getElementById('lift').classList.toggle('down', ready && down);
    set('power', ready ? fmtPower(power) : '');
    set('reReal', fmtSci(reReal));
    set('reSim', fmtSci(reSim));
    set('kmh', fmtRoadSpeed(U));
    set('mach', (U / 343).toFixed(2));
    set('perf', `${fpsEMA.toFixed(0)} fps · ${msPerStep.toFixed(1)} ms/step · ${stepsPerFrame}×`);
    set('grid', `${sim.nx}×${sim.ny}×${sim.nz} = ${(sim.n / 1000).toFixed(0)}k cells`);
    set('tsteps', sim.steps.toLocaleString());
  }
  /** Cd/Cl against time, so you can see whether it has settled or is shedding. */
  function drawTrace(st) {
    const c = document.getElementById('trace');
    const g = c.getContext('2d');
    const w = c.width, h = c.height;
    const dark = S.theme !== 'irix';
    g.clearRect(0, 0, w, h);
    g.fillStyle = dark ? 'rgba(255,255,255,.035)' : 'rgba(0,0,0,.10)';
    g.fillRect(0, 0, w, h);
    if (!st) return;

    const n = ordered(hCd, scratch); const cd = scratch.slice(0, n);
    const nl = ordered(hCl, scratch); const cl = scratch.slice(0, nl);

    let lo = Infinity, hi = -Infinity;
    for (let i = 0; i < n; i++) { if (cd[i] < lo) lo = cd[i]; if (cd[i] > hi) hi = cd[i]; }
    for (let i = 0; i < nl; i++) { if (cl[i] < lo) lo = cl[i]; if (cl[i] > hi) hi = cl[i]; }
    if (!(hi > lo)) hi = lo + 1;
    const pad = (hi - lo) * 0.12;
    lo -= pad; hi += pad;
    const Y = v => h - (v - lo) / (hi - lo) * h;

    if (lo < 0 && hi > 0) {
      g.strokeStyle = dark ? 'rgba(255,255,255,.18)' : 'rgba(0,0,0,.25)';
      g.beginPath(); g.moveTo(0, Y(0)); g.lineTo(w, Y(0)); g.stroke();
    }
    const plot = (a, len, colour) => {
      if (len < 2) return;
      g.strokeStyle = colour; g.lineWidth = 1; g.beginPath();
      for (let i = 0; i < len; i++) {
        const x = i / (len - 1) * w, y = Y(a[i]);
        i ? g.lineTo(x, y) : g.moveTo(x, y);
      }
      g.stroke();
    };
    plot(cl, nl, dark ? 'rgba(255,180,84,.7)' : 'rgba(120,20,70,.7)');
    plot(cd, n, dark ? 'rgba(79,208,255,.9)' : 'rgba(20,16,90,.85)');

    const band = 1.96 * st.cd.sem;                  // mean and its 95% band
    g.fillStyle = dark ? 'rgba(79,208,255,.16)' : 'rgba(20,16,90,.14)';
    g.fillRect(0, Y(st.cd.mean + band), w,
      Math.max(1, Y(st.cd.mean - band) - Y(st.cd.mean + band)));
    g.strokeStyle = dark ? 'rgba(79,208,255,.55)' : 'rgba(20,16,90,.55)';
    g.beginPath(); g.moveTo(0, Y(st.cd.mean)); g.lineTo(w, Y(st.cd.mean)); g.stroke();
  }

  /** A run report a technician can archive: conditions, results, raw series. */
  function buildRunCSV() {
    const st = runStats();
    const planform = body.def.ref === 'planform';
    const A_cells = refCells();
    const mPerCell = S.lengthM / body.scale;
    const A_real = A_cells * mPerCell * mPerCell;
    const U = S.windMs, q = 0.5 * RHO_AIR * U * U;
    const L = [];
    const p = (k, v, u) => L.push([k, v, u === undefined ? '' : u].join(','));

    L.push('Wind Tunnel run report,,');
    L.push('parameter,value,unit');
    p('generated', new Date().toISOString());
    p('body', JSON.stringify(body.def.name));
    p('body_source', body.def.custom ? 'uploaded' : 'bundled');
    if (model && body.key === model.key && model.credit) {
      p('model_title', JSON.stringify(model.credit.title || ''));
      p('model_author', JSON.stringify(model.credit.author || ''));
      p('model_license', JSON.stringify(model.credit.license || ''));
      p('model_source', JSON.stringify(model.credit.source || ''));
    }
    p('grid_nx', sim.nx, 'cells'); p('grid_ny', sim.ny, 'cells'); p('grid_nz', sim.nz, 'cells');
    p('grid_total', sim.n, 'cells');
    p('solid_cells', sim.solidCells, 'cells');
    p('reference_area_type', planform ? 'planform' : 'frontal');
    p('reference_area', A_cells, 'cells');
    p('reference_area_real', A_real.toFixed(4), 'm^2');
    p('blockage', (100 * sim.frontalArea / (sim.ny * sim.nz)).toFixed(2), '%');
    p('body_length', S.lengthM, 'm');
    p('cells_per_body_length', body.scale.toFixed(2), 'cells');
    p('wind_speed', U, 'm/s');
    p('wind_speed_kmh', (U * 3.6).toFixed(1), 'km/h');
    p('mach', (U / 343).toFixed(4));
    p('dynamic_pressure', q.toFixed(2), 'Pa');
    p('angle_of_attack', S.aoa, 'deg');
    p('yaw', S.yaw, 'deg');
    p('ground_plane', body.def.ground ? 'yes' : 'no');
    p('lattice_u0', sim.u0);
    p('lattice_tau', sim.tau.toFixed(6));
    p('lattice_nu', sim.nu.toExponential(4));
    p('les_smagorinsky', sim.les ? 'on' : 'off');
    p('smagorinsky_C', sim.csm);
    p('reynolds_simulated', Math.round(sim.u0 * body.scale / sim.nu));
    p('reynolds_full_scale', Math.round(U * S.lengthM / NU_AIR));
    p('timesteps_total', sim.steps);
    p('samples', st ? st.n : 0);
    p('convergence', st ? st.state : 'not started');

    L.push(',,');
    L.push('result,mean,ci95,unit');
    if (st) {
      const row = (k, m, sem, u) => L.push([k, m.toFixed(5), (1.96 * sem).toFixed(5), u].join(','));
      row('Cd', st.cd.mean, st.cd.sem, '-');
      row('Cl', st.cl.mean, st.cl.sem, '-');
      row('Cy', st.cy.mean, st.cy.sem, '-');
      row('drag', st.cd.mean * q * A_real, st.cd.sem * q * A_real, 'N');
      row('lift', st.cl.mean * q * A_real, st.cl.sem * q * A_real, 'N');
      row('side', st.cy.mean * q * A_real, st.cy.sem * q * A_real, 'N');
      row('power', st.cd.mean * q * A_real * U, st.cd.sem * q * A_real * U, 'W');
    }
    L.push(',,');
    L.push('# ci95 is from the blocking method (16 blocks) which accounts for,,');
    L.push('# correlation between consecutive timesteps; sigma/sqrt(N) would,,');
    L.push('# understate it badly in a shedding wake.,,');
    L.push(',,');
    L.push('sample,Cd,Cl,Cy');
    const n = ordered(hCd, scratch); const a = scratch.slice(0, n);
    const nl = ordered(hCl, scratch); const b = scratch.slice(0, nl);
    const ny = ordered(hCy, scratch); const c2 = scratch.slice(0, ny);
    for (let i = 0; i < n; i++) {
      L.push(i + ',' + a[i].toFixed(5) + ',' + b[i].toFixed(5) + ',' + c2[i].toFixed(5));
    }
    return L.join('\n');
  }

  function exportRun() {
    const slug = (body.def.name || 'body').replace(/[^\w-]+/g, '_').slice(0, 40);
    const stamp = new Date().toISOString().replace(/[:.]/g, '-').slice(0, 19);
    const blob = new Blob([buildRunCSV()], { type: 'text/csv;charset=utf-8' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `windtunnel_${slug}_${S.windMs}ms_aoa${S.aoa}_${stamp}.csv`;
    document.body.appendChild(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 2000);
  }

  // ---- units ------------------------------------------------------------
  // Everything is computed in SI internally; only the display converts, so
  // switching systems can never drift the physics.
  const MPH = 2.2369363, FT = 3.2808399, LBF = 0.22480894, HP = 1 / 745.6999,
    FT2 = 10.763910, PSF = 0.020885434;
  const imp = () => S.units === 'imperial';

  function fmtWind(ms) {                    // wind-speed slider
    return imp() ? Math.round(ms * MPH) + ' mph' : ms + ' m/s';
  }
  function fmtLength(m) {                   // body-length slider
    return imp() ? (m * FT).toFixed(1) + ' ft' : m + ' m';
  }
  function fmtRoadSpeed(ms) {               // "Speed" readout
    return imp() ? Math.round(ms * MPH) + ' mph' : Math.round(ms * 3.6) + ' km/h';
  }
  function fmtForce(N) {
    if (imp()) {
      const v = N * LBF;
      return Math.abs(v) >= 1000 ? (v / 1000).toFixed(2) + ' klbf' : v.toFixed(1) + ' lbf';
    }
    return Math.abs(N) >= 1000 ? (N / 1000).toFixed(2) + ' kN' : N.toFixed(1) + ' N';
  }
  function fmtPower(W) {
    return imp() ? (W * HP).toFixed(1) + ' hp' : (W / 1000).toFixed(1) + ' kW';
  }
  function fmtArea(m2) {
    return imp() ? (m2 * FT2).toFixed(2) + ' ft²' : m2.toFixed(2) + ' m²';
  }
  function fmtPressure(Pa) {
    // lb/ft^2, not psi — psi puts wind-tunnel q at 0.08 and throws away
    // every useful digit
    return imp() ? (Pa * PSF).toFixed(1) + ' lb/ft²' : (Pa / 1000).toFixed(2) + ' kPa';
  }
  function fmtSci(v) {
    if (v < 1e4) return v.toFixed(0);
    const e = Math.floor(Math.log10(v));
    return (v / Math.pow(10, e)).toFixed(1) + '×10' + sup(e);
  }
  const SUP = '⁰¹²³⁴⁵⁶⁷⁸⁹';
  function sup(n) { return String(n).split('').map(c => SUP[+c]).join(''); }

  // ---- main loop ---------------------------------------------------------
  function animate() {
    requestAnimationFrame(animate);
    const now = performance.now();
    const dt = Math.min(0.05, (now - lastT) / 1000);
    lastT = now;
    fpsEMA += (1 / Math.max(1e-3, dt) - fpsEMA) * 0.08;

    if (needsVoxelize) {
      // a rolling road and spinning wheels are one physical setup: a car on a
      // road rather than a model bolted to a static tunnel floor
      const rolling = S.rollingRoad && S.ground;
      const wheelVel = rolling ? body.makeWheelVel(sim.u0) : null;
      sim.voxelize(body.makeTest(), body.box(), S.ground, wheelVel, rolling ? sim.u0 : 0);
      needsVoxelize = false;
      if (pool) pool.syncParams();        // movingWalls may have changed
      updateRoadNote();
    }

    if (S.running && pool) {
      // Multi-core: fire a batch and let it run while this frame renders the
      // previous state. Rendering no longer waits on the solver at all.
      if (!batchInFlight && !needsVoxelize) {
        batchInFlight = true;
        const n = Math.max(1, Math.min(8, stepsPerFrame));
        const t0 = performance.now();
        pool.step(n).then(series => {
          for (const f of series) {
            sim.force[0] = f[0]; sim.force[1] = f[1]; sim.force[2] = f[2];
            accumulate();
          }
          const el = (performance.now() - t0) / n;
          msPerStep += (el - msPerStep) * 0.15;
          stepsPerFrame = msPerStep < 3 ? 8 : msPerStep < 6 ? 4 : msPerStep < 12 ? 2 : 1;
          batchInFlight = false;
        }).catch(err => {
          batchInFlight = false;
          console.error('worker pool failed, falling back to single core', err);
          setCompute('single');
        });
      }
    } else if (S.running) {
      const t0 = performance.now();
      for (let s = 0; s < stepsPerFrame; s++) { sim.step(); accumulate(); }
      const el = (performance.now() - t0) / stepsPerFrame;
      msPerStep += (el - msPerStep) * 0.15;
      stepsPerFrame = msPerStep < 6 ? 3 : msPerStep < 12 ? 2 : 1;
    }

    if (S.smoke) updateParticles(Math.min(2.0, dt * 60));
    if (S.ribbons) updateRibbons(Math.min(2.0, dt * 60));
    if (S.slice) updateSlice();
    if (S.streamlines && (sim.steps % 4 === 0 || !S.running)) updateStreamlines();

    if (++statTick % 12 === 0) updateStats();
    renderer.render(scene, camera);
  }

  // ---- uploaded models -----------------------------------------------------
  // Raw triangles are kept so the model can be re-oriented without re-reading
  // the file; every re-orient re-runs the voxeliser.
  let model = null;   // {name, raw, up, fwd, flip, ground, key}
  let modelSeq = 0;

  function statusModel(msg) { document.getElementById('modelInfo').textContent = msg; }

  /** First http(s) URL inside a string, or null. Rejects javascript: etc. */
  function httpUrl(s) {
    if (!s) return null;
    const m = /https?:\/\/[^\s)>\]]+/.exec(s);
    if (!m) return null;
    try {
      const u = new URL(m[0]);
      return (u.protocol === 'http:' || u.protocol === 'https:') ? u.href : null;
    } catch (e) { return null; }
  }

  /**
   * Render the credits a model file carries. Licences like CC-BY require the
   * author to be named wherever the work is shown, so this stays on screen
   * for as long as the model is loaded.
   *
   * The text comes from an untrusted downloaded file, so every value goes in
   * as textContent and links are only built from validated http(s) URLs.
   */
  function showCredit(c) {
    const box = document.getElementById('modelCredit');
    box.textContent = '';
    box.removeAttribute('title');
    if (!c) { box.classList.add('hide'); return; }

    // "Ameer Studio (https://…)" reads better as just the name
    const clean = (s) => s.replace(/\s*\(https?:\/\/[^\s)]*\)\s*/g, '').trim();

    const add = (text, url, cls) => {
      if (!text) return false;
      if (box.childNodes.length) {
        const sep = document.createElement('span');
        sep.className = 'csep'; sep.textContent = '·';
        box.appendChild(sep);
      }
      let el;
      if (url) {
        el = document.createElement('a');
        el.href = url; el.target = '_blank'; el.rel = 'noopener noreferrer';
      } else {
        el = document.createElement('span');
      }
      if (cls) el.className = cls;
      el.textContent = text;
      box.appendChild(el);
      return true;
    };

    const src = httpUrl(c.source);
    let any = false;
    any = add(c.title && clean(c.title), src, 'ctitle') || any;
    any = add(c.author && 'by ' + clean(c.author), httpUrl(c.author)) || any;
    any = add(c.license && clean(c.license), httpUrl(c.license)) || any;

    if (!any) {                                   // unstructured header text
      const t = (c.copyright || c.notes || '').replace(/\s+/g, ' ').trim();
      if (t) {
        any = add(t.length > 110 ? t.slice(0, 107) + '…' : t, null);
        box.title = t;                            // full text on hover
      }
    } else if (!c.title && src) {
      add('source', src);
    }

    box.classList.toggle('hide', !any);
  }

  function loadModelFile(file) {
    const reader = new FileReader();
    statusModel('Reading ' + file.name + ' …');
    document.getElementById('modelBox').classList.remove('hide');
    reader.onerror = () => statusModel('Could not read the file.');
    reader.onload = () => {
      let raw;
      try {
        raw = MeshIO.parse(file.name, reader.result);
      } catch (e) {
        statusModel('✗ ' + e.message);
        return;
      }
      // one upload slot: drop any previous one now that parsing succeeded
      Object.keys(Shapes.BODIES).forEach(k => {
        if (Shapes.BODIES[k].custom) delete Shapes.BODIES[k];
      });
      let credit = null;
      try { credit = MeshIO.meta(file.name, reader.result); } catch (e) { credit = null; }
      showCredit(credit);
      const b = Voxel.bounds(raw);
      const size = [b[3] - b[0], b[4] - b[1], b[5] - b[2]];
      const guess = Shapes.guessAxes(size, raw);
      model = {
        name: file.name, raw, size, credit,
        up: guess.up, fwd: guess.forward, flip: false,
        ground: true, key: 'upload'
      };
      modelSeq++;
      document.getElementById('mUp').value = 'auto';
      document.getElementById('mFwd').value = 'auto';
      document.getElementById('mFlip').classList.remove('on');
      document.getElementById('mGround').classList.add('on');
      rebuildModel();
    };
    reader.readAsArrayBuffer(file);
  }

  /** (Re)orient + voxelise the current upload and select it. */
  function rebuildModel() {
    if (!model) return;
    const t0 = performance.now();
    statusModel('Voxelising ' + (model.raw.length / 9 | 0).toLocaleString() + ' triangles…');
    // orient() works in place, so always start from a fresh copy
    const P = model.raw.slice();
    let def;
    try {
      def = Shapes.addModel(model.key, model.name.replace(/\.[^.]+$/, ''), P, {
        forward: model.fwd, up: model.up, flip: model.flip,
        ground: model.ground, lengthM: model.ground ? 4.5 : 10, res: 96
      });
    } catch (e) {
      statusModel('✗ ' + e.message);
      return;
    }
    const ms = performance.now() - t0;
    const st = def.vox.stats;
    const open = st.openness;
    statusModel(
      `${st.tris.toLocaleString()} triangles · ${st.solid.toLocaleString()} voxels · ${ms.toFixed(0)} ms\n` +
      `source ${model.size.map(v => v.toPrecision(3)).join(' × ')} · ` +
      `length=${'XYZ'[model.fwd]}${model.flip ? '−' : '+'} up=${'XYZ'[model.up]}\n` +
      (open > 0.35
        ? '⚠ mesh looks open or self-intersecting — the solid shape may be wrong'
        : open > 0.12 ? 'mesh is slightly leaky but usable' : 'mesh looks watertight'));

    refreshShapeChips();
    selectShape(model.key, true);
  }

  // ---- UI wiring ----------------------------------------------------------
  let refreshShapeChips = () => { }, selectShape = () => { };

  function ui() {
    const $ = id => document.getElementById(id);

    // shape buttons
    const list = $('shapes');
    selectShape = (k) => {
      S.shape = k; S.aoa = 0; S.yaw = 0;
      $('aoa').value = 0; $('aoav').textContent = '0°';
      $('yaw').value = 0; $('yawv').textContent = '0°';
      [...list.children].forEach(c => c.classList.toggle('on', c.dataset.k === k));
      // credits belong to the upload, not to the bundled bodies
      showCredit(model && k === model.key ? model.credit : null);
      buildBody(true);
      sim.reset();
      resetAverages();
      buildSlice();
    };
    refreshShapeChips = () => {
      list.innerHTML = '';
      Object.keys(Shapes.BODIES).forEach(k => {
        const b = document.createElement('button');
        b.textContent = Shapes.BODIES[k].name;
        b.dataset.k = k;
        b.classList.toggle('on', k === S.shape);
        b.onclick = () => selectShape(k);
        list.appendChild(b);
      });
    };
    refreshShapeChips();

    // ---- model upload ----
    $('loadModel').onclick = $('dropHelp').onclick = () => $('modelFile').click();
    $('modelFile').onchange = e => {
      if (e.target.files[0]) loadModelFile(e.target.files[0]);
      e.target.value = '';                       // let the same file reload
    };
    // "up" and "length" must stay on different axes — whichever the user just
    // set wins, and the other one moves to the best remaining axis
    const others = (skip) => [0, 1, 2].filter(a => a !== skip);
    $('mUp').onchange = e => {
      if (!model) return;
      model.up = e.target.value === 'auto'
        ? Shapes.guessAxes(model.size, model.raw).up : +e.target.value;
      if (model.fwd === model.up) {                    // longest remaining axis
        model.fwd = others(model.up).sort((a, b) => model.size[b] - model.size[a])[0];
      }
      rebuildModel();
    };
    $('mFwd').onchange = e => {
      if (!model) return;
      model.fwd = e.target.value === 'auto'
        ? Shapes.guessAxes(model.size, model.raw).forward : +e.target.value;
      if (model.up === model.fwd) {                    // shortest remaining axis
        model.up = others(model.fwd).sort((a, b) => model.size[a] - model.size[b])[0];
      }
      rebuildModel();
    };
    // Yaw the model in 90° steps about its up axis. With up fixed there are
    // only two horizontal axes, so the four orientations are (A,B) × (nose
    // forward / reversed) — pressing four times returns to where you started.
    $('mRot').onclick = () => {
      if (!model) return;
      const h = others(model.up);
      const cycle = [[h[0], false], [h[1], false], [h[0], true], [h[1], true]];
      let k = cycle.findIndex(s => s[0] === model.fwd && s[1] === model.flip);
      const next = cycle[((k < 0 ? 0 : k) + 1) % 4];
      model.fwd = next[0]; model.flip = next[1];
      $('mFlip').classList.toggle('on', model.flip);
      rebuildModel();
    };
    $('mFlip').onclick = () => {
      if (!model) return;
      model.flip = !model.flip;
      $('mFlip').classList.toggle('on', model.flip);
      rebuildModel();
    };
    $('mGround').onclick = () => {
      if (!model) return;
      model.ground = !model.ground;
      $('mGround').classList.toggle('on', model.ground);
      rebuildModel();
    };

    // drag & drop anywhere on the window
    const dz = $('dropzone');
    let depth = 0;
    window.addEventListener('dragenter', e => {
      e.preventDefault(); if (++depth === 1) dz.classList.remove('hide');
    });
    window.addEventListener('dragover', e => e.preventDefault());
    window.addEventListener('dragleave', () => {
      if (--depth <= 0) { depth = 0; dz.classList.add('hide'); }
    });
    window.addEventListener('drop', e => {
      e.preventDefault(); depth = 0; dz.classList.add('hide');
      const f = e.dataTransfer && e.dataTransfer.files && e.dataTransfer.files[0];
      if (f) loadModelFile(f);
    });

    const relabel = [];   // re-run on a unit switch
    const slider = (id, fn, fmt) => {
      const el = $(id), out = $(id + 'v');
      const upd = () => { const v = +el.value; if (out) out.textContent = fmt(v); fn(v); };
      el.addEventListener('input', upd); upd();
      relabel.push(() => { if (out) out.textContent = fmt(+el.value); });
    };

    slider('wind', v => { S.windMs = v; }, fmtWind);
    slider('lengthM', v => { S.lengthM = v; }, fmtLength);
    slider('aoa', v => { S.aoa = v; body.pitch = -v * Math.PI / 180; syncBodyTransform(); }, v => v + '°');
    slider('yaw', v => { S.yaw = v; body.yaw = v * Math.PI / 180; syncBodyTransform(); }, v => v + '°');
    slider('re', v => { S.reKnob = v / 100; applySolverParams(); resetAverages(); }, v => v + '%');
    slider('smokeN', v => { S.nParticles = v; buildParticles(sim.nx, sim.ny, sim.nz); }, v => v.toLocaleString());
    slider('vspeed', v => { S.flowSpeed = v; }, v => v + '×');

    const toggle = (id, fn, init) => {
      const el = $(id); el.classList.toggle('on', init);
      el.onclick = () => { const v = !el.classList.contains('on'); el.classList.toggle('on', v); fn(v); };
    };
    toggle('tSmoke', v => { S.smoke = v; points.visible = v; }, S.smoke);
    toggle('tRibbons', v => {
      S.ribbons = v; ribbons.visible = v;
      if (v) seedRibbons();
    }, S.ribbons);
    toggle('tSlice', v => { S.slice = v; sliceMesh.visible = v; }, S.slice);
    toggle('tLines', v => { S.streamlines = v; lineMesh.visible = v; }, S.streamlines);
    toggle('tBody', v => { S.showBody = v; bodyGroup.visible = v; }, S.showBody);
    toggle('tLes', v => { S.les = v; sim.les = v; }, S.les);
    toggle('tRoad', v => {
      S.rollingRoad = v;
      needsVoxelize = true; resetAverages();
      updateRoadNote();
    }, S.rollingRoad);
    lineMesh.visible = S.streamlines;

    const setUnits = (u) => {
      S.units = u;
      localStorage.setItem('wt-units', u);
      $('uMetric').classList.toggle('on', u === 'metric');
      $('uImperial').classList.toggle('on', u === 'imperial');
      relabel.forEach(f => f());
      updateStats();
    };
    $('uMetric').onclick = () => setUnits('metric');
    $('uImperial').onclick = () => setUnits('imperial');
    setUnits(S.units);

    $('tModern').onclick = () => applyTheme('modern');
    $('tIrix').onclick = () => applyTheme('irix');
    applyTheme(S.theme);

    $('sliceAxis').onchange = e => { S.sliceAxis = e.target.value; buildSlice(); };
    $('sliceField').onchange = e => { S.sliceField = e.target.value; drawLegend(); };
    drawLegend();
    $('cSingle').onclick = () => setCompute('single');
    $('cMulti').onclick = () => setCompute('multi');
    $('cBench').onclick = benchmark;
    setCompute(S.compute);

    $('resSel').onchange = e => {
      S.res = e.target.value;
      buildWorld();                       // drops the pool; re-spawn on the new grid
      setCompute(S.compute);
    };

    $('pause').onclick = () => {
      S.running = !S.running;
      $('pause').textContent = S.running ? '❚❚  Pause' : '▶  Run';
    };
    $('reset').onclick = () => {
      sim.reset(); resetAverages();
      for (let i = 0; i < S.nParticles; i++) respawn(i, true);
      seedRibbons();
    };
    $('exportRun').onclick = exportRun;
    $('help').onclick = () => $('helpBox').classList.toggle('hide');
    $('closeHelp').onclick = () => $('helpBox').classList.add('hide');

    window.addEventListener('keydown', e => {
      if (e.key === ' ') { e.preventDefault(); $('pause').click(); }
      if (e.key === 'r') $('reset').click();
    });
  }

  // ---- go ------------------------------------------------------------------
  function start() {
    if (!window.THREE) {
      document.body.innerHTML =
        '<div style="color:#ddd;font:15px system-ui;padding:40px">' +
        'three.js failed to load. Make sure <code>vendor/three.min.js</code> sits next to index.html.</div>';
      return;
    }
    initScene();
    buildWorld();
    ui();
    animate();
    window.WT = {
      S, orbit,
      get sim() { return sim; },
      get body() { return body; },
      get avg() { return { settle, fCount, force: forceEMA.slice() }; },
      get stats() { return runStats(); },
      get ribbons() { return ribbons; },
      stepRibbons: updateRibbons,
      buildRunCSV, exportRun
    };
  }
  window.addEventListener('DOMContentLoaded', start);
})();
