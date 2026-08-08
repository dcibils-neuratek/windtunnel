# Wind Tunnel

A 3D virtual wind tunnel in the browser. The air is solved with a real
**lattice Boltzmann (D3Q19) Navier–Stokes solver** — no faked vector fields, no
scripted particle paths. Drag a car, an airliner, a wing or a sphere around the
test section and watch the wake react.

Built by **Diego Cibils** with the help of **Claude**, for learning and
experimentation. It is a teaching and exploration tool, not a
certification-grade CFD package — read [Honest limitations](#honest-limitations)
before trusting a number.

## Run it

Just open `index.html` — everything is vendored, no build step, no network needed.

If you plan to edit the source, serve it instead. Browsers cache `.js`/`.css`
hard enough that you will reload and get the old file; `serve.py` is
`http.server` plus `Cache-Control: no-store`:

```bash
python serve.py
```

then visit <http://localhost:8777>.

## Three ways to see the flow

* **Smoke** — a cloud of tracer particles, coloured by speed. Good for seeing
  where air is fast, slow or trapped.
* **Ribbons** — smoke-wand **streaklines**, the thick continuous bands you see
  in wind-tunnel footage.
* **Streamlines** — thin curves traced through the *instantaneous* field.

Ribbons and streamlines are not the same thing, and the difference matters. A
streamline is tangent to the velocity field frozen at one instant. A streakline
is the trail of everything ever released from one fixed nozzle — which is
exactly what a real smoke wand emits. In steady flow the two coincide; in an
unsteady wake they do not, and it is the *streakline* that visibly rolls up
into vortices. That roll-up is what makes tunnel footage readable.

Each wand keeps a polyline: every point is advected each frame, then a fresh
point is pushed in at the nozzle. They are drawn as camera-facing quads, since
WebGL will not give lines real thickness, tapering from a point at the nozzle
and dissolving at the tail. 35 wands × 88 points costs about 1.2 ms a frame.

## Reading the slice

Three different fields, two different colour maps — there is a legend under the
slice selectors that repaints itself from the same palette the slice uses.

* **Speed** (sequential: navy → cyan → yellow → red). Navy is *stopped* air,
  yellow-green is the freestream (the legend marks it with a tick), red is
  ~1.6× the wind speed. Blue here means slow, **not** low pressure.
* **Pressure** (diverging: blue → black → red). Blue is negative C_p — suction.
  Red is positive C_p — compression. Black is ambient.
* **Vorticity** (diverging). Blue and red are opposite *directions* of spin.

Speed and pressure look like inverses of each other wherever the flow is
attached — that is just Bernoulli. Where they stop being inverses is the wake:
air there is both slow *and* low pressure, because separated flow dissipates
energy rather than trading it. That failure is exactly what pressure drag is.

The pressure view spatially averages its four in-plane neighbours before
colouring. BGK lattice Boltzmann carries a well-known odd/even checkerboard in
density; that is a lattice artefact rather than pressure, and smoothing it is
display-only — the solver and the reported forces are untouched.

## Testing your own models

Drop a file anywhere on the window, or use **Load model…**.

| format | notes |
|---|---|
| `.obj` | text, polygons fan-triangulated, materials ignored |
| `.stl` | binary and ASCII |
| `.glb` | binary glTF 2.0, node transforms applied, geometry only |

Not supported: `.fbx` (proprietary, version-fragmented), `.dae` (bloated, every
exporter disagrees), `.gltf` (references external files — export `.glb`), and
Draco/meshopt-compressed `.glb` (re-export with compression off; you get a clear
error message rather than a silent failure).

A wind tunnel needs triangles and nothing else — no materials, UVs or animation
— which is why the simple formats are the *right* ones here, not a compromise.

### The pipeline

```
file → parse → triangle soup → auto-orient → normalise
     → voxelise ONCE into a body-space occupancy grid
     → drag/rotate = transform + O(1) lookup per lattice cell
```

Voxelising on every drag frame would be hopeless for a 200k-triangle model, so
it happens once at load. After that a moving model costs exactly what an
analytic one does, independent of triangle count.

**Voxelisation** is ray parity along all three axes with a majority vote, OR'd
with a rasterised surface shell. Parity alone is exact for a clean watertight
mesh but degenerates on the self-intersecting, hole-ridden meshes people
actually download; voting across three independent axes recovers most of that.
The fraction of cells where the axes disagree is a real watertightness metric,
and the panel reports it — `mesh looks watertight`, `slightly leaky but usable`,
or a warning that the solid shape may be wrong. It never silently simulates
something that isn't your model.

The result is deliberately **conservative**: every cell the surface passes
through is solid. That costs about half a cell of outward bias and guarantees
that features thinner than one cell — spoilers, mirrors, wings — block flow
instead of vanishing.

**Orientation** is guessed from *silhouette area*, not bounding-box extents.
Ranking by size fails on aircraft, where span and length are within a few
percent — a coin flip. But a vehicle is shaped to present the least area along
its direction of travel and the most seen from above, and for a closed surface
the silhouette along `d` is exactly `½·Σ|Aᵢ·d|` over triangle area vectors:
one cheap pass, no rasterisation. Falls back to extents when the three
projections are too close to call (a sphere has no meaningful "forward").

No heuristic is perfect, so there are manual controls: **⟳ Rotate 90°** cycles
the four yaw orientations about the up axis (four presses returns to start),
**Flip 180°** reverses the nose, the *Up axis* / *Nose axis* selectors force a
specific axis, and *On the road* rests the model on the floor instead of
centring it. The resolved axes are always shown so you can see what was picked.

Model **size in the tunnel** comes from the frontal area measured off the voxel
grid, not an assumed bounding-box fill — a thin-winged airliner fills ~20% of
its box and a brick fills 100%, and guessing wrong throws the blockage off.

**Attribution** is read out of the file and shown for as long as the model is
loaded. glTF has a standard place for it — `asset.copyright` and
`asset.extras`, which is what Sketchfab populates — so a CC-BY model arrives
carrying its own title, author, licence and source URL. OBJ and ASCII STL only
have comment headers and binary STL an 80-byte one, so those are parsed
best-effort for `Title:` / `Author:` / `License:` / `Source:` lines. Nothing is
shown when the file says nothing.

That text is untrusted input from a downloaded file, so it is rendered
exclusively via `textContent`, and links are only built from URLs that parse as
`http`/`https` (a `javascript:` source is shown as plain text, never as an
anchor) with `rel="noopener noreferrer"`.

Set the real vehicle length with the **Body length** slider — that is what
converts coefficients into newtons.

### Validation

`testmodels.html` generates meshes with known exact volumes, runs them through
the real parse → orient → voxelise path, and checks the result against analytic
truth. All three parsers reproduce the same geometry bit-for-bit, and an
uploaded mesh cube lands within 1.1% of the built-in analytic cube
(C_d 1.371 vs 1.386, frontal area identical at 144 cells).

## Look and feel

Toggle at the top of the control panel, remembered in `localStorage`:

* **Modern** — dark translucent panels, rounded corners, cyan accents.
* **IRIX** — Silicon Graphics Indigo Magic / 4Dwm: Motif 3D relief (light
  top-left, dark bottom-right, inverted when pressed), square corners,
  purple-grey panels, teal entry fields, magenta action buttons, indigo window
  title bars with the 4Dwm minimise bar and box widgets. The 3D scene retints
  too — indigo background, lavender tunnel frame.

It is a pure theme: one `body.irix` class plus a colour table for the WebGL
scene. No layout or behaviour differs between the two.

## Units

**Metric (default)** — m/s, m, km/h, N/kN, kW, m², kPa.
**Imperial** — mph, ft, lbf, hp, ft², lb/ft².

Toggle in the *Flow* section; the choice is remembered in `localStorage`.
Everything is computed in SI internally and only the display converts, so
switching can never perturb the physics. Coefficients (`C_d`, `C_l`, `C_y`),
Mach and Reynolds number are dimensionless and read the same either way.

## Controls

| | |
|---|---|
| **Drag the body** | move it upstream/downstream and vertically |
| **Shift-drag the body** | move it sideways |
| **Left-drag empty space** | orbit |
| **Right-drag** | pan |
| **Scroll** | zoom |
| **Space / R** | pause / reset the flow |

## What is actually being computed

Each of the ~186 000 lattice cells carries 19 particle-distribution functions
`f_q`. Every timestep they stream to their neighbours and relax toward the
Maxwell–Boltzmann equilibrium:

```
f_q(x + c_q Δt, t + Δt) = f_q(x, t) − (1/τ) · [ f_q − f_q^eq ]

f_q^eq = w_q ρ [ 1 + 3(c_q·u) + 4.5(c_q·u)² − 1.5 u² ]
```

Chapman–Enskog analysis shows this recovers the **incompressible Navier–Stokes
equations** at low Mach number, with kinematic viscosity
`ν = c_s²(τ − ½)` and `c_s² = ⅓`. Mass and momentum are conserved exactly, and
pressure follows from density as `p = ρc_s²` — so unlike a projection-method
solver there is no pressure Poisson equation to invert.

* **Walls** — halfway bounce-back gives a true no-slip surface. The body is
  voxelised from the *same* analytic primitives used to build the render mesh,
  so the solver feels exactly the shape on screen. It is re-voxelised live as
  you drag.
* **Forces** — momentum exchange: every population that bounces off the surface
  deposits `2·c_q·f_q` of momentum. Summing over all boundary links gives drag,
  lift and side force directly, and `C_d = 2F / (ρu²A)`.
* **Turbulence** — a Smagorinsky sub-grid (LES) model adds eddy viscosity where
  the strain rate is high. Without it the wake blows up as soon as τ → ½.
* **Boundaries** — velocity inlet, zero-gradient outlet, free-stream far field,
  optional no-slip road for ground vehicles.

## Rolling road and rotating wheels

Optional (off by default). A real car meets still air over a moving road; a
model bolted to a static tunnel floor does not, which is why full-scale tunnels
run a rolling belt. Turning it on makes the floor sweep past at wind speed and
spins the wheels with it.

Both are **moving no-slip walls**, implemented as the standard bounce-back
correction `f_q = f_q̄ + 2·w_q·ρ·(c_q·u_w)/c_s²`, with the wall velocity stored
per solid cell. Wheels roll without slip, so `ω = u₀/R` about the axle: the
contact patch sits at `+u₀` matching the road and the crown at `−u₀`, meaning
the crown meets the oncoming air at **twice** the free-stream speed. That is
the whole reason rotating wheels matter. Verified directly: contact `+0.85 u₀`,
crown `−0.85 u₀` at 0.85R, hub zero.

The boundary condition is verified by the floor profile in an empty tunnel —
static floor gives `0.21, 0.62, 0.93, 1.07…` (a clear boundary layer), rolling
road gives `1.02, 0.92, 1.00, 1.00…` (none).

**It needs a resolved ride height.** Measured on the same sedan: with a 5.8-cell
underbody gap the rolling road changes C_d by **−2%**, which matches real
tunnel correlation. At the default 1.8-cell gap it changes it by **+242%** —
a one-cell gap cannot support a shear layer. The panel reports the ride height
in cells and warns below four. Wheel rotation also needs identifiable wheel
geometry, so it applies to the bundled vehicles; uploads get the rolling road
only.

## Single-core vs multi-core

A toggle in the Solver section. Double-buffered lattice Boltzmann is
embarrassingly parallel: every worker reads `f` and writes a **disjoint** part
of `g`, so the domain splits into z-slabs with no locking and no halo exchange
— one barrier per step is the entire synchronisation.

The lattice lives in `SharedArrayBuffer`s, so workers operate on it in place and
switching modes costs nothing but spawning threads; the flow state carries
straight over. A batch of steps costs exactly two `postMessage`s per worker,
everything between is `Atomics`.

The kernel is **not duplicated**: workers call the same `stepSlab()` the
single-threaded path calls, so the two cannot drift. Verified — across 120 steps
on 7 workers the velocity and density fields are **bit-identical** to the
single-threaded solver (max difference exactly 0); forces differ by 3.6e-15,
which is float summation order.

Boundaries are z-sliced across workers too, so there is no serial section for
one thread to run alone while the others wait.

**⏱ Benchmark & tune** times 1 core against several worker counts on your
machine and keeps the best. Worth pressing once: the optimum is
machine-specific, since `hardwareConcurrency` counts SMT siblings and this
kernel is memory-bandwidth bound, so running a thread on every logical core can
be *slower* than half of them.

⚠ Multi-core needs `SharedArrayBuffer`, which browsers gate behind
cross-origin isolation. `serve.py` sends the required COOP/COEP headers.
Opening `index.html` straight off disk still works, but stays single-core —
a `file://` URL has no headers at all, and the toggle says so.

## Measurement

Coefficients are reported as a **mean with a 95% confidence interval**, plus a
live trace of C_d and C_l against time and a convergence state
(`averaging` → `noisy` / `drifting` → `converged`). The point is to make it
obvious whether a number has finished moving, instead of reading a digit that
is still drifting.

Error bars come from the **blocking method** (16 blocks) rather than `σ/√N`.
Consecutive timesteps in a shedding wake are heavily correlated, so the naive
formula understates the true error badly; averaging into blocks longer than the
correlation time decorrelates them, and the scatter between block means is
honest.

Samples are **pair-averaged over consecutive timesteps** before being recorded.
The momentum-exchange force carries a strong odd/even parity oscillation —
measured lag-1 autocorrelation of −0.98 — a known lattice Boltzmann artefact
whose mean is right but whose per-step value alternates. Cancelling it leaves
the mean untouched and cuts the scatter ~16× (σ 1.50 → 0.09); without it the
error bar measures the artefact rather than the flow.

**● Start capture** marks a deliberate measurement window. Nothing is recorded
for export until you press it, and it stops the moment you press it again — or
automatically if you change the test conditions, since a capture spanning a
change of body, angle or wind speed would be meaningless.

**⤓ Export** writes a full report: every test condition (grid, τ, ν, LES
setting, both Reynolds numbers, reference area and type, blockage, α, yaw, wind
speed, body length), the results with confidence intervals in both coefficient
and newton form, the model's licence credits if it was an upload, and the raw
C_d / C_l / C_y series. Statistics in the report are computed over the captured
window, not the live one.

With no capture it falls back to exporting the rolling live window, so the
button always produces something; the CSV records which of the two it was in
`sample_source`.

A small rolling buffer always runs regardless — the live coefficients, error
bars and trace chart are computed from it — but that is only the last few
thousand steps, an arbitrary window rather than one you chose.

## Validation

`test.html` runs the standard bluff bodies and prints measured `C_d` / `C_l`
against textbook values. Run it after changing anything in `js/lbm.js`.
Measured on a 96×44×44 lattice at the default Reynolds setting:

| case | C_d sim | C_d full-scale | note |
|---|---|---|---|
| cylinder | 1.25 | 1.17 | ✔ |
| cube | 1.36 | 1.05 | ✔ |
| sphere | 0.96 | 0.47 | matches the *low-Re* value (~1.0 at Re≈150) |
| sedan | 1.11 | 0.32 | ranks below the truck ✔ |
| truck | 1.54 | 0.75 | ✔ |
| teardrop | 1.03 | 0.04 | streamlining does not pay off at Re≈10³ |

Wing lift curve (planform area, NACA 4418):
`C_l = −0.26 / +0.01 / +0.42` at `α = −8° / 0° / +8°` — a slope of ≈2.4 per
radian against ≈1.9 from lifting-line theory at this aspect ratio.

## Honest limitations

* A real car at 100 km/h sits at **Re ≈ 10⁷**; a 96×44×44 browser grid resolves
  **Re ≈ 10³**. Both numbers are shown in the panel. The qualitative physics —
  separation, recirculation, vortex shedding, wake momentum deficit, pressure
  recovery — is right. Absolute `C_d` is *indicative, not certified*: the
  boundary layer is far too thick at this resolution, so bluff-body drag reads
  high. **Compare shapes against each other**, not against a wind-tunnel report.
* Tunnel blockage is kept near 10% of the test-section area, but that is still
  well above real-tunnel practice and inflates `C_d` further.
* Smoke advection is time-scaled for visibility — particle speed on screen is
  not the physical air speed. Streamlines are traced through the instantaneous
  field.

## Layout

```
index.html        markup + panels
test.html         solver validation suite
testmodels.html   mesh pipeline validation suite
serve.py          no-cache dev server
css/style.css     both themes
js/lbm.js         D3Q19 lattice Boltzmann solver (the physics)
js/meshio.js      .obj / .stl / .glb readers -> triangle soup
js/voxel.js       triangle soup -> solid occupancy grid
js/shapes.js      test bodies (analytic + uploaded) -> mesh + point test
js/main.js        scene, visualisation, UI, interaction
vendor/           three.js r147
```
