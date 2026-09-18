# Architecture

This document explains what we are building, the decisions that are already made, and
why. It is written for anyone joining the project — no Rust required to follow the
reasoning, though the code examples are Rust.

## The goal

Zemax OpticStudio is five products in one: a sequential design engine, an analysis suite,
an optimiser with tolerancing, a non-sequential engine for illumination and stray light,
and a physical-optics propagator. The open-source world covers fragments of this well —
`rayoptics` and `Optiland` for sequential design, `prysm` for physical optics, `POPPY` for
astronomical diffraction — and covers the *integration* not at all.

The gap is not the physics. The physics is in textbooks. The gap is the application: a
tool where changing a radius updates the spot diagram, the MTF and the layout at once, and
where an optimiser is one click away. That is what we are building, starting with imaging
lens design.

## Decisions already made

| Decision | Choice |
|---|---|
| First user | Sequential imaging design — camera lenses, objectives, telescopes |
| Kernel | Rust, scalar-generic and differentiable |
| Scripting | Python via PyO3, exposing the same API the UI uses |
| Application | Tauri v2 desktop shell, web-tech UI |
| UX stance | Familiar Zemax layout, modernised behaviour |
| Zemax `.zmx` import | Early — it is the adoption lever |
| Licence | Apache-2.0 |

## Decisions from Elias

Answers to the opening questions, recorded here because they are binding.

| Question | Answer |
|---|---|
| Reference designs | **W. Smith, *Modern Lens Design*** — the source his practical work used. Malacara's *Handbook of Optical Design*, Kingslake's *Lens Design Fundamentals* and Geary's *Introduction to Lens Design* give analytical treatments of the Cooke triplet worth checking against |
| Agreement tolerance | **±0.5%** against published values |
| Pupil model | **Vignetting factors, as in Zemax** — not physical apertures |
| Sign conventions | **Zemax's**, throughout |
| Aspheres | Needs the **r² term**, and **Forbes Q-type** kept |
| Glasses | **Schott** mainly, other vendors useful, **no obsolete glasses**. **Temperature-dependent index is required**, not optional |
| Merit function | Wizard-built base (usually RMS spot) **plus hand-added operands**. Wants the full Zemax operand set, and is open to new ones — he suggested wavefront quality on a surface other than the image |
| File format | **`.zmx` only** |

The analyses he actually uses, which is the M1 list: layout (2D and 3D), spot diagrams,
MTF, Seidel diagram and coefficients, ray fan, footprint, field curvature and distortion,
grid distortion, longitudinal aberration, lateral colour, chromatic focal shift, wavefront
map, interferogram, PSF.

### First external validation

Two Cooke triplets from expired patents — English Patent 155,640 (1919) and German Patent
287,089 (1913), tabulated in W. Smith, *Modern Optical Engineering*, figures 12.13 and
12.14 — are now reference designs in the test suite. Both are published at focal length
100 units with plots of spherical aberration, field curvature and distortion.

| | our focal length | published | rounding budget (RMS) |
|---|---|---|---|
| EP 155,640 | 99.637 | 100 | ±0.276 |
| DE 287,089 | 99.089 | 100 | ±1.534 |

Both agree to within the precision of the printed data. That claim is not a judgement
call: the focal length is differentiated with respect to all six radii and five
thicknesses at once, and the ±0.05 implied by one-decimal rounding is propagated through.
The wide-field design's shorter radii make it eleven times more sensitive — its first
radius moves the focal length by 20.8 units per unit of radius, against 3.6 for the other
— which is exactly why it agrees less well. A discrepancy that tracks sensitivity is
rounding; one that does not would be a bug.

Distortion agrees too: +0.030% at 20° and +0.712% at 30°, against plots drawn on a ±1%
scale showing very little and roughly one percent respectively.

Note that this is *Modern Optical Engineering*, not the *Modern Lens Design* Elias named.
Same author, different book: the former is a textbook with a handful of worked designs,
the latter a catalogue of about a hundred with performance data. The catalogue remains
the better source.

### The open discrepancy

Elias analysed a 50 mm f/5 triplet and got RMS spot radii of 13.819 µm on axis and
32.843 µm at 20°, with distortion of +0.1389% at the primary line. This kernel gives
13.5 µm, 23.0 µm and +0.1153%.

The 20° difference is **not** explained by any setting. Removing every aperture so nothing
is vignetted moves it by 1.4%; sampling is converged from 531 rays to 15,491; using one
wavelength instead of three moves it the wrong way. Nothing reaches 32.8 µm.

Distortion settles it. It is a **chief-ray** property — one ray through the centre of the
pupil — so it cannot depend on aperture, vignetting factors, ray aiming or pupil sampling.
Two tools that disagree about distortion are not disagreeing about settings; they are
describing different lenses. Since the sample prescription here was constructed to be
physically sound rather than transcribed from a publication, the likeliest explanation by
far is that it is not the same triplet.

**Resolving this needs his prescription**, not more analysis. Which is also the argument
for pulling `.zmx` import forward: every file he has becomes an exact comparison.

## Layering

```
crates/optic-core        geometry, surfaces, materials, sequential trace   <- we are here
crates/optic-analysis    spot, ray fans, OPD, Zernike, PSF, MTF
crates/optic-opt         merit operands, Levenberg-Marquardt, tolerancing
crates/optic-io          document format, .zmx import, catalogs, CAD export
crates/optic-py          PyO3 bindings
crates/optic-rpc         the command and event surface the UI speaks
app/                     Tauri + React/TypeScript
```

`optic-core` contains no I/O, no formatting and no UI. Every analysis is a pure function
of a `System`. That is what makes results cacheable, optimisation reproducible, and the
whole kernel testable without a running application.

The UI is a *client* of the same command surface Python uses. If you can click it, you can
script it — enforced by construction rather than by discipline.

## The three load-bearing decisions

### 1. Forward-mode automatic differentiation

Every geometric routine is generic over a `Scalar` trait. Instantiate it with `f64` and
you get an ordinary ray trace. Instantiate it with `Dual<N>` and you get the trace *plus*
exact derivatives with respect to `N` design variables, from the same pass.

```rust
let sys = cooke_with(Dual::<2>::variable(0.0, 0),   // curvature of surface 1
                     Dual::<2>::variable(0.0, 1));  // thickness of surface 2
let y = trace(&sys, wl, ray).image_point().unwrap().y;
let [dy_dc, dy_dt] = *y.grad();                      // exact, not finite-differenced
```

Zemax computes merit-function Jacobians by finite differences, which forces a step-size
compromise between truncation and rounding error — one that is hardest to get right for
high-order aspheric coefficients. We get exact derivatives instead, with no step to tune.
Both approaches scale as `O(n)` in the number of variables; what dual numbers save is a
constant factor, because control flow, the value part and shared work are done once in a
single pass rather than `n + 1` times. The result is that damped least squares — which
*is* Levenberg–Marquardt — gets a cleaner, somewhat cheaper Jacobian, and the same
machinery hands us tolerance sensitivities.

The price is that this only works because every routine on the path from a variable to a
residual is generic over `Scalar`. A finite-difference optimizer can treat the merit
function as a black box — user DLLs, macros, arbitrary operands — and we cannot. Any
iterative solve or non-smooth step (ray aiming, clipping, vignetting) also has to be
differentiated deliberately rather than inherited for free.

Forward mode rather than reverse is deliberate. Lens design has few variables (tens:
curvatures, thicknesses, conics, aspheric coefficients) and many residuals (thousands: ray
errors across pupil, field and wavelength). Levenberg–Marquardt needs the whole `m × n`
Jacobian; forward mode produces all of it in `~n` trace-passes, where reverse mode would
cost `~m`.

One subtlety already handled: aspheric intersection needs Newton iteration, and naively
differentiating through the iterations makes the gradient depend on the iteration count.
Instead, once the value has converged we reset `t` to a constant and take a single Newton
step in dual arithmetic, which reproduces `dt/dp = −(∂F/∂p)/(∂F/∂t)` exactly. There is a
test asserting that the iteration count cannot leak into the Jacobian.

### 2. The system is a document

An optical system is a plain, serialisable value — text on disk, git-diffable. Everything
else is a pure function of it. This gives undo/redo, version control, headless CI,
reproducible papers and diffable design reviews essentially for free, and it is a real
differentiator: the opacity of Zemax's file format is a standing complaint.

Two departures from the Zemax data model, both deliberate:

- **The object is not row 0 of the surface list.** An object at infinity would put
  `INFINITY` into the scalar type, and infinity times a zero derivative is `NaN` — one
  poisoned entry propagates through the whole Jacobian. Modelling the object as its own
  sum type keeps every number in the trace finite. The editor can still *present* it as
  row 0.
- **A surface owns the medium that follows it**, which is how prescriptions are written
  and read.

### 3. Validation is the product

Nobody adopts a new optical design tool on charm. The kernel is checked against closed
forms and physical invariants, never against its own past output alone:

- the thick-lens equation, analytically, for a singlet
- published `n_d` and `V_d` for every catalog glass, so a transcription error fails the
  build rather than silently biasing every design
- Snell's law and incident/refracted/normal coplanarity at every interface
- conservation of the Lagrange invariant through the system
- **third-order convergence**: halving aperture and field together divides the
  real-versus-paraxial disagreement by exactly eight
- autodiff gradients against central differences

Tests that merely pin today's numbers are labelled baselines and kept separate from these.
They are not claims about published designs.

## What is deliberately not built yet, but already has room

- **Tilts and decentres.** The tracer already visits every surface through a `Transform`,
  so coordinate breaks compose into the existing structure without touching the trace.
  This is the single most invasive feature to retrofit, hence the early indirection.
- **Mirrors.** Signed indices track reflection parity through the paraxial trace. Untested;
  no sample uses one yet.
- **Solves.** Built: marginal-ray-height, chief-ray-height and pickups on thicknesses,
  evaluated to a fixed point before any ray is traced, with derivatives flowing through
  them so the optimiser sees a solved thickness's true sensitivity. Still to come:
  solves on curvature and on glass, and f-number solves.
- **Real ray aiming.** Rays are currently aimed at the *paraxial* entrance pupil, so a
  system with strong pupil aberration will not have its pupil filled uniformly. The
  interface does not change when real aiming arrives.
- **Glass as a continuous variable.** Real glass is discrete. The standard approach is to
  optimise in `(n_d, V_d)` space with a penalty keeping the design inside the glass map,
  then substitute the nearest real glasses and re-optimise. `Material::ModelGlass` exists
  for this.

## Catalogs and licensing

Manufacturer `.agf` files are freely downloadable but their redistribution terms are
unclear. The repository therefore ships no vendor catalog. The built-in glasses are a
handful of Schott N-series entries whose coefficients are verified against published
constants by the test suite; anything larger will come from a fetcher plus a parser, or
from refractiveindex.info, which is CC0.

## Roadmap

- **M0** — sequential trace, conic and aspheric surfaces, catalog glasses, spot diagrams,
  ray fans, layout. CLI and Python only. *In progress.*
- **M1** — optimisation: variables, merit operands (the Zemax set), Levenberg–Marquardt
  with analytic Jacobians. Vignetting factors. The analysis list above, starting with
  Seidel, MTF, ray fans and field curvature. Aspheric r² term and Forbes Q-type.
- **M2** — the desktop application, with live-updating analysis windows. `.zmx` import
  and export are **done** and available in the browser demo, so any Zemax file can be
  opened here and any design here checked in Zemax.
- **M3** — tolerancing (sensitivity and Monte Carlo), coatings, polarisation. Thermal
  moves earlier if temperature work blocks him.
- **M4** — non-sequential mode: sources, detectors, scattering, stray light.
- **M5** — physical optics propagation, Gaussian beams, CAD and ISO 10110 export.

M0 and M1 hold the intellectual work. M2 is where it becomes software people use.

## Working on this

The kernel's correctness is the whole project's credibility, so the physics and the
implementation are reviewed separately.

**Physics and validation** — `tests/`, `samples.rs`, `material.rs`, and the prescription
suite. Reviewing these needs optics, not Rust: the question is always "is this the right
equation, the right sign convention, the right published value?" New reference designs
with literature values are the highest-value contribution available right now.

**Kernel implementation** — `math/`, `surface.rs`, `trace.rs`, `paraxial.rs`.

Run `cargo test` before pushing; `cargo clippy --all-targets` is expected to be clean.
