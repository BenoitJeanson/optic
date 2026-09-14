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

Zemax computes merit-function Jacobians by finite differences, which costs one extra trace
per variable and forces a step-size compromise between truncation and rounding error. We
get exact derivatives instead. This makes damped least squares — which *is*
Levenberg–Marquardt — faster and more robust, and hands us tolerance sensitivities for
free.

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
- **Solves and pickups.** Marginal-ray-height solves, f-number solves, parameter pickups.
  These turn the document into a small dependency graph that must be evaluated before
  every trace. Retrofitting that is painful, so it goes into the data model early.
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
- **M1** — optimisation: variables, solves, merit operands, Levenberg–Marquardt with
  analytic Jacobians. Zernike and Seidel coefficients, OPD, MTF.
- **M2** — the desktop application, with live-updating analysis windows. `.zmx` import.
- **M3** — tolerancing (sensitivity and Monte Carlo), coatings, polarisation, thermal.
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
