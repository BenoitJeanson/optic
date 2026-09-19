# optic

[![CI](https://github.com/BenoitJeanson/optic/actions/workflows/ci.yml/badge.svg)](https://github.com/BenoitJeanson/optic/actions/workflows/ci.yml)
[![Demo](https://github.com/BenoitJeanson/optic/actions/workflows/pages.yml/badge.svg)](https://benoitjeanson.github.io/optic/)
[![Licence](https://img.shields.io/badge/licence-Apache--2.0-blue.svg)](LICENSE)

An open-source optical design environment: sequential ray tracing, aberration analysis
and lens optimisation, with a desktop application that a Zemax user can sit down in front
of without retraining.

**Status: early. Milestone 0 in progress.** The kernel traces real and paraxial rays
through centred systems and differentiates them exactly. There is no optimiser and no
desktop application yet. See [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) for where this
is going and why.

## Try it in your browser

The kernel compiles to WebAssembly, so the demo runs entirely on your machine with
nothing installed and nothing sent to a server: **https://benoitjeanson.github.io/optic/**

Edit a radius, drag the defocus slider, switch between the built-in systems, and the
layout, spot diagrams and first-order data update as you type.

## What works today

- Sequential real ray tracing through spherical, conic and even-aspheric surfaces
- Paraxial solver: focal length, back focal distance, entrance and exit pupils, f-number
- Sellmeier, Schott, model-glass and fixed-index dispersion, with a small verified catalog
- **Exact derivatives** of any traced quantity with respect to any design parameter,
  from a single trace, via forward-mode automatic differentiation
- Thickness solves — marginal ray height, chief ray height and pickups — resolved to a
  fixed point before tracing, with derivatives flowing through them
- **Zemax `.zmx` import and export**, including UTF-16 files, six-digit glass codes, and
  an explicit list of anything the file used that we do not yet model
- Zemax-style **vignetting factors** per field — two decentres, two compressions and a
  rotation of the pupil — honoured when rays are launched, carried through `.zmx` in
  both directions, and deliberately ignored by distortion, which stays a property of
  the lens rather than of how the pupil was sampled
- Spot diagrams with three pupil sampling patterns, distortion per field and wavelength,
  scale layout drawings with vignetting shown, and a browser demo

## Try it

Requires Rust 1.78 or newer.

```bash
cargo run --release --example report   # first-order data and spot sizes for the samples
cargo test --workspace                 # the verification suite

./scripts/build-web.sh                 # build the demo
python3 -m http.server -d web 8080     # then open http://localhost:8080
```

## Verification

165 tests, in three layers, plus 39 headless checks of the browser demo. **Unit tests** live beside the code they cover, one module at a
time, and can reach private functions: dual-number calculus, vector and transform algebra,
dispersion formulae, sag and intersection geometry, the paraxial marching step, and the
tracer's refraction, reflection, clipping and path-length bookkeeping. **Integration
tests** in `tests/` exercise whole systems.

Both layers are checked against closed forms and physical invariants, not against the
kernel's own past output:

| Check | Result |
|---|---|
| Singlet focal length vs. the thick-lens equation | agrees to 1e-9 |
| Catalog `n_d` and `V_d` vs. published values | 5 glasses, within 1e-5 and 0.01 |
| Paraxial image height vs. `EFL·tan θ` | agrees to 1e-9 |
| Snell's law and coplanarity at every interface | residual < 1e-12 |
| Lagrange invariant across the system | conserved to 1e-10 relative |
| Real vs. paraxial rays as aperture and field shrink | converges at exactly third order |
| Autodiff gradients vs. central differences | agree to 1e-5 relative |
| Autofocus solve vs. the computed back focal distance | agree to 1e-10 |
| **EP 155,640 (1919) triplet** vs. its published focal length | within the source's rounding |
| **DE 287,089 (1913) triplet** vs. its published focal length | within the source's rounding |
| Distortion of both, vs. their published plots | right sign and magnitude |
| `.zmx` round trip | focal length preserved exactly |
| Derivatives through a solved thickness | agree with central differences to 1e-5 |
| Thin lens vs. the lensmaker's equation | 4 configurations, to 1e-9 |
| Concave mirror focal length vs. `R/2` | to 1e-10 |
| Plane-parallel plate displacement vs. `t(tan A - tan A')` | to 1e-12 |
| A chief ray aimed at the entrance pupil | crosses the axis at the stop, to 1e-12 |
| Paraxial trace reversed through `backward` | returns the launch state, to 1e-12 |

Two published triplets are reproduced to within the precision of the source. Because the
prescriptions are printed to one decimal place, agreement is judged against a budget
computed from the data itself: the focal length is differentiated with respect to all
eleven rounded quantities at once, and the disagreement must fit inside what ±0.05 on
each can produce. That is the differentiable kernel earning its keep on a question it
was not built for.

The third-order convergence result is the load-bearing one: halving both aperture and
field divides the real-versus-paraxial disagreement by eight, which is what aberration
theory demands and which nothing but a correct tracer will produce.

The demo has its own headless smoke test (`tests/web`), which runs the real `app.js`
against a real DOM. It catches what unit tests cannot: a mistyped element id, a listener
on the wrong event, a canvas that never gets drawn. The demo does not deploy unless it
passes.

CI runs everything on Linux, macOS and Windows, with and without default features, plus
formatting, clippy, a minimum-supported-Rust-version check, and a WebAssembly build.

## Layout

```
crates/optic-core    the kernel: geometry, materials, tracing, paraxial optics, solves
crates/optic-io      reading and writing prescriptions: Zemax .zmx in and out
crates/optic-wasm    a JSON analysis API over it, and the C ABI the browser calls
web/                 the demo page -- plain ES modules, no build step
tests/web            headless smoke test for the page
docs/                architecture and decisions
```

## Licence

Apache-2.0. See [LICENSE](LICENSE).
