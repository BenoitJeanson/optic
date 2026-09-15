# optic

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
- Spot diagrams, scale layout drawings with vignetting shown, and a browser demo

## Try it

```bash
cargo run --release --example report   # first-order data and spot sizes for the samples
cargo test --workspace                 # the verification suite

./scripts/build-web.sh                 # build the demo
python3 -m http.server -d web 8080     # then open http://localhost:8080
```

## Verification

125 tests, in two layers. **Unit tests** live beside the code they cover, one module at a
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
| Thin lens vs. the lensmaker's equation | 4 configurations, to 1e-9 |
| Concave mirror focal length vs. `R/2` | to 1e-10 |
| Plane-parallel plate displacement vs. `t(tan A - tan A')` | to 1e-12 |
| A chief ray aimed at the entrance pupil | crosses the axis at the stop, to 1e-12 |
| Paraxial trace reversed through `backward` | returns the launch state, to 1e-12 |

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
crates/optic-core    the kernel: geometry, materials, tracing, paraxial optics
crates/optic-wasm    a JSON analysis API over it, and the C ABI the browser calls
web/                 the demo page -- plain ES modules, no build step
tests/web            headless smoke test for the page
docs/                architecture and decisions
```

## Licence

Apache-2.0. See [LICENSE](LICENSE).
