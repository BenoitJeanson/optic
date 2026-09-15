# optic

An open-source optical design environment: sequential ray tracing, aberration analysis
and lens optimisation, with a desktop application that a Zemax user can sit down in front
of without retraining.

**Status: early. Milestone 0 in progress.** The kernel traces real and paraxial rays
through centred systems and differentiates them exactly. There is no optimiser and no UI
yet. See [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) for where this is going and why.

## What works today

- Sequential real ray tracing through spherical, conic and even-aspheric surfaces
- Paraxial solver: focal length, back focal distance, entrance and exit pupils, f-number
- Sellmeier, Schott, model-glass and fixed-index dispersion, with a small verified catalog
- **Exact derivatives** of any traced quantity with respect to any design parameter,
  from a single trace, via forward-mode automatic differentiation

## Try it

```bash
cargo run --release --example report   # first-order data and spot sizes for the samples
cargo test                             # the verification suite
```

## Verification

104 tests, in two layers. **Unit tests** live beside the code they cover, one module at a
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

CI runs the suite on Linux, macOS and Windows, with and without default features, plus
formatting, clippy, a minimum-supported-Rust-version check, and a WebAssembly build.

## Licence

Apache-2.0. See [LICENSE](LICENSE).
