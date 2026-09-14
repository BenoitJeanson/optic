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

The kernel is checked against closed forms and physical invariants, not against itself:

| Check | Result |
|---|---|
| Singlet focal length vs. the thick-lens equation | agrees to 1e-9 |
| Catalog `n_d` and `V_d` vs. published values | 5 glasses, within 1e-5 and 0.01 |
| Paraxial image height vs. `EFL·tan θ` | agrees to 1e-9 |
| Snell's law and coplanarity at every interface | residual < 1e-12 |
| Lagrange invariant across the system | conserved to 1e-10 relative |
| Real vs. paraxial rays as aperture and field shrink | converges at exactly third order |
| Autodiff gradients vs. central differences | agree to 1e-5 relative |

The third-order convergence result is the load-bearing one: halving both aperture and
field divides the real-versus-paraxial disagreement by eight, which is what aberration
theory demands and which nothing but a correct tracer will produce.

## Licence

Apache-2.0. See [LICENSE](LICENSE).
