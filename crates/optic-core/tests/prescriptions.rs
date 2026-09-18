//! Reference designs from the literature.
//!
//! Every other test in this project checks the kernel against itself or against a closed
//! form. These check it against lenses somebody else designed, computed and published —
//! the only evidence that carries weight with a working optical designer.
//!
//! The two triplets here come from expired patents, tabulated in W. Smith, *Modern
//! Optical Engineering*, figures 12.13 and 12.14. Both are published at focal length
//! 100 units, with plots of spherical aberration, field curvature and distortion.
//!
//! **Comparing against them needs care about precision.** The prescriptions are quoted
//! to one decimal place, so each radius and thickness carries up to ±0.05 of unstated
//! value. Rather than assert a tolerance picked by eye, each test differentiates the
//! focal length with respect to all eleven rounded quantities and asks whether the
//! disagreement fits inside what that rounding can produce. Exact derivatives are what
//! make that practical: one trace gives the whole gradient.

use optic_core::surface::Profile;
use optic_core::{lines, samples, trace, Dual, Field, Object, Paraxial, Scalar, Surface, System};

/// Six radii and five thicknesses: everything the source rounds.
type D = Dual<11>;

const ROUNDING: f64 = 0.05;

/// Rebuild a system with every rounded quantity seeded as a variable.
fn with_rounding_variables(sys: &System<f64>) -> System<D> {
    let mut surfaces = Vec::with_capacity(sys.surfaces.len());
    for (i, s) in sys.surfaces.iter().enumerate() {
        let profile = match s.profile {
            Profile::Plane => Profile::Plane,
            Profile::Conic { curvature, conic } => {
                // The source rounds radii, not curvatures, so the radius is the variable.
                let radius = D::variable(1.0 / curvature, i.min(5));
                Profile::Conic {
                    curvature: D::from_f64(1.0) / radius,
                    conic: D::from_f64(conic),
                }
            }
            _ => unreachable!("these designs are all-spherical"),
        };
        surfaces.push(Surface {
            profile,
            thickness: if i < 5 {
                D::variable(s.thickness, 6 + i)
            } else {
                D::from_f64(s.thickness)
            },
            thickness_solve: Default::default(),
            material: s.material.clone(),
            semi_diameter: s.semi_diameter.map(D::from_f64),
            is_stop: s.is_stop,
            label: s.label.clone(),
        });
    }
    let mut out = System::new(Object::Infinity, surfaces);
    out.aperture = sys.aperture;
    out.fields = sys.fields.clone();
    out.wavelengths = sys.wavelengths.clone();
    out
}

/// Focal-length uncertainty implied by the published rounding: RMS, and worst case.
fn focal_length_budget(sys: &System<f64>) -> (f64, f64) {
    let efl = Paraxial::compute(&with_rounding_variables(sys), lines::D).efl;
    let grad = efl.grad();
    let rms = grad
        .iter()
        .map(|g| (g * ROUNDING).powi(2))
        .sum::<f64>()
        .sqrt();
    let worst = grad.iter().map(|g| g.abs() * ROUNDING).sum();
    (rms, worst)
}

fn distortion_at(sys: &System<f64>, field: Field) -> f64 {
    let par = Paraxial::compute(sys, lines::D);
    let paraxial = par.image_height(sys, lines::D, field);
    let point = trace(
        sys,
        lines::D,
        optic_core::launch(sys, &par, field, 0.0, 0.0),
    )
    .image_point()
    .expect("the chief ray reaches the image");
    let real = (point.x * point.x + point.y * point.y).sqrt();
    100.0 * (real - paraxial) / paraxial
}

#[test]
fn published_triplets_reproduce_their_focal_length() {
    for sys in [
        samples::smith_triplet_moderate::<f64>(),
        samples::smith_triplet_wide::<f64>(),
    ] {
        let efl = Paraxial::compute(&sys, lines::D).efl;
        let error = (efl - 100.0).abs();
        let (rms, worst) = focal_length_budget(&sys);

        assert!(
            error <= worst,
            "{}: focal length {efl:.4} differs from the published 100 by {error:.4}, \
             more than the published rounding can explain (worst case {worst:.4})",
            sys.title
        );
        // A sanity floor: if the budget ever collapsed, the test above would pass
        // vacuously.
        assert!(
            rms > 0.0 && worst.is_finite() && worst < 10.0,
            "{}",
            sys.title
        );
    }
}

#[test]
fn the_rounding_budget_explains_why_one_design_agrees_better_than_the_other() {
    // The wide-field design has much smaller radii, so the same +/-0.05 of rounding
    // moves its focal length far more. If this relationship ever inverted, the error
    // would not be coming from the source data.
    let (moderate_rms, _) = focal_length_budget(&samples::smith_triplet_moderate::<f64>());
    let (wide_rms, _) = focal_length_budget(&samples::smith_triplet_wide::<f64>());
    assert!(
        wide_rms > moderate_rms * 2.0,
        "expected the short-radius design to be far more sensitive: {moderate_rms:.3} vs {wide_rms:.3}"
    );
}

#[test]
fn published_triplets_show_the_distortion_their_plots_show() {
    // Both designs are plotted on a +/-1% distortion scale: the moderate one essentially
    // flat out to 20 degrees, the wide-field one climbing towards 1% at 30 degrees. The
    // stop position is not printed, and distortion depends on it, so these are ranges
    // rather than measurements.
    let moderate = samples::smith_triplet_moderate::<f64>();
    let d20 = distortion_at(&moderate, Field::angle(20.0));
    assert!(
        d20.abs() < 0.25,
        "EP 155,640 is plotted as nearly distortion-free at 20 degrees, got {d20:+.3}%"
    );

    let wide = samples::smith_triplet_wide::<f64>();
    let d30 = distortion_at(&wide, Field::angle(30.0));
    assert!(
        (0.2..1.5).contains(&d30),
        "DE 287,089 is plotted at roughly +1% at 30 degrees, got {d30:+.3}%"
    );
    assert!(d30 > d20, "the wide-field design should distort more");
}

#[test]
fn published_triplets_trace_over_their_full_field() {
    for sys in [
        samples::smith_triplet_moderate::<f64>(),
        samples::smith_triplet_wide::<f64>(),
    ] {
        let par = Paraxial::compute(&sys, lines::D);
        assert!(
            par.efl > 0.0 && par.fno > 1.0 && par.fno < 10.0,
            "{}",
            sys.title
        );

        for field in sys.fields.clone() {
            for (px, py) in [(0.0, 0.0), (0.0, 0.95), (0.0, -0.95), (0.95, 0.0)] {
                let traced = trace(
                    &sys,
                    lines::D,
                    optic_core::launch(&sys, &par, field, px, py),
                );
                assert!(
                    traced.is_complete(),
                    "{}: ray ({px}, {py}) at {:.0} degrees was blocked: {:?}",
                    sys.title,
                    field.radius(),
                    traced.blocked
                );
            }
        }
    }
}

#[test]
fn model_glass_from_a_six_digit_code_reproduces_the_code() {
    // The designs are specified by glass code, so the code must round-trip exactly or
    // every conclusion drawn from them is built on the wrong glass.
    for (code, nd, vd) in [
        ("613585", 1.613, 58.5),
        ("621362", 1.621, 36.2),
        ("549458", 1.549, 45.8),
    ] {
        let glass = optic_core::glass_code(code).expect("valid code");
        assert!((glass.nd() - nd).abs() < 1e-9, "{code}: n_d {}", glass.nd());
        assert!(
            (glass.abbe() - vd).abs() < 1e-6,
            "{code}: V_d {}",
            glass.abbe()
        );
    }
    assert!(
        optic_core::glass_code("N-BK7").is_none(),
        "a catalogue name is not a code"
    );
    assert!(
        optic_core::glass_code("61358").is_none(),
        "five digits is not a code"
    );
    assert!(
        optic_core::glass_code("6135850").is_none(),
        "seven digits is not a code"
    );
    assert!(
        optic_core::glass_code("000000").is_none(),
        "an index of 1.0 is not a glass"
    );
    // The decoder checks that each number is in range, not that the pair is achievable:
    // no real glass has n_d 1.999 at V_d 99.9, but encoding the index-dispersion
    // trade-off is a catalogue's job, not a decoder's.
    assert!(optic_core::glass_code("999999").is_some());
}
