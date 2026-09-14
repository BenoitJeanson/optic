//! Kernel verification.
//!
//! Two kinds of test live here, and the distinction matters.
//!
//! * **Invariants and closed forms** are proofs: the thick-lens equation, Snell's law at
//!   every interface, conservation of the Lagrange invariant, agreement between real and
//!   paraxial rays as the aperture shrinks, derivatives against central differences.
//!   These are right or the kernel is wrong.
//! * **Baselines** merely pin today's numbers so a refactor cannot move them silently.
//!   They are not claims that the numbers match any published design; validating the
//!   sample prescriptions against the literature is separate work.

use optic_core::{
    catalog, launch, lines, samples, surface::Profile, trace, Dual, Field, Material, Paraxial,
    Scalar, System, Vec3,
};

const TOL: f64 = 1e-9;

/// Build the Cooke triplet with two of its parameters offset, generic over the scalar
/// type so the same construction serves both the dual-number and finite-difference paths.
fn cooke_with<S: Scalar>(d_curvature: S, d_thickness: S) -> System<S> {
    let mut sys = samples::cooke_triplet::<S>();
    if let Profile::Conic { curvature, .. } = &mut sys.surfaces[0].profile {
        *curvature += d_curvature;
    }
    sys.surfaces[1].thickness += d_thickness;
    sys
}

// ---------------------------------------------------------------- closed forms

#[test]
fn catalog_reproduces_published_constants() {
    for (name, glass, nd, vd) in catalog::PUBLISHED {
        assert!(
            (glass.nd() - nd).abs() < 1e-5,
            "{name}: n_d {} differs from published {nd}",
            glass.nd()
        );
        assert!(
            (glass.abbe() - vd).abs() < 0.01,
            "{name}: V_d {} differs from published {vd}",
            glass.abbe()
        );
    }
}

#[test]
fn singlet_focal_length_matches_thick_lens_equation() {
    let sys = samples::singlet::<f64>();
    let par = Paraxial::compute(&sys, lines::D);

    let n = catalog::N_BK7.index(lines::D);
    let (r1, r2, d) = (100.0, -100.0, 5.0);
    let expected = 1.0 / ((n - 1.0) * (1.0 / r1 - 1.0 / r2 + (n - 1.0) * d / (n * r1 * r2)));

    assert!(
        (par.efl - expected).abs() < 1e-9,
        "EFL {} vs thick-lens {expected}",
        par.efl
    );
}

#[test]
fn paraxial_image_height_follows_efl_tangent() {
    let sys = samples::cooke_triplet::<f64>();
    let par = Paraxial::compute(&sys, lines::D);
    for deg in [0.0, 5.0, 14.0, 20.0] {
        let field = Field::angle(deg);
        let h = par.image_height(&sys, lines::D, field);
        let expected = par.efl * deg.to_radians().tan();
        assert!(
            (h - expected).abs() < 1e-9,
            "field {deg}: image height {h} vs EFL*tan {expected}"
        );
    }
}

#[test]
fn refocusing_is_idempotent() {
    // The back focal distance must be measured from the last powered surface. If it is
    // measured from the image plane instead, refocusing oscillates about the focus
    // rather than landing on it -- and the second application hides the first's error.
    let mut sys = samples::cooke_triplet::<f64>();
    let once = sys.surfaces[sys.image_index() - 1].thickness;
    samples::focus(&mut sys, lines::D);
    let twice = sys.surfaces[sys.image_index() - 1].thickness;
    assert!((once - twice).abs() < TOL, "{once} then {twice}");
}

// ------------------------------------------------------------------ invariants

#[test]
fn snells_law_holds_at_every_surface() {
    let sys = samples::cooke_triplet::<f64>();
    let par = Paraxial::compute(&sys, lines::D);
    let wl = lines::D;

    for field in [Field::angle(0.0), Field::angle(20.0)] {
        for (px, py) in [(0.0, 0.0), (0.7, 0.0), (0.0, -0.9), (0.5, 0.5)] {
            let start = launch(&sys, &par, field, px, py);
            let traced = trace(&sys, wl, start);
            assert!(traced.is_complete(), "ray blocked: {:?}", traced.blocked);

            let mut incoming = start.dir;
            for hit in &traced.hits {
                let i = hit.surface;
                let normal = sys.surfaces[i].profile.normal(hit.local).unwrap();
                let n1 = sys.medium_before(i).index(wl);
                let n2 = if sys.surfaces[i].material.is_reflective() {
                    n1
                } else {
                    sys.medium_after(i).index(wl)
                };

                let sin1 = incoming.cross(normal).norm();
                let sin2 = hit.dir.cross(normal).norm();
                assert!(
                    (n1 * sin1 - n2 * sin2).abs() < 1e-12,
                    "surface {i}: n1 sin1 = {} vs n2 sin2 = {}",
                    n1 * sin1,
                    n2 * sin2
                );

                // Incident ray, refracted ray and normal must be coplanar.
                let triple = incoming.cross(hit.dir).dot(normal);
                assert!(triple.abs() < 1e-12, "surface {i}: non-coplanar, {triple}");

                incoming = hit.dir;
            }
        }
    }
}

#[test]
fn lagrange_invariant_is_conserved() {
    let sys = samples::cooke_triplet::<f64>();
    let par = Paraxial::compute(&sys, lines::D);
    let chief = optic_core::paraxial::chief_ray(&sys, lines::D, Field::angle(20.0), par.ep_z);

    let mut reference: Option<f64> = None;
    for (m, c) in par.marginal.iter().zip(chief.iter()) {
        let h = m.n * (c.y * m.u - m.y * c.u);
        match reference {
            None => reference = Some(h),
            Some(h0) => assert!(
                (h - h0).abs() < 1e-10 * h0.abs().max(1.0),
                "Lagrange invariant drifted: {h0} -> {h}"
            ),
        }
    }
    assert!(reference.unwrap().abs() > 1e-6, "invariant is degenerate");
}

#[test]
fn real_rays_converge_on_paraxial_at_third_order() {
    // Shrinking the aperture alone does not make a ray paraxial: the chief ray still
    // crosses every surface far off-axis, and what remains is distortion. The paraxial
    // limit is approached only when aperture *and* field shrink together -- and then
    // aberration theory says the residual is third order, so halving the scale must
    // divide the error by eight. Anything else means the real trace and the paraxial
    // trace disagree about first-order optics.
    let mut sys = samples::cooke_triplet::<f64>();
    let mut previous = f64::NAN;
    let mut ratios = Vec::new();

    for k in 0..6 {
        let eps = 0.5f64.powi(k);
        sys.aperture = optic_core::Aperture::EntrancePupilDiameter(10.0 * eps);
        let field = Field::angle(20.0 * eps);
        let par = Paraxial::compute(&sys, lines::D);
        let real = trace(&sys, lines::D, launch(&sys, &par, field, 1.0, 0.0))
            .image_point()
            .expect("marginal ray reached the image")
            .y;
        let error = (real - par.image_height(&sys, lines::D, field)).abs();
        if k > 0 {
            ratios.push(error / previous);
        }
        previous = error;
    }

    for &r in &ratios[2..] {
        assert!(
            (r - 0.125).abs() < 0.015,
            "error ratio {r} is not the 1/8 of third-order convergence; ratios {ratios:?}"
        );
    }
}

#[test]
fn distortion_is_a_chief_ray_property() {
    // Distortion is defined at the chief ray, so it must not depend on aperture.
    let mut sys = samples::cooke_triplet::<f64>();
    let field = Field::angle(20.0);
    let mut values = Vec::new();

    for epd in [1.0, 1e-2, 1e-4] {
        sys.aperture = optic_core::Aperture::EntrancePupilDiameter(epd);
        let par = Paraxial::compute(&sys, lines::D);
        let real = trace(&sys, lines::D, launch(&sys, &par, field, 1.0, 0.0))
            .image_point()
            .unwrap()
            .y;
        let paraxial = par.image_height(&sys, lines::D, field);
        values.push(100.0 * (real - paraxial) / paraxial);
    }

    let spread = values.iter().cloned().fold(f64::MIN, f64::max)
        - values.iter().cloned().fold(f64::MAX, f64::min);
    assert!(spread < 1e-3, "distortion varies with aperture: {values:?}");
    // Baseline, recorded 2026-09-14: a well-behaved triplet, a fraction of a percent.
    assert!(
        (values[2] - 0.1153).abs() < 1e-3,
        "distortion {} vs baseline 0.1153%",
        values[2]
    );
}

#[test]
fn aspheric_intersections_lie_on_the_surface() {
    let profile: Profile<f64> = Profile::EvenAsphere {
        curvature: 1.0 / 30.0,
        conic: -0.7,
        coeffs: vec![1.2e-6, -3.4e-9, 7.0e-12],
    };
    for (ox, oy, dy) in [(0.0, 0.0, 0.0), (3.0, -2.0, 0.05), (-6.0, 6.0, -0.12)] {
        let origin = Vec3::new(ox, oy, -20.0);
        let dir = Vec3::new(0.02, dy, 1.0).normalized();
        let t = profile.intersect(origin, dir).expect("intersection exists");
        let p = origin + dir * t;
        let residual = p.z - profile.sag(p.x * p.x + p.y * p.y).unwrap();
        assert!(residual.abs() < 1e-12, "off-surface by {residual}");
    }
}

// ------------------------------------------------------- automatic differentiation

#[test]
fn dual_gradients_match_central_differences() {
    let wl = lines::D;
    let field = Field::angle(14.0);
    let (px, py) = (0.6, -0.4);

    // One trace in dual arithmetic yields both partials at once.
    let sys_d = cooke_with(Dual::<2>::variable(0.0, 0), Dual::<2>::variable(0.0, 1));
    let par_d = Paraxial::compute(&sys_d, wl);
    let traced = trace(&sys_d, wl, launch(&sys_d, &par_d, field, px, py));
    let y = traced.image_point().expect("ray reached the image").y;
    let [d_curvature, d_thickness] = *y.grad();

    let sample = |dc: f64, dt: f64| -> f64 {
        let sys = cooke_with(dc, dt);
        let par = Paraxial::compute(&sys, wl);
        trace(&sys, wl, launch(&sys, &par, field, px, py))
            .image_point()
            .expect("ray reached the image")
            .y
    };

    let hc = 1e-7;
    let fd_curvature = (sample(hc, 0.0) - sample(-hc, 0.0)) / (2.0 * hc);
    let ht = 1e-6;
    let fd_thickness = (sample(0.0, ht) - sample(0.0, -ht)) / (2.0 * ht);

    let rel = |a: f64, b: f64| (a - b).abs() / a.abs().max(b.abs()).max(1e-12);
    assert!(
        rel(d_curvature, fd_curvature) < 1e-5,
        "d(y)/d(curvature): AD {d_curvature} vs FD {fd_curvature}"
    );
    assert!(
        rel(d_thickness, fd_thickness) < 1e-5,
        "d(y)/d(thickness): AD {d_thickness} vs FD {fd_thickness}"
    );
    assert!(d_curvature.abs() > 1e-6, "gradient is trivially zero");
}

#[test]
fn aspheric_newton_does_not_leak_iteration_count_into_derivatives() {
    // The intersection is refined iteratively, but its derivative comes from the
    // implicit function theorem. Perturbing the seed must not move the gradient.
    let make = |c: Dual<1>| Profile::EvenAsphere {
        curvature: c,
        conic: Dual::constant(-0.6),
        coeffs: vec![Dual::constant(2.0e-6), Dual::constant(-5.0e-10)],
    };
    let origin = Vec3::new(
        Dual::constant(2.0),
        Dual::constant(-3.0),
        Dual::constant(-15.0),
    );
    let dir = Vec3::new(
        Dual::constant(0.01),
        Dual::constant(0.03),
        Dual::constant(1.0),
    )
    .normalized();

    let c0 = 1.0 / 25.0;
    let t = make(Dual::<1>::variable(c0, 0))
        .intersect(origin, dir)
        .unwrap();
    let ad = t.grad()[0];

    let value_at = |c: f64| -> f64 {
        let profile: Profile<f64> = Profile::EvenAsphere {
            curvature: c,
            conic: -0.6,
            coeffs: vec![2.0e-6, -5.0e-10],
        };
        let o = Vec3::new(2.0, -3.0, -15.0);
        let d = Vec3::new(0.01, 0.03, 1.0).normalized();
        profile.intersect(o, d).unwrap()
    };
    let h = 1e-8;
    let fd = (value_at(c0 + h) - value_at(c0 - h)) / (2.0 * h);

    assert!(
        (ad - fd).abs() / fd.abs() < 1e-6,
        "dt/dc: AD {ad} vs FD {fd}"
    );
}

#[test]
fn dual_arithmetic_obeys_the_chain_rule() {
    // f(x) = sqrt(x) * sin(x) / (1 + x^2), differentiated by hand.
    let x0 = 1.7;
    let x = Dual::<1>::variable(x0, 0);
    let f = x.sqrt() * x.sin() / (Dual::constant(1.0) + x.powi(2));
    let expected = {
        let (s, c, r) = (x0.sin(), x0.cos(), x0.sqrt());
        let d = 1.0 + x0 * x0;
        ((0.5 / r * s + r * c) * d - r * s * (2.0 * x0)) / (d * d)
    };
    assert!((f.re - x0.sqrt() * x0.sin() / (1.0 + x0 * x0)).abs() < TOL);
    assert!(
        (f.grad()[0] - expected).abs() < 1e-12,
        "{} vs {expected}",
        f.grad()[0]
    );
}

// -------------------------------------------------------------------- baselines

/// RMS spot radius in micrometres over a square grid clipped to the pupil.
fn rms_spot(sys: &System<f64>, par: &Paraxial<f64>, field: Field, wl: f64, n: usize) -> f64 {
    let mut pts = Vec::new();
    for i in 0..n {
        for j in 0..n {
            let px = -1.0 + 2.0 * (i as f64 + 0.5) / n as f64;
            let py = -1.0 + 2.0 * (j as f64 + 0.5) / n as f64;
            if px * px + py * py > 1.0 {
                continue;
            }
            if let Some(p) = trace(sys, wl, launch(sys, par, field, px, py)).image_point() {
                pts.push((p.x, p.y));
            }
        }
    }
    let k = pts.len() as f64;
    let cx = pts.iter().map(|p| p.0).sum::<f64>() / k;
    let cy = pts.iter().map(|p| p.1).sum::<f64>() / k;
    let var = pts
        .iter()
        .map(|p| (p.0 - cx).powi(2) + (p.1 - cy).powi(2))
        .sum::<f64>()
        / k;
    var.sqrt() * 1000.0
}

#[test]
fn cooke_triplet_baseline() {
    let sys = samples::cooke_triplet::<f64>();
    let par = Paraxial::compute(&sys, lines::D);

    // First-order: regression baselines, recorded 2026-09-14.
    for (label, actual, expected) in [
        ("EFL", par.efl, 50.021332),
        ("BFD", par.bfd, 42.436702),
        ("EP z", par.ep_z, 14.827966),
        ("f/#", par.fno, 5.002133),
    ] {
        assert!(
            (actual - expected).abs() < 1e-5,
            "{label}: {actual} vs baseline {expected}"
        );
    }

    // A Cooke triplet at f/5 should sit in the tens of micrometres across the field.
    // The bound is loose on purpose: it is a sanity range, not a fitted value.
    for (deg, limit) in [(0.0, 20.0), (14.0, 30.0), (20.0, 35.0)] {
        let rms = rms_spot(&sys, &par, Field::angle(deg), lines::D, 24);
        assert!(
            rms < limit,
            "field {deg}: RMS spot {rms} um exceeds {limit}"
        );
        assert!(
            rms > 1.0,
            "field {deg}: RMS spot {rms} um implausibly small"
        );
    }
}

#[test]
fn vacuum_and_mirror_media_are_consistent() {
    assert_eq!(Material::Vacuum.index(lines::D), 1.0);
    assert!(Material::Mirror.is_reflective());
    assert!(!catalog::N_BK7.is_reflective());
    assert_eq!(catalog::by_name("n-bk7"), Some(catalog::N_BK7));
    assert_eq!(catalog::by_name("NBK7"), Some(catalog::N_BK7));
    assert_eq!(catalog::by_name("nonesuch"), None);
}
