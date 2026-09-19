//! Real (finite) sequential ray tracing.
//!
//! Every surface is visited through its [`Transform`], so the tracer is already written
//! for tilted and decentred systems even though nothing yet produces a non-identity
//! rotation.
//!
//! Optical path length accumulates from the ray's launch plane, which is the entrance
//! pupil. Referencing OPL there is what makes the wavefront analyses in `optic-analysis`
//! a subtraction rather than a separate trace.

use crate::math::{Scalar, Transform, Vec3};
use crate::paraxial::Paraxial;
use crate::surface::MissReason;
use crate::system::{Field, Object, System};

/// A geometric ray in the global frame.
#[derive(Clone, Copy, Debug)]
pub struct Ray<S: Scalar> {
    pub pos: Vec3<S>,
    /// Unit direction.
    pub dir: Vec3<S>,
    /// Accumulated optical path length.
    pub opl: S,
}

impl<S: Scalar> Ray<S> {
    pub fn new(pos: Vec3<S>, dir: Vec3<S>) -> Self {
        Self {
            pos,
            dir: dir.normalized(),
            opl: S::zero(),
        }
    }
}

/// Where a ray met one surface.
#[derive(Clone, Copy, Debug)]
pub struct Hit<S: Scalar> {
    pub surface: usize,
    /// Intersection point in the surface's local frame. This is what aperture checks,
    /// footprint diagrams and sag tests want.
    pub local: Vec3<S>,
    /// Intersection point in the global frame.
    pub global: Vec3<S>,
    /// Unit direction *after* the surface.
    pub dir: Vec3<S>,
    /// Optical path length from the launch plane to this surface.
    pub opl: S,
    /// Angle of incidence, radians. Needed for coatings and for Fresnel losses.
    pub incidence: S,
}

/// The result of pushing one ray through the system.
#[derive(Clone, Debug)]
pub struct TracedRay<S: Scalar> {
    pub hits: Vec<Hit<S>>,
    /// `Some` if the ray failed before reaching the image plane.
    pub blocked: Option<(usize, MissReason)>,
}

impl<S: Scalar> TracedRay<S> {
    /// Whether the ray reached the image plane.
    pub fn is_complete(&self) -> bool {
        self.blocked.is_none()
    }

    /// Intersection with the final surface, if the ray got there.
    pub fn image_point(&self) -> Option<Vec3<S>> {
        if self.blocked.is_some() {
            None
        } else {
            self.hits.last().map(|h| h.local)
        }
    }
}

/// Push one ray through the whole system at one wavelength.
pub fn trace<S: Scalar>(sys: &System<S>, wl: f64, mut ray: Ray<S>) -> TracedRay<S> {
    let transforms = sys.transforms();
    let mut hits = Vec::with_capacity(sys.surfaces.len());

    for (i, xform) in transforms.iter().enumerate() {
        match trace_one(sys, wl, i, xform, &mut ray) {
            Ok(hit) => hits.push(hit),
            Err(reason) => {
                return TracedRay {
                    hits,
                    blocked: Some((i, reason)),
                }
            }
        }
    }
    TracedRay {
        hits,
        blocked: None,
    }
}

fn trace_one<S: Scalar>(
    sys: &System<S>,
    wl: f64,
    i: usize,
    xform: &Transform<S>,
    ray: &mut Ray<S>,
) -> Result<Hit<S>, MissReason> {
    let surf = &sys.surfaces[i];
    let o = xform.point_to_local(ray.pos);
    let d = xform.dir_to_local(ray.dir);

    let t = surf.profile.intersect(o, d)?;
    let p = o + d * t;

    if let Some(sd) = surf.semi_diameter {
        if p.radius().value() > sd.value() {
            return Err(MissReason::ClippedByAperture);
        }
    }

    let n_before = S::from_f64(sys.medium_before(i).index(wl));
    let normal = surf.profile.normal(p)?;

    // Orient the normal against the incoming ray so the Snell branch is unambiguous.
    let mut nrm = normal;
    let mut cos_i = -d.dot(nrm);
    if cos_i.value() < 0.0 {
        nrm = -nrm;
        cos_i = -cos_i;
    }

    let new_dir = if surf.material.is_reflective() {
        d + nrm * (S::two() * cos_i)
    } else {
        let n_after = S::from_f64(sys.medium_after(i).index(wl));
        let mu = n_before / n_after;
        let sin2_t = mu * mu * (S::one() - cos_i * cos_i);
        if sin2_t.value() > 1.0 {
            return Err(MissReason::TotalInternalReflection);
        }
        let cos_t = (S::one() - sin2_t).sqrt();
        d * mu + nrm * (mu * cos_i - cos_t)
    };

    ray.opl += n_before * t;
    ray.pos = xform.point_to_global(p);
    ray.dir = xform.dir_to_global(new_dir).normalized();

    Ok(Hit {
        surface: i,
        local: p,
        global: ray.pos,
        dir: ray.dir,
        opl: ray.opl,
        incidence: cos_i.acos(),
    })
}

/// Build the ray that leaves `field` and crosses the entrance pupil at normalised
/// coordinates `(px, py)`, each in `[-1, 1]`.
///
/// This is paraxial pupil aiming: the ray is aimed at the *paraxial* entrance pupil, so
/// in a system with strong pupil aberration the real pupil is not filled uniformly. Real
/// ray aiming (iterating on the stop intersection) is a later refinement; the interface
/// does not change when it arrives.
///
/// The field's vignetting factors are applied here, and only here. `(px, py)` therefore
/// always means "the full pupil" to the caller, and whatever part of it this field
/// actually uses is settled in one place rather than at every call site.
pub fn launch<S: Scalar>(
    sys: &System<S>,
    par: &Paraxial<S>,
    field: Field,
    px: f64,
    py: f64,
) -> Ray<S> {
    let (px, py) = field.vignette().apply(px, py);
    let r = par.epd / S::two();
    let pupil = Vec3::new(S::from_f64(px) * r, S::from_f64(py) * r, par.ep_z);

    match (sys.object, field) {
        (Object::Infinity, Field::Angle { x, y, .. }) => {
            let dir = Vec3::new(
                S::from_f64(x.to_radians().tan()),
                S::from_f64(y.to_radians().tan()),
                S::one(),
            );
            Ray::new(pupil, dir)
        }
        (Object::Finite { distance }, Field::Height { x, y, .. }) => {
            let src = Vec3::new(S::from_f64(x), S::from_f64(y), -distance);
            Ray::new(src, pupil - src)
        }
        // A mismatched field type is a document error; degrade to the axial ray rather
        // than producing silently meaningless output.
        (Object::Infinity, Field::Height { .. }) => Ray::new(pupil, Vec3::axis()),
        (Object::Finite { distance }, Field::Angle { .. }) => {
            let src = Vec3::new(S::zero(), S::zero(), -distance);
            Ray::new(src, pupil - src)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::material::{lines, Material};
    use crate::system::{Aperture, Surface};

    const N_PLATE: f64 = 1.5;
    const T_PLATE: f64 = 8.0;

    /// A plane-parallel plate in vacuum, with the image plane 20 mm past it.
    fn plate(n: f64) -> System<f64> {
        System::new(
            Object::Infinity,
            vec![
                Surface::plane(T_PLATE, Material::Fixed { n }).stop(),
                Surface::plane(20.0, Material::Vacuum),
                Surface::plane(0.0, Material::Vacuum),
            ],
        )
        .with_aperture(Aperture::EntrancePupilDiameter(4.0))
    }

    fn ray_at(theta: f64) -> Ray<f64> {
        Ray::new(
            Vec3::new(0.0, 0.0, -5.0),
            Vec3::new(0.0, theta.sin(), theta.cos()),
        )
    }

    #[test]
    fn vignetting_factors_shrink_the_beam_that_is_launched() {
        use crate::system::Vignette;
        let sys = crate::samples::cooke_triplet::<f64>();
        let par = Paraxial::compute(&sys, lines::D);
        let field = Field::angle(10.0);

        let full = launch(&sys, &par, field, 0.0, 1.0);
        let squeezed = launch(
            &sys,
            &par,
            field.vignetted(Vignette {
                cy: 0.5,
                ..Vignette::NONE
            }),
            0.0,
            1.0,
        );

        // Same field, so the ray still travels in the same direction; it just enters
        // through half the pupil height.
        assert!((squeezed.pos.y - full.pos.y * 0.5).abs() < 1e-12);
        assert!((squeezed.dir - full.dir).norm() < 1e-15);

        // A decentred pupil moves the ray the caller thinks of as the chief ray.
        let shifted = launch(
            &sys,
            &par,
            field.vignetted(Vignette {
                dy: 0.5,
                ..Vignette::NONE
            }),
            0.0,
            0.0,
        );
        assert!((shifted.pos.y - full.pos.y * 0.5).abs() < 1e-12);
    }

    #[test]
    fn a_new_ray_is_normalised() {
        let r = Ray::new(Vec3::zero(), Vec3::new(0.0, 3.0, 4.0));
        assert!((r.dir.norm() - 1.0).abs() < 1e-15);
        assert!((r.dir - Vec3::new(0.0, 0.6, 0.8)).norm() < 1e-15);
        assert_eq!(r.opl, 0.0);
    }

    #[test]
    fn a_plate_displaces_a_ray_without_deviating_it() {
        // Closed form: the emergent ray is parallel to the incident one, and the
        // transverse offset is t (tan A - tan A'), with A' from Snell's law.
        let sys = plate(N_PLATE);
        for theta in [0.1f64, 0.3, 0.5] {
            let start = ray_at(theta);
            let traced = trace(&sys, lines::D, start);
            assert!(traced.is_complete());

            let exit = traced.hits.last().unwrap();
            assert!(
                (exit.dir - start.dir).norm() < 1e-14,
                "the plate deviated the ray at {theta} rad"
            );

            let theta_p = (theta.sin() / N_PLATE).asin();
            let undeviated = 33.0 * theta.tan();
            let expected = undeviated - T_PLATE * (theta.tan() - theta_p.tan());
            let got = traced.image_point().unwrap().y;
            assert!((got - expected).abs() < 1e-12, "{got} vs {expected}");
        }
    }

    #[test]
    fn a_plate_of_unit_index_does_nothing_at_all() {
        let traced = trace(&plate(1.0), lines::D, ray_at(0.4));
        let expected = 33.0 * 0.4f64.tan();
        assert!((traced.image_point().unwrap().y - expected).abs() < 1e-13);
    }

    #[test]
    fn optical_path_length_counts_each_medium_at_its_own_index() {
        let traced = trace(&plate(N_PLATE), lines::D, ray_at(0.0));
        // 5 mm of vacuum, 8 mm of glass, then 20 mm of vacuum.
        let expected = 5.0 + N_PLATE * T_PLATE + 20.0;
        let got = traced.hits.last().unwrap().opl;
        assert!((got - expected).abs() < 1e-13, "OPL {got} vs {expected}");

        // It must accumulate monotonically, surface by surface.
        let mut previous = 0.0;
        for h in &traced.hits {
            assert!(h.opl > previous);
            previous = h.opl;
        }
    }

    #[test]
    fn angles_of_incidence_are_recorded() {
        let theta = 0.4;
        let traced = trace(&plate(N_PLATE), lines::D, ray_at(theta));
        assert!((traced.hits[0].incidence - theta).abs() < 1e-14);
        // Leaving the plate, the internal angle is the refracted one.
        let theta_p = (theta.sin() / N_PLATE).asin();
        assert!((traced.hits[1].incidence - theta_p).abs() < 1e-14);
        assert!((traced.hits[0].incidence - 0.0).abs() > 1e-9);
    }

    #[test]
    fn total_internal_reflection_stops_the_ray() {
        // Starting inside a dense medium, past the critical angle of 33.7 degrees.
        let mut sys: System<f64> = System::new(
            Object::Infinity,
            vec![
                Surface::plane(10.0, Material::Vacuum).stop(),
                Surface::plane(0.0, Material::Vacuum),
            ],
        );
        sys.object_medium = Material::Fixed { n: 1.8 };

        let steep = trace(&sys, lines::D, ray_at(0.785)); // 45 degrees
        assert_eq!(
            steep.blocked,
            Some((0, MissReason::TotalInternalReflection))
        );
        assert!(!steep.is_complete());
        assert!(steep.image_point().is_none());
        assert!(
            steep.hits.is_empty(),
            "no hit is recorded for a blocked surface"
        );

        let shallow = trace(&sys, lines::D, ray_at(0.3)); // 17 degrees, well under
        assert!(shallow.is_complete());
    }

    #[test]
    fn a_clear_aperture_blocks_what_falls_outside_it() {
        let mut sys = plate(N_PLATE);
        sys.surfaces[0].semi_diameter = Some(1.0);

        let inside = trace(
            &sys,
            lines::D,
            Ray::new(Vec3::new(0.0, 0.5, -5.0), Vec3::axis()),
        );
        assert!(inside.is_complete());

        let outside = trace(
            &sys,
            lines::D,
            Ray::new(Vec3::new(0.0, 5.0, -5.0), Vec3::axis()),
        );
        assert_eq!(outside.blocked, Some((0, MissReason::ClippedByAperture)));

        // Hits before the blocking surface are still reported, for footprint plots.
        sys.surfaces[0].semi_diameter = None;
        sys.surfaces[1].semi_diameter = Some(1.0);
        let late = trace(
            &sys,
            lines::D,
            Ray::new(Vec3::new(0.0, 5.0, -5.0), Vec3::axis()),
        );
        assert_eq!(late.blocked, Some((1, MissReason::ClippedByAperture)));
        assert_eq!(late.hits.len(), 1);
    }

    #[test]
    fn a_flat_mirror_reverses_an_axial_ray() {
        let sys: System<f64> = System::new(
            Object::Infinity,
            vec![
                Surface::plane(-10.0, Material::Mirror).stop(),
                Surface::plane(0.0, Material::Vacuum),
            ],
        );
        let traced = trace(
            &sys,
            lines::D,
            Ray::new(Vec3::new(0.0, 2.0, -5.0), Vec3::axis()),
        );
        assert!(traced.is_complete());
        assert_eq!(traced.hits[0].dir.value(), [0.0, 0.0, -1.0]);
        // Reflected straight back, the ray returns to the height it arrived at.
        assert!((traced.image_point().unwrap().y - 2.0).abs() < 1e-14);
    }

    #[test]
    fn launching_from_infinity_sets_the_direction_from_the_field_angle() {
        let sys = crate::samples::cooke_triplet::<f64>();
        let par = Paraxial::compute(&sys, lines::D);

        let r = launch(&sys, &par, Field::angle_xy(0.0, 10.0), 0.0, 0.0);
        let expected = Vec3::new(0.0, 10f64.to_radians().tan(), 1.0).normalized();
        assert!((r.dir - expected).norm() < 1e-15);
        // The pupil centre, at the entrance pupil plane.
        assert_eq!(r.pos.value(), [0.0, 0.0, par.ep_z]);

        // Normalised pupil coordinates span the pupil diameter.
        let edge = launch(&sys, &par, Field::angle(0.0), 0.0, 1.0);
        assert!((edge.pos.y - par.epd / 2.0).abs() < 1e-15);
        let corner = launch(&sys, &par, Field::angle(0.0), -1.0, 0.0);
        assert!((corner.pos.x + par.epd / 2.0).abs() < 1e-15);
    }

    #[test]
    fn launching_from_a_finite_object_aims_from_the_object_point() {
        let sys: System<f64> = System::new(
            Object::Finite { distance: 200.0 },
            vec![
                Surface::new(60.0, 5.0, Material::Fixed { n: 1.5 }).stop(),
                Surface::new(-60.0, 100.0, Material::Vacuum),
                Surface::plane(0.0, Material::Vacuum),
            ],
        )
        .with_aperture(Aperture::EntrancePupilDiameter(10.0))
        .with_fields(vec![Field::height(-8.0)]);

        let par = Paraxial::compute(&sys, lines::D);
        let r = launch(&sys, &par, Field::height_xy(0.0, -8.0), 0.0, 1.0);
        assert_eq!(r.pos.value(), [0.0, -8.0, -200.0]);

        // It must point at the top of the entrance pupil.
        let target = Vec3::new(0.0, par.epd / 2.0, par.ep_z);
        assert!(r.dir.cross(target - r.pos).norm() < 1e-14);
        assert!(trace(&sys, lines::D, r).is_complete());
    }

    #[test]
    fn rays_are_traced_through_each_surfaces_own_frame() {
        // Every hit records its local point; on a centred system that differs from the
        // global point by exactly the vertex position.
        let sys = crate::samples::cooke_triplet::<f64>();
        let par = Paraxial::compute(&sys, lines::D);
        let traced = trace(
            &sys,
            lines::D,
            launch(&sys, &par, Field::angle(10.0), 0.3, 0.4),
        );
        for (h, z) in traced.hits.iter().zip(sys.vertices()) {
            assert!((h.global.z - h.local.z - z).abs() < 1e-12);
            assert_eq!(h.global.x, h.local.x);
            assert_eq!(h.global.y, h.local.y);
        }
    }

    #[test]
    fn every_hit_direction_is_a_unit_vector() {
        let sys = crate::samples::cooke_triplet::<f64>();
        let par = Paraxial::compute(&sys, lines::D);
        for (px, py) in [(0.0, 0.0), (0.9, 0.0), (-0.5, 0.6)] {
            let traced = trace(
                &sys,
                lines::D,
                launch(&sys, &par, Field::angle(14.0), px, py),
            );
            for h in &traced.hits {
                assert!((h.dir.norm() - 1.0).abs() < 1e-14, "surface {}", h.surface);
            }
        }
    }

    #[test]
    fn a_traced_ray_reports_one_hit_per_surface() {
        let sys = crate::samples::cooke_triplet::<f64>();
        let par = Paraxial::compute(&sys, lines::D);
        let traced = trace(
            &sys,
            lines::D,
            launch(&sys, &par, Field::angle(0.0), 0.0, 0.5),
        );
        assert_eq!(traced.hits.len(), sys.surfaces.len());
        for (i, h) in traced.hits.iter().enumerate() {
            assert_eq!(h.surface, i);
        }
    }
}
