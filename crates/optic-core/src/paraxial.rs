//! First-order optics: the y-u trace, focal length, and the pupils.
//!
//! Real rays cannot be launched until we know where the entrance pupil is, so this runs
//! first for every trace. Angles are true angles, not reduced; the reduced form appears
//! only inside the refraction step.
//!
//! Sign convention: distances and angles are positive in the `+z` sense. A mirror
//! negates the index of the following medium, which is why thicknesses after a mirror
//! are written negative, exactly as in a conventional prescription.

use crate::math::Scalar;
use crate::system::{Aperture, Object, System};

/// Height and angle of a paraxial ray at one surface.
#[derive(Clone, Copy, Debug)]
pub struct ParaxialState<S: Scalar> {
    /// Ray height at the surface.
    pub y: S,
    /// Ray angle in the medium *following* the surface.
    pub u: S,
    /// Index of the medium following the surface, signed by reflection parity.
    pub n: S,
}

/// First-order properties of a system at one wavelength.
#[derive(Clone, Debug)]
pub struct Paraxial<S: Scalar> {
    /// Effective focal length, referred to image space.
    pub efl: S,
    /// Back focal distance: last surface vertex to the paraxial focus.
    pub bfd: S,
    /// Entrance pupil diameter.
    pub epd: S,
    /// Entrance pupil position, relative to the first surface vertex.
    pub ep_z: S,
    /// Exit pupil position, relative to the image plane.
    pub xp_z: S,
    /// Image-space working f-number.
    pub fno: S,
    /// Marginal (axial edge-of-pupil) ray at each surface.
    pub marginal: Vec<ParaxialState<S>>,
    /// Chief (field edge, pupil centre) ray at each surface.
    pub chief: Vec<ParaxialState<S>>,
}

/// Signed indices of the medium before and after each surface, tracking reflections.
fn signed_indices<S: Scalar>(sys: &System<S>, wl: f64) -> Vec<(S, S)> {
    let mut parity = 1.0;
    let mut out = Vec::with_capacity(sys.surfaces.len());
    for i in 0..sys.surfaces.len() {
        let n_before = sys.medium_before(i).index(wl) * parity;
        if sys.surfaces[i].material.is_reflective() {
            parity = -parity;
        }
        let after = sys.medium_after(i);
        let n_after = if after.is_reflective() {
            // A mirror leaves the medium unchanged but reverses the direction.
            sys.medium_before(i).index(wl) * parity
        } else {
            after.index(wl) * parity
        };
        out.push((S::from_f64(n_before), S::from_f64(n_after)));
    }
    out
}

/// March a paraxial ray from object space through to the image plane.
///
/// `y0`/`u0` are the height at the first surface vertex plane and the angle in object
/// space. Returns the state at every surface.
pub fn forward<S: Scalar>(sys: &System<S>, wl: f64, y0: S, u0: S) -> Vec<ParaxialState<S>> {
    let idx = signed_indices(sys, wl);
    let mut states = Vec::with_capacity(sys.surfaces.len());
    let (mut y, mut u) = (y0, u0);

    for (i, surf) in sys.surfaces.iter().enumerate() {
        let (n, np) = idx[i];
        let c = surf.profile.paraxial_curvature();
        // n' u' = n u - y c (n' - n)
        u = (n * u - y * c * (np - n)) / np;
        states.push(ParaxialState { y, u, n: np });
        y += u * surf.thickness;
    }
    states
}

/// March a paraxial ray backwards from just before surface `from` into object space.
///
/// Returns the height at the first surface vertex plane and the angle in object space.
fn backward<S: Scalar>(sys: &System<S>, wl: f64, from: usize, y: S, u: S) -> (S, S) {
    let idx = signed_indices(sys, wl);
    let (mut y, mut u) = (y, u);
    for i in (0..from).rev() {
        let (n, np) = idx[i];
        // Undo the transfer from surface i to surface i+1, then undo the refraction.
        y -= u * sys.surfaces[i].thickness;
        u = (np * u + y * sys.surfaces[i].profile.paraxial_curvature() * (np - n)) / n;
    }
    (y, u)
}

impl<S: Scalar> Paraxial<S> {
    /// Compute the first-order properties of `sys` at `wl` micrometres.
    pub fn compute(sys: &System<S>, wl: f64) -> Self {
        // Effective focal length is defined by a collimated unit-height input,
        // independent of where the object actually is.
        let probe = forward(sys, wl, S::one(), S::zero());
        let k = sys.image_index();
        let u_img = probe[k].u;
        let efl = -S::one() / u_img;
        // Back focal distance is measured from the last *refracting* surface, not from
        // wherever the image plane currently sits -- otherwise refocusing oscillates
        // about the true focus instead of converging on it.
        let last_powered = if k == 0 { 0 } else { k - 1 };
        let bfd = -probe[last_powered].y / probe[last_powered].u;

        let (ep_z, ep_scale) = entrance_pupil(sys, wl);
        let epd = resolve_epd(sys, wl, efl, ep_scale);

        // Marginal ray: from the axial object point through the edge of the entrance pupil.
        let (my0, mu0) = match sys.object {
            Object::Infinity => (epd / S::two(), S::zero()),
            Object::Finite { distance } => {
                let u = (epd / S::two()) / (distance + ep_z);
                (u * distance, u)
            }
        };
        let marginal = forward(sys, wl, my0, mu0);

        // Chief ray for the widest field in the set: this is the one that sizes
        // apertures and locates the exit pupil.
        let widest = sys
            .fields
            .iter()
            .copied()
            .max_by(|a, b| a.radius().total_cmp(&b.radius()))
            .unwrap_or(crate::system::Field::angle(0.0));
        let chief = chief_ray(sys, wl, widest, ep_z);

        // Exit pupil: where the chief ray crosses the axis in image space.
        let c_img = chief[sys.image_index()];
        let xp_z = if c_img.u.value().abs() > 0.0 {
            -c_img.y / c_img.u
        } else {
            S::from_f64(f64::NEG_INFINITY)
        };

        let m_img = marginal[sys.image_index()];
        let fno = if m_img.u.value().abs() > 0.0 {
            (S::one() / (S::two() * m_img.u)).abs()
        } else {
            S::from_f64(f64::INFINITY)
        };

        Self {
            efl,
            bfd,
            epd,
            ep_z,
            xp_z,
            fno,
            marginal,
            chief,
        }
    }
}

/// Locate the entrance pupil: the image of the stop in object space.
///
/// Returns its axial position relative to the first surface vertex, and the ratio of
/// entrance pupil radius to stop radius.
fn entrance_pupil<S: Scalar>(sys: &System<S>, wl: f64) -> (S, S) {
    let Some(stop) = sys.stop_index() else {
        // No stop marked: the first surface is the pupil.
        return (S::zero(), S::one());
    };
    if stop == 0 {
        return (S::zero(), S::one());
    }
    // A ray leaving the stop centre, and one leaving its edge parallel to the axis.
    let (y_c, u_c) = backward(sys, wl, stop, S::zero(), S::one());
    let (y_e, u_e) = backward(sys, wl, stop, S::one(), S::zero());

    let z = if u_c.value().abs() > 0.0 {
        -y_c / u_c
    } else {
        S::zero()
    };
    let scale = y_e + u_e * z;
    (z, scale.abs())
}

/// Turn whatever aperture the user specified into an entrance pupil diameter.
fn resolve_epd<S: Scalar>(sys: &System<S>, wl: f64, efl: S, ep_scale: S) -> S {
    match sys.aperture {
        Aperture::EntrancePupilDiameter(d) => S::from_f64(d),
        Aperture::ImageSpaceFNumber(f) => (efl / S::from_f64(f)).abs(),
        Aperture::StopDiameter(d) => S::from_f64(d) * ep_scale,
        Aperture::ObjectSpaceNA(na) => match sys.object {
            Object::Finite { distance } => {
                let n = S::from_f64(sys.object_medium.index(wl));
                let u = S::from_f64(na) / n;
                // Small-angle: half-diameter = u * (object to entrance pupil distance).
                let (ep_z, _) = entrance_pupil(sys, wl);
                (u * (distance + ep_z)).abs() * S::two()
            }
            Object::Infinity => S::from_f64(f64::NAN),
        },
    }
}

/// Paraxial chief ray for one field: through the centre of the entrance pupil.
pub fn chief_ray<S: Scalar>(
    sys: &System<S>,
    wl: f64,
    field: crate::system::Field,
    ep_z: S,
) -> Vec<ParaxialState<S>> {
    let (y0, u0) = match (sys.object, field) {
        (Object::Infinity, crate::system::Field::Angle { x, y }) => {
            let r = (x * x + y * y).sqrt();
            if r == 0.0 {
                (S::zero(), S::zero())
            } else {
                let u = S::from_f64(r.to_radians().tan());
                (-u * ep_z, u)
            }
        }
        (Object::Finite { distance }, crate::system::Field::Height { x, y }) => {
            let h = S::from_f64((x * x + y * y).sqrt());
            let u = -h / (distance + ep_z);
            (h + u * distance, u)
        }
        _ => (S::zero(), S::zero()),
    };
    forward(sys, wl, y0, u0)
}

impl<S: Scalar> Paraxial<S> {
    /// Paraxial image height for one field point.
    pub fn image_height(&self, sys: &System<S>, wl: f64, field: crate::system::Field) -> S {
        chief_ray(sys, wl, field, self.ep_z)[sys.image_index()].y
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::material::{catalog, lines, Material};
    use crate::surface::Profile;
    use crate::system::{Aperture, Field, Surface, Wavelength};

    /// A thin lens in vacuum: two surfaces with no glass thickness between them.
    fn thin_lens(r1: f64, r2: f64, n: f64) -> System<f64> {
        System::new(
            Object::Infinity,
            vec![
                Surface::new(r1, 0.0, Material::Fixed { n }).stop(),
                Surface::new(r2, 100.0, Material::Vacuum),
                Surface::plane(0.0, Material::Vacuum),
            ],
        )
        .with_aperture(Aperture::EntrancePupilDiameter(4.0))
    }

    #[test]
    fn a_thin_lens_obeys_the_lensmakers_equation() {
        for (r1, r2, n) in [
            (100.0, -100.0, 1.5),
            (50.0, f64::INFINITY, 1.5),
            (f64::INFINITY, -75.0, 1.62),
            (-60.0, 60.0, 1.5), // a negative lens
        ] {
            let expected = 1.0 / ((n - 1.0) * (1.0 / r1 - 1.0 / r2));
            let efl = Paraxial::compute(&thin_lens(r1, r2, n), lines::D).efl;
            assert!(
                (efl - expected).abs() < 1e-9,
                "R1={r1} R2={r2}: EFL {efl} vs {expected}"
            );
        }
    }

    #[test]
    fn a_negative_lens_has_a_negative_focal_length() {
        assert!(Paraxial::compute(&thin_lens(-60.0, 60.0, 1.5), lines::D).efl < 0.0);
    }

    #[test]
    fn a_concave_mirror_focuses_at_half_its_radius() {
        // The reflective path flips the sign of the following index; nothing else in
        // the system exercises it, so it is checked in isolation against f = R/2.
        for radius in [-100.0, -250.0] {
            let sys: System<f64> = System::new(
                Object::Infinity,
                vec![
                    Surface::new(radius, radius / 2.0, Material::Mirror).stop(),
                    Surface::plane(0.0, Material::Vacuum),
                ],
            )
            .with_aperture(Aperture::EntrancePupilDiameter(20.0));

            let par = Paraxial::compute(&sys, lines::D);
            assert!(
                (par.efl - radius / 2.0).abs() < 1e-10,
                "R={radius}: EFL {} vs {}",
                par.efl,
                radius / 2.0
            );
            assert!(
                (par.bfd - radius / 2.0).abs() < 1e-10,
                "R={radius}: BFD {}",
                par.bfd
            );
        }
    }

    #[test]
    fn reflection_flips_the_sign_of_every_index_downstream() {
        let sys: System<f64> = System::new(
            Object::Infinity,
            vec![
                Surface::new(-100.0, -50.0, Material::Mirror).stop(),
                Surface::plane(0.0, Material::Vacuum),
            ],
        );
        let idx = signed_indices(&sys, lines::D);
        assert_eq!(idx[0].0, 1.0); // before the mirror, forward
        assert_eq!(idx[0].1, -1.0); // after it, folded
        assert_eq!(idx[1].0, -1.0); // and it stays folded
        assert_eq!(idx[1].1, -1.0);
    }

    #[test]
    fn dispersion_shortens_the_focal_length_in_the_blue() {
        let sys: System<f64> = System::new(
            Object::Infinity,
            vec![
                Surface::new(80.0, 5.0, catalog::N_BK7).stop(),
                Surface::new(-80.0, 80.0, Material::Vacuum),
                Surface::plane(0.0, Material::Vacuum),
            ],
        )
        .with_wavelengths(vec![Wavelength::new(lines::F), Wavelength::new(lines::C)]);

        let blue = Paraxial::compute(&sys, lines::F).efl;
        let red = Paraxial::compute(&sys, lines::C).efl;
        assert!(
            blue < red,
            "axial colour has the wrong sign: {blue} vs {red}"
        );
        assert!((red - blue) / red > 1e-3, "axial colour implausibly small");
    }

    #[test]
    fn tracing_backwards_undoes_tracing_forwards() {
        let sys = crate::samples::cooke_triplet::<f64>();
        let (y0, u0) = (0.7, -0.013);
        let states = forward(&sys, lines::D, y0, u0);

        for from in 1..sys.surfaces.len() {
            let (y, u) = backward(&sys, lines::D, from, states[from].y, states[from - 1].u);
            assert!((y - y0).abs() < 1e-12, "from {from}: height {y} vs {y0}");
            assert!((u - u0).abs() < 1e-12, "from {from}: angle {u} vs {u0}");
        }
    }

    #[test]
    fn a_stop_at_the_front_is_its_own_entrance_pupil() {
        let sys = thin_lens(100.0, -100.0, 1.5);
        let (z, scale) = entrance_pupil(&sys, lines::D);
        assert_eq!(z, 0.0);
        assert_eq!(scale, 1.0);
    }

    #[test]
    fn aiming_at_the_entrance_pupil_lands_on_the_stop() {
        // This is the definition of the entrance pupil, and the only check of it worth
        // making: a chief ray aimed at the pupil centre must cross the axis exactly at
        // the stop, and a marginal ray at the pupil edge must graze the stop's rim.
        // Where the pupil physically sits -- it is a virtual image, and in this triplet
        // it falls behind the stop -- says nothing about whether it is correct.
        let sys = crate::samples::cooke_triplet::<f64>();
        let par = Paraxial::compute(&sys, lines::D);
        let stop = sys.stop_index().unwrap();

        for deg in [5.0, 14.0, 20.0] {
            let chief = chief_ray(&sys, lines::D, Field::angle(deg), par.ep_z);
            assert!(
                chief[stop].y.abs() < 1e-12,
                "chief ray at {deg} deg misses the stop centre by {}",
                chief[stop].y
            );
        }

        let (_, scale) = entrance_pupil(&sys, lines::D);
        let expected_stop_radius = par.epd / 2.0 / scale;
        assert!(
            (par.marginal[stop].y.abs() - expected_stop_radius).abs() < 1e-12,
            "marginal ray height at the stop {} vs expected rim {expected_stop_radius}",
            par.marginal[stop].y
        );
        assert!(scale > 0.0 && scale.is_finite());
    }

    #[test]
    fn every_aperture_specification_resolves_to_the_same_pupil() {
        let base = crate::samples::cooke_triplet::<f64>();
        let reference = Paraxial::compute(&base, lines::D);

        // Restating the same aperture three different ways must not change the system.
        let mut by_fno = base.clone();
        by_fno.aperture = Aperture::ImageSpaceFNumber(reference.efl / reference.epd);
        assert!((Paraxial::compute(&by_fno, lines::D).epd - reference.epd).abs() < 1e-9);

        let (_, scale) = entrance_pupil(&base, lines::D);
        let mut by_stop = base.clone();
        by_stop.aperture = Aperture::StopDiameter(reference.epd / scale);
        assert!((Paraxial::compute(&by_stop, lines::D).epd - reference.epd).abs() < 1e-9);
    }

    #[test]
    fn the_working_f_number_follows_the_focal_ratio() {
        let sys = crate::samples::cooke_triplet::<f64>();
        let par = Paraxial::compute(&sys, lines::D);
        // At infinite conjugates the working f/# is EFL / EPD to first order.
        assert!(
            (par.fno - par.efl / par.epd).abs() < 1e-6,
            "f/# {}",
            par.fno
        );
    }

    #[test]
    fn the_chief_ray_is_flat_on_axis() {
        let sys = crate::samples::cooke_triplet::<f64>();
        let par = Paraxial::compute(&sys, lines::D);
        let axial = chief_ray(&sys, lines::D, Field::angle(0.0), par.ep_z);
        for s in &axial {
            assert_eq!(s.y, 0.0);
            assert_eq!(s.u, 0.0);
        }
        assert_eq!(par.image_height(&sys, lines::D, Field::angle(0.0)), 0.0);
    }

    #[test]
    fn image_height_grows_with_field_angle() {
        let sys = crate::samples::cooke_triplet::<f64>();
        let par = Paraxial::compute(&sys, lines::D);
        let mut previous = 0.0;
        for deg in [1.0, 5.0, 10.0, 20.0] {
            let h = par.image_height(&sys, lines::D, Field::angle(deg));
            assert!(h > previous, "height did not grow at {deg} degrees");
            previous = h;
        }
    }

    #[test]
    fn a_mismatched_field_type_degrades_to_the_axis() {
        // A finite-conjugate field angle, or an infinite-conjugate height, is a document
        // error. It must not silently produce a plausible-looking wrong answer.
        let sys = crate::samples::cooke_triplet::<f64>();
        let states = chief_ray(&sys, lines::D, Field::height(5.0), 10.0);
        assert!(states.iter().all(|s| s.y == 0.0 && s.u == 0.0));
    }

    #[test]
    fn aspheric_terms_above_r2_do_not_disturb_first_order_optics() {
        let mut sys = thin_lens(100.0, -100.0, 1.5);
        let before = Paraxial::compute(&sys, lines::D).efl;

        // alpha-1 is zero here, so only r^4 and r^6 are present: fourth order and above,
        // invisible to the paraxial trace.
        sys.surfaces[0].profile = Profile::EvenAsphere {
            curvature: 1.0 / 100.0,
            conic: -2.5,
            coeffs: vec![0.0, 1e-6, -3e-9],
        };
        assert!((Paraxial::compute(&sys, lines::D).efl - before).abs() < 1e-12);

        // An r^2 term, by contrast, must move the focal length.
        sys.surfaces[0].profile = Profile::EvenAsphere {
            curvature: 1.0 / 100.0,
            conic: -2.5,
            coeffs: vec![1e-4],
        };
        let after = Paraxial::compute(&sys, lines::D).efl;
        assert!(
            (after - before).abs() > 0.1,
            "an r^2 coefficient left the focal length at {before}"
        );
    }
}
