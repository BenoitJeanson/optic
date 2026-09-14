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
