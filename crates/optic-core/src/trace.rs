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
pub fn launch<S: Scalar>(
    sys: &System<S>,
    par: &Paraxial<S>,
    field: Field,
    px: f64,
    py: f64,
) -> Ray<S> {
    let r = par.epd / S::two();
    let pupil = Vec3::new(S::from_f64(px) * r, S::from_f64(py) * r, par.ep_z);

    match (sys.object, field) {
        (Object::Infinity, Field::Angle { x, y }) => {
            let dir = Vec3::new(
                S::from_f64(x.to_radians().tan()),
                S::from_f64(y.to_radians().tan()),
                S::one(),
            );
            Ray::new(pupil, dir)
        }
        (Object::Finite { distance }, Field::Height { x, y }) => {
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
