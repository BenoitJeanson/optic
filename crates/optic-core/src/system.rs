//! The optical system: the document everything else is a pure function of.
//!
//! Two departures from the Zemax data model, both deliberate:
//!
//! * The object is not row 0 of the surface list. An object at infinity would put
//!   `INFINITY` into the scalar type, and infinity times a zero derivative is `NaN` --
//!   one poisoned entry propagates through the entire Jacobian. Modelling the object
//!   as its own sum type keeps every number in the trace finite. The editor can still
//!   *present* it as row 0.
//! * A surface owns the medium that *follows* it, which is how prescriptions are
//!   written and read.

use crate::material::Material;
use crate::math::{Scalar, Transform};
use crate::solve::ThicknessSolve;
use crate::surface::Profile;

/// Where the object lives.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(tag = "at", rename_all = "snake_case"))]
pub enum Object<S: Scalar> {
    /// Collimated input. Fields are angles.
    Infinity,
    /// A finite conjugate `distance` in front of the first surface. Fields are heights.
    Finite { distance: S },
}

/// One row of the prescription.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Surface<S: Scalar> {
    pub profile: Profile<S>,
    /// Axial distance to the next surface.
    ///
    /// When [`Surface::thickness_solve`] is not `Fixed`, this holds the *solved* value:
    /// the solve writes here, so everything downstream reads a thickness without caring
    /// whether a person or a constraint put it there.
    pub thickness: S,
    /// How `thickness` is determined.
    #[cfg_attr(feature = "serde", serde(default = "ThicknessSolve::default"))]
    pub thickness_solve: ThicknessSolve<S>,
    /// The medium *after* this surface.
    pub material: Material,
    /// Clear semi-aperture. `None` means "as large as needed", resolved by ray tracing.
    pub semi_diameter: Option<S>,
    /// Marks the aperture stop. Exactly one surface should carry it.
    pub is_stop: bool,
    pub label: String,
}

impl<S: Scalar> Surface<S> {
    /// A spherical refracting surface.
    pub fn new(radius: f64, thickness: f64, material: Material) -> Self {
        Self {
            profile: Profile::sphere(radius),
            thickness: S::from_f64(thickness),
            thickness_solve: ThicknessSolve::Fixed,
            material,
            semi_diameter: None,
            is_stop: false,
            label: String::new(),
        }
    }

    /// A flat surface, typically a stop or the image plane.
    pub fn plane(thickness: f64, material: Material) -> Self {
        Self::new(f64::INFINITY, thickness, material)
    }

    pub fn stop(mut self) -> Self {
        self.is_stop = true;
        self
    }

    pub fn labelled(mut self, label: &str) -> Self {
        self.label = label.to_string();
        self
    }

    pub fn with_semi_diameter(mut self, sd: f64) -> Self {
        self.semi_diameter = Some(S::from_f64(sd));
        self
    }

    /// Let a constraint determine this thickness instead of a typed value.
    pub fn solved_by(mut self, solve: ThicknessSolve<S>) -> Self {
        self.thickness_solve = solve;
        self
    }

    /// Autofocus: place the next surface where the paraxial marginal ray crosses the axis.
    pub fn autofocus(self) -> Self {
        self.solved_by(ThicknessSolve::MarginalRayHeight { height: S::zero() })
    }
}

/// One wavelength in the analysis set.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Wavelength {
    /// Micrometres.
    pub um: f64,
    pub weight: f64,
}

impl Wavelength {
    pub fn new(um: f64) -> Self {
        Self { um, weight: 1.0 }
    }
}

/// How much of the pupil a field point actually uses.
///
/// A real lens loses light at the edge of the field to mounts, barrels and the rims of
/// its own elements. Modelling every one of those is possible but tedious, so the trade
/// is to describe the surviving beam directly: shrink and shift the pupil per field
/// until it matches the light that gets through. Zemax calls the five numbers vignetting
/// factors and applies them to the normalised pupil coordinates before launching a ray,
/// which is what we do here.
///
/// All zero means the whole pupil, and is the default.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Vignette {
    /// Pupil decentre in x, in normalised pupil coordinates.
    #[cfg_attr(feature = "serde", serde(default))]
    pub dx: f64,
    /// Pupil decentre in y.
    #[cfg_attr(feature = "serde", serde(default))]
    pub dy: f64,
    /// Pupil compression in x: 0 keeps the full width, 0.5 keeps half of it.
    #[cfg_attr(feature = "serde", serde(default))]
    pub cx: f64,
    /// Pupil compression in y.
    #[cfg_attr(feature = "serde", serde(default))]
    pub cy: f64,
    /// Rotation of the vignetted pupil, in degrees.
    #[cfg_attr(feature = "serde", serde(default))]
    pub angle: f64,
}

impl Vignette {
    /// The whole pupil: no decentre, no compression, no rotation.
    pub const NONE: Vignette = Vignette {
        dx: 0.0,
        dy: 0.0,
        cx: 0.0,
        cy: 0.0,
        angle: 0.0,
    };

    /// Whether this field uses the pupil in full.
    pub fn is_full(&self) -> bool {
        *self == Vignette::NONE
    }

    /// Map a point on the full pupil onto the part of it this field actually uses.
    ///
    /// Compress and decentre in x and y, then rotate: the order Zemax applies, and the
    /// order that makes the rotation act on the vignetted pupil rather than on the
    /// full one.
    pub fn apply(&self, px: f64, py: f64) -> (f64, f64) {
        let x = self.dx + px * (1.0 - self.cx);
        let y = self.dy + py * (1.0 - self.cy);
        if self.angle == 0.0 {
            return (x, y);
        }
        let (s, c) = self.angle.to_radians().sin_cos();
        (x * c - y * s, x * s + y * c)
    }
}

/// One field point.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(tag = "kind", rename_all = "snake_case"))]
pub enum Field {
    /// Incidence angle in degrees. Only meaningful for an object at infinity.
    Angle {
        x: f64,
        y: f64,
        #[cfg_attr(feature = "serde", serde(default))]
        vignette: Vignette,
    },
    /// Object height in lens units. Only meaningful for a finite object.
    Height {
        x: f64,
        y: f64,
        #[cfg_attr(feature = "serde", serde(default))]
        vignette: Vignette,
    },
}

impl Field {
    pub fn angle(y: f64) -> Self {
        Field::Angle {
            x: 0.0,
            y,
            vignette: Vignette::NONE,
        }
    }
    pub fn height(y: f64) -> Self {
        Field::Height {
            x: 0.0,
            y,
            vignette: Vignette::NONE,
        }
    }

    /// A field angle off the meridional plane.
    pub fn angle_xy(x: f64, y: f64) -> Self {
        Field::Angle {
            x,
            y,
            vignette: Vignette::NONE,
        }
    }

    /// An object height off the meridional plane.
    pub fn height_xy(x: f64, y: f64) -> Self {
        Field::Height {
            x,
            y,
            vignette: Vignette::NONE,
        }
    }

    /// This field point with vignetting factors attached.
    pub fn vignetted(self, v: Vignette) -> Self {
        match self {
            Field::Angle { x, y, .. } => Field::Angle { x, y, vignette: v },
            Field::Height { x, y, .. } => Field::Height { x, y, vignette: v },
        }
    }

    /// This field point using the whole pupil.
    ///
    /// Anything that must not depend on how the pupil is sampled asks for this — the
    /// chief ray behind a distortion figure, most of all.
    pub fn unvignetted(self) -> Self {
        self.vignetted(Vignette::NONE)
    }

    /// The part of the pupil this field uses.
    pub fn vignette(&self) -> Vignette {
        match *self {
            Field::Angle { vignette, .. } | Field::Height { vignette, .. } => vignette,
        }
    }

    /// Distance of this field point from the axis, in its own units.
    pub fn radius(&self) -> f64 {
        match *self {
            Field::Angle { x, y, .. } | Field::Height { x, y, .. } => (x * x + y * y).sqrt(),
        }
    }

    /// Unit vector pointing from the axis toward this field point.
    ///
    /// A rotationally symmetric system images every field at the same distance from the
    /// axis, so [`radius`](Self::radius) is all the paraxial trace needs — but the image
    /// point also has a *direction*, and that is what tells +20 degrees from -20 degrees,
    /// and a field in X from the same field in Y. The axial field has no direction and
    /// reports `(0, 0)`, so anything scaled by it lands on the axis where it belongs.
    pub fn direction(&self) -> (f64, f64) {
        match *self {
            Field::Angle { x, y, .. } | Field::Height { x, y, .. } => {
                let r = (x * x + y * y).sqrt();
                if r == 0.0 {
                    (0.0, 0.0)
                } else {
                    (x / r, y / r)
                }
            }
        }
    }
}

/// How the system's aperture is specified.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(tag = "kind", rename_all = "snake_case"))]
pub enum Aperture {
    /// Diameter of the entrance pupil.
    EntrancePupilDiameter(f64),
    /// Working f-number in image space.
    ImageSpaceFNumber(f64),
    /// Numerical aperture in object space. For finite conjugates.
    ObjectSpaceNA(f64),
    /// Physical diameter of the stop surface.
    StopDiameter(f64),
}

/// A complete sequential optical system.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct System<S: Scalar> {
    pub title: String,
    pub object: Object<S>,
    /// First refracting surface through to the image plane, in order.
    pub surfaces: Vec<Surface<S>>,
    pub wavelengths: Vec<Wavelength>,
    pub fields: Vec<Field>,
    pub aperture: Aperture,
    /// The medium in object space, before the first surface.
    pub object_medium: Material,
}

impl<S: Scalar> System<S> {
    pub fn new(object: Object<S>, surfaces: Vec<Surface<S>>) -> Self {
        Self {
            title: String::new(),
            object,
            surfaces,
            wavelengths: vec![Wavelength::new(crate::material::lines::D)],
            fields: vec![Field::angle(0.0)],
            aperture: Aperture::EntrancePupilDiameter(1.0),
            object_medium: Material::Vacuum,
        }
    }

    pub fn titled(mut self, title: &str) -> Self {
        self.title = title.to_string();
        self
    }

    pub fn with_aperture(mut self, aperture: Aperture) -> Self {
        self.aperture = aperture;
        self
    }

    pub fn with_fields(mut self, fields: Vec<Field>) -> Self {
        self.fields = fields;
        self
    }

    pub fn with_wavelengths(mut self, wavelengths: Vec<Wavelength>) -> Self {
        self.wavelengths = wavelengths;
        self
    }

    /// Index of the last surface, which is the image plane.
    pub fn image_index(&self) -> usize {
        self.surfaces.len() - 1
    }

    /// Index of the aperture stop, if one is marked.
    pub fn stop_index(&self) -> Option<usize> {
        self.surfaces.iter().position(|s| s.is_stop)
    }

    /// The medium immediately before surface `i`.
    pub fn medium_before(&self, i: usize) -> &Material {
        if i == 0 {
            &self.object_medium
        } else {
            &self.surfaces[i - 1].material
        }
    }

    /// The medium immediately after surface `i`.
    pub fn medium_after(&self, i: usize) -> &Material {
        &self.surfaces[i].material
    }

    /// Axial position of each surface vertex, with surface 0 at `z = 0`.
    pub fn vertices(&self) -> Vec<S> {
        let mut z = S::zero();
        let mut out = Vec::with_capacity(self.surfaces.len());
        for s in &self.surfaces {
            out.push(z);
            z += s.thickness;
        }
        out
    }

    /// Placement of each surface in the global frame.
    ///
    /// Centred for now; coordinate breaks will compose into these without changing
    /// anything downstream, because the tracer already works through this indirection.
    pub fn transforms(&self) -> Vec<Transform<S>> {
        self.vertices()
            .into_iter()
            .map(Transform::along_axis)
            .collect()
    }

    /// Signed direction of propagation after `i` reflections. `+1` forward, `-1` folded.
    pub fn is_reflective_at(&self, i: usize) -> bool {
        self.surfaces[i].material.is_reflective()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::material::catalog;

    fn doublet() -> System<f64> {
        System::new(
            Object::Infinity,
            vec![
                Surface::new(50.0, 4.0, catalog::N_BK7),
                Surface::new(-30.0, 2.0, catalog::F2).stop(),
                Surface::new(-80.0, 90.0, Material::Vacuum),
                Surface::plane(0.0, Material::Vacuum).labelled("image"),
            ],
        )
    }

    #[test]
    fn vertices_accumulate_the_thicknesses() {
        assert_eq!(doublet().vertices(), vec![0.0, 4.0, 6.0, 96.0]);
    }

    #[test]
    fn transforms_place_each_vertex_on_the_axis() {
        let sys = doublet();
        let t = sys.transforms();
        assert_eq!(t.len(), sys.surfaces.len());
        for (x, z) in t.iter().zip(sys.vertices()) {
            assert_eq!(x.point_to_global(crate::math::Vec3::zero()).z, z);
            assert!(x.is_identity_rotation(), "centred systems must not rotate");
        }
    }

    #[test]
    fn a_surface_owns_the_medium_that_follows_it() {
        let sys = doublet();
        assert_eq!(sys.medium_before(0), &Material::Vacuum); // object space
        assert_eq!(sys.medium_after(0), &catalog::N_BK7);
        assert_eq!(sys.medium_before(1), &catalog::N_BK7);
        assert_eq!(sys.medium_after(1), &catalog::F2);
        assert_eq!(sys.medium_before(2), &catalog::F2);
        assert_eq!(sys.medium_after(2), &Material::Vacuum);
    }

    #[test]
    fn the_object_medium_is_configurable() {
        let mut sys = doublet();
        sys.object_medium = Material::Fixed { n: 1.33 };
        assert_eq!(sys.medium_before(0), &Material::Fixed { n: 1.33 });
    }

    #[test]
    fn the_stop_and_the_image_are_located_by_index() {
        let sys = doublet();
        assert_eq!(sys.stop_index(), Some(1));
        assert_eq!(sys.image_index(), 3);

        let mut no_stop = doublet();
        no_stop.surfaces[1].is_stop = false;
        assert_eq!(no_stop.stop_index(), None);
    }

    #[test]
    fn plane_surfaces_really_are_planar() {
        let s = Surface::<f64>::plane(3.0, Material::Vacuum);
        assert!(matches!(s.profile, crate::surface::Profile::Plane));
        assert_eq!(s.thickness, 3.0);
        assert!(!s.is_stop);
        assert_eq!(s.semi_diameter, None);
    }

    #[test]
    fn builders_are_chainable_and_do_not_disturb_each_other() {
        let s = Surface::<f64>::new(10.0, 1.0, Material::Vacuum)
            .stop()
            .labelled("aperture")
            .with_semi_diameter(4.5);
        assert!(s.is_stop);
        assert_eq!(s.label, "aperture");
        assert_eq!(s.semi_diameter, Some(4.5));

        let sys = doublet()
            .titled("test")
            .with_aperture(Aperture::ImageSpaceFNumber(2.8))
            .with_fields(vec![Field::angle(0.0), Field::angle(10.0)])
            .with_wavelengths(vec![Wavelength::new(0.6)]);
        assert_eq!(sys.title, "test");
        assert_eq!(sys.aperture, Aperture::ImageSpaceFNumber(2.8));
        assert_eq!(sys.fields.len(), 2);
        assert_eq!(sys.wavelengths.len(), 1);
    }

    #[test]
    fn a_fresh_system_has_sensible_defaults() {
        let sys = doublet();
        assert_eq!(sys.wavelengths.len(), 1);
        assert_eq!(sys.wavelengths[0].um, crate::material::lines::D);
        assert_eq!(sys.wavelengths[0].weight, 1.0);
        assert_eq!(sys.fields, vec![Field::angle(0.0)]);
        assert_eq!(sys.object_medium, Material::Vacuum);
    }

    #[test]
    fn field_radius_is_measured_from_the_axis() {
        assert_eq!(Field::angle_xy(3.0, 4.0).radius(), 5.0);
        assert_eq!(Field::height_xy(-3.0, 4.0).radius(), 5.0);
        assert_eq!(Field::angle(7.0).radius(), 7.0);
        assert_eq!(Field::height(0.0).radius(), 0.0);
    }

    #[test]
    fn the_full_pupil_is_the_default_and_changes_nothing() {
        let v = Vignette::default();
        assert!(v.is_full());
        assert_eq!(v.apply(0.7, -0.3), (0.7, -0.3));
        assert!(Field::angle(5.0).vignette().is_full());
    }

    #[test]
    fn compression_shrinks_the_pupil_and_decentre_moves_it() {
        let squeezed = Vignette {
            cy: 0.5,
            ..Vignette::NONE
        };
        assert_eq!(squeezed.apply(1.0, 1.0), (1.0, 0.5));

        // A decentred pupil takes its centre with it, which is the whole point: the ray
        // at (0, 0) is the one analyses call the chief ray.
        let shifted = Vignette {
            dy: 0.25,
            ..Vignette::NONE
        };
        assert_eq!(shifted.apply(0.0, 0.0), (0.0, 0.25));
        assert_eq!(shifted.apply(0.0, 1.0), (0.0, 1.25));
    }

    #[test]
    fn the_vignetting_angle_turns_the_reduced_pupil_not_the_full_one() {
        // Squeeze in y, then turn a quarter turn: the squeeze must end up in x.
        let v = Vignette {
            cy: 0.5,
            angle: 90.0,
            ..Vignette::NONE
        };
        let (x, y) = v.apply(0.0, 1.0);
        assert!((x + 0.5).abs() < 1e-15 && y.abs() < 1e-15, "({x}, {y})");
        let (x, y) = v.apply(1.0, 0.0);
        assert!(x.abs() < 1e-15 && (y - 1.0).abs() < 1e-15, "({x}, {y})");
    }

    #[test]
    fn a_field_carries_its_vignetting_and_can_be_stripped_of_it() {
        let v = Vignette {
            cx: 0.3,
            dy: -0.2,
            ..Vignette::NONE
        };
        let f = Field::angle(12.0).vignetted(v);
        assert_eq!(f.vignette(), v);
        assert_eq!(f.radius(), 12.0);
        assert_eq!(f.unvignetted(), Field::angle(12.0));
    }

    #[test]
    fn field_direction_separates_what_radius_merges() {
        assert_eq!(Field::angle_xy(3.0, 4.0).direction(), (0.6, 0.8));
        assert_eq!(Field::angle(7.0).direction(), (0.0, 1.0));
        assert_eq!(Field::angle(-7.0).direction(), (0.0, -1.0));
        assert_eq!(Field::angle_xy(7.0, 0.0).direction(), (1.0, 0.0));

        // The axis has no direction to report, and must not invent one.
        assert_eq!(Field::angle(0.0).direction(), (0.0, 0.0));
        assert_eq!(Field::height(0.0).direction(), (0.0, 0.0));
    }

    #[test]
    fn reflectivity_is_read_from_the_following_medium() {
        let mut sys = doublet();
        assert!(!sys.is_reflective_at(0));
        sys.surfaces[0].material = Material::Mirror;
        assert!(sys.is_reflective_at(0));
    }
}
