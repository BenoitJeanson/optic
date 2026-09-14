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
    pub thickness: S,
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

/// One field point.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(tag = "kind", rename_all = "snake_case"))]
pub enum Field {
    /// Incidence angle in degrees. Only meaningful for an object at infinity.
    Angle { x: f64, y: f64 },
    /// Object height in lens units. Only meaningful for a finite object.
    Height { x: f64, y: f64 },
}

impl Field {
    pub fn angle(y: f64) -> Self {
        Field::Angle { x: 0.0, y }
    }
    pub fn height(y: f64) -> Self {
        Field::Height { x: 0.0, y }
    }

    /// Distance of this field point from the axis, in its own units.
    pub fn radius(&self) -> f64 {
        match *self {
            Field::Angle { x, y } | Field::Height { x, y } => (x * x + y * y).sqrt(),
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
