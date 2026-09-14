//! `optic-core` -- the sequential ray tracing kernel.
//!
//! The crate holds no I/O, no formatting and no UI. Everything here is a pure function
//! of a [`System`], which is what makes analyses cacheable, optimisation reproducible
//! and the whole thing testable without a running application.
//!
//! Everything geometric is generic over [`Scalar`]. Trace with `f64` for a plain result,
//! or with [`Dual`] to get exact derivatives with respect to design variables in the
//! same pass -- the basis for gradient-based optimisation and tolerance sensitivities.
//!
//! ```
//! use optic_core::{samples, Paraxial};
//!
//! let sys = samples::cooke_triplet::<f64>();
//! let par = Paraxial::compute(&sys, 0.5875618);
//! assert!((par.efl - 50.0).abs() < 0.5);
//! assert!(par.fno > 4.9 && par.fno < 5.1);
//! ```
//!
//! Lengths are millimetres by convention and wavelengths are micrometres, matching
//! every published dispersion formula.

#![deny(rust_2018_idioms)]
#![warn(missing_debug_implementations)]

pub mod material;
pub mod math;
pub mod paraxial;
pub mod samples;
pub mod surface;
pub mod system;
pub mod trace;

pub use material::{catalog, lines, Material};
pub use math::{Dual, Scalar, Transform, Vec3};
pub use paraxial::{Paraxial, ParaxialState};
pub use surface::{MissReason, Profile};
pub use system::{Aperture, Field, Object, Surface, System, Wavelength};
pub use trace::{launch, trace, Hit, Ray, TracedRay};
