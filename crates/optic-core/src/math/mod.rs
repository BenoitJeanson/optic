//! Scalar-generic numerics: the layer that makes the tracer differentiable.

mod dual;
mod scalar;
mod transform;
mod vec3;

pub use dual::Dual;
pub use scalar::Scalar;
pub use transform::Transform;
pub use vec3::Vec3;
