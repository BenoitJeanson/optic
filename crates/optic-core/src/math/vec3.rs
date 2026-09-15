//! Three-vectors, generic over the scalar type.

use super::Scalar;
use core::ops::{Add, Div, Mul, Neg, Sub};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vec3<S: Scalar> {
    pub x: S,
    pub y: S,
    pub z: S,
}

impl<S: Scalar> Vec3<S> {
    #[inline]
    pub fn new(x: S, y: S, z: S) -> Self {
        Self { x, y, z }
    }

    #[inline]
    pub fn splat(v: S) -> Self {
        Self { x: v, y: v, z: v }
    }

    #[inline]
    pub fn zero() -> Self {
        Self::splat(S::zero())
    }

    /// Unit vector along `+z`, the optical axis.
    #[inline]
    pub fn axis() -> Self {
        Self::new(S::zero(), S::zero(), S::one())
    }

    #[inline]
    pub fn from_f64(x: f64, y: f64, z: f64) -> Self {
        Self::new(S::from_f64(x), S::from_f64(y), S::from_f64(z))
    }

    #[inline]
    pub fn dot(self, o: Self) -> S {
        self.x * o.x + self.y * o.y + self.z * o.z
    }

    #[inline]
    pub fn cross(self, o: Self) -> Self {
        Self::new(
            self.y * o.z - self.z * o.y,
            self.z * o.x - self.x * o.z,
            self.x * o.y - self.y * o.x,
        )
    }

    #[inline]
    pub fn norm_squared(self) -> S {
        self.dot(self)
    }

    #[inline]
    pub fn norm(self) -> S {
        self.norm_squared().sqrt()
    }

    /// Radial distance from the optical axis.
    #[inline]
    pub fn radius(self) -> S {
        (self.x * self.x + self.y * self.y).sqrt()
    }

    #[inline]
    pub fn normalized(self) -> Self {
        self * self.norm().recip()
    }

    /// Drop derivative information. Used where a value-only decision is needed.
    #[inline]
    pub fn value(self) -> [f64; 3] {
        [self.x.value(), self.y.value(), self.z.value()]
    }
}

impl<S: Scalar> Add for Vec3<S> {
    type Output = Self;
    #[inline]
    fn add(self, o: Self) -> Self {
        Self::new(self.x + o.x, self.y + o.y, self.z + o.z)
    }
}

impl<S: Scalar> Sub for Vec3<S> {
    type Output = Self;
    #[inline]
    fn sub(self, o: Self) -> Self {
        Self::new(self.x - o.x, self.y - o.y, self.z - o.z)
    }
}

impl<S: Scalar> Neg for Vec3<S> {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        Self::new(-self.x, -self.y, -self.z)
    }
}

impl<S: Scalar> Mul<S> for Vec3<S> {
    type Output = Self;
    #[inline]
    fn mul(self, k: S) -> Self {
        Self::new(self.x * k, self.y * k, self.z * k)
    }
}

impl<S: Scalar> Div<S> for Vec3<S> {
    type Output = Self;
    #[inline]
    fn div(self, k: S) -> Self {
        Self::new(self.x / k, self.y / k, self.z / k)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::Dual;

    fn v(x: f64, y: f64, z: f64) -> Vec3<f64> {
        Vec3::new(x, y, z)
    }

    #[test]
    fn dot_and_cross_obey_their_identities() {
        let (a, b) = (v(1.0, 2.0, 3.0), v(-4.0, 5.0, 6.0));
        assert_eq!(a.dot(b), -4.0 + 10.0 + 18.0);
        // Cross product is orthogonal to both and anticommutative.
        let c = a.cross(b);
        assert!(c.dot(a).abs() < 1e-14 && c.dot(b).abs() < 1e-14);
        let back = b.cross(a);
        assert!((c + back).norm() < 1e-14);
    }

    #[test]
    fn parallel_vectors_have_no_cross_product() {
        let a = v(2.0, -1.0, 0.5);
        assert!(a.cross(a * 3.0).norm() < 1e-15);
    }

    #[test]
    fn normalisation_produces_unit_length() {
        let n = v(3.0, -4.0, 12.0).normalized();
        assert!((n.norm() - 1.0).abs() < 1e-15);
        assert!((n.x - 3.0 / 13.0).abs() < 1e-15);
    }

    #[test]
    fn norm_and_radius_measure_different_things() {
        let a = v(3.0, 4.0, 10.0);
        assert!((a.norm() - (125.0f64).sqrt()).abs() < 1e-14);
        assert!((a.radius() - 5.0).abs() < 1e-15); // radius ignores z
        assert!((a.norm_squared() - 125.0).abs() < 1e-13);
    }

    #[test]
    fn arithmetic_is_componentwise() {
        let (a, b) = (v(1.0, 2.0, 3.0), v(4.0, 5.0, 6.0));
        assert_eq!((a + b).value(), [5.0, 7.0, 9.0]);
        assert_eq!((b - a).value(), [3.0, 3.0, 3.0]);
        assert_eq!((a * 2.0).value(), [2.0, 4.0, 6.0]);
        assert_eq!((a / 2.0).value(), [0.5, 1.0, 1.5]);
        assert_eq!((-a).value(), [-1.0, -2.0, -3.0]);
    }

    #[test]
    fn the_axis_is_the_z_unit_vector() {
        assert_eq!(Vec3::<f64>::axis().value(), [0.0, 0.0, 1.0]);
        assert_eq!(Vec3::<f64>::zero().value(), [0.0, 0.0, 0.0]);
        assert_eq!(
            Vec3::<f64>::from_f64(1.0, 2.0, 3.0).value(),
            [1.0, 2.0, 3.0]
        );
    }

    #[test]
    fn vectors_carry_derivatives_through_their_operations() {
        // d/dx of |(x, 3, 4)| at x = 0 is 0, and of (x, 3, 4).dot((1,0,0)) is 1.
        let x = Dual::<1>::variable(0.0, 0);
        let a = Vec3::new(x, Dual::constant(3.0), Dual::constant(4.0));
        assert!((a.norm().re - 5.0).abs() < 1e-14);
        assert!(a.norm().grad()[0].abs() < 1e-14);
        assert_eq!(a.dot(Vec3::axis()).grad()[0], 0.0);
    }
}
