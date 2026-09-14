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
