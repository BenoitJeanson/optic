//! The scalar abstraction that makes the whole tracer differentiable.
//!
//! Every geometric routine in this crate is generic over [`Scalar`]. Instantiating it
//! with `f64` gives an ordinary ray trace; instantiating it with [`Dual`](super::Dual)
//! gives the trace *and* exact derivatives of every output with respect to up to `N`
//! design variables, in a single pass.
//!
//! We use forward-mode AD deliberately. Lens design has few variables (curvatures,
//! thicknesses, conics, aspheric coefficients: tens) and many residuals (ray errors
//! across pupil, field and wavelength: thousands). Levenberg-Marquardt needs the full
//! `m x n` Jacobian; forward mode produces all of it in `~n` trace-passes, whereas
//! reverse mode would cost `~m`.

use core::fmt::Debug;
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// A real number that arithmetic in this crate is generic over.
pub trait Scalar:
    Copy
    + Clone
    + Debug
    + PartialEq
    + PartialOrd
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
    + AddAssign
    + SubAssign
    + MulAssign
    + DivAssign
{
    /// Lift a constant. The derivative part, if any, is zero.
    fn from_f64(v: f64) -> Self;

    /// The real part, discarding any derivative information.
    fn value(self) -> f64;

    fn sqrt(self) -> Self;
    fn abs(self) -> Self;
    fn sin(self) -> Self;
    fn cos(self) -> Self;
    fn tan(self) -> Self;
    fn asin(self) -> Self;
    fn acos(self) -> Self;
    fn atan2(self, x: Self) -> Self;
    fn powi(self, n: i32) -> Self;

    #[inline]
    fn zero() -> Self {
        Self::from_f64(0.0)
    }
    #[inline]
    fn one() -> Self {
        Self::from_f64(1.0)
    }
    #[inline]
    fn two() -> Self {
        Self::from_f64(2.0)
    }
    #[inline]
    fn recip(self) -> Self {
        Self::one() / self
    }
    #[inline]
    fn square(self) -> Self {
        self * self
    }
    #[inline]
    fn is_finite(self) -> bool {
        self.value().is_finite()
    }
    #[inline]
    fn min_of(self, other: Self) -> Self {
        if self < other {
            self
        } else {
            other
        }
    }
    #[inline]
    fn max_of(self, other: Self) -> Self {
        if self > other {
            self
        } else {
            other
        }
    }
}

impl Scalar for f64 {
    #[inline]
    fn from_f64(v: f64) -> Self {
        v
    }
    #[inline]
    fn value(self) -> f64 {
        self
    }
    #[inline]
    fn sqrt(self) -> Self {
        f64::sqrt(self)
    }
    #[inline]
    fn abs(self) -> Self {
        f64::abs(self)
    }
    #[inline]
    fn sin(self) -> Self {
        f64::sin(self)
    }
    #[inline]
    fn cos(self) -> Self {
        f64::cos(self)
    }
    #[inline]
    fn tan(self) -> Self {
        f64::tan(self)
    }
    #[inline]
    fn asin(self) -> Self {
        f64::asin(self)
    }
    #[inline]
    fn acos(self) -> Self {
        f64::acos(self)
    }
    #[inline]
    fn atan2(self, x: Self) -> Self {
        f64::atan2(self, x)
    }
    #[inline]
    fn powi(self, n: i32) -> Self {
        f64::powi(self, n)
    }
}
