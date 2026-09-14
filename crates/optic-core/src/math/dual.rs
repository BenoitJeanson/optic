//! Forward-mode automatic differentiation over `N` simultaneous design variables.
//!
//! These are the hottest loops in the crate -- every arithmetic operation in a trace
//! passes through one. They are written as indexed loops over a const-sized array on
//! purpose: it is the form LLVM unrolls and vectorises most reliably, and the iterator
//! equivalent reads far worse for no gain.
#![allow(clippy::needless_range_loop)]

use super::Scalar;
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// A real number carrying its partial derivatives with respect to `N` variables.
///
/// Comparison operators look at [`Dual::re`] only: ordering a value against its own
/// derivative makes no sense, and every branch in the tracer (which root, which side of
/// an aperture, has Newton converged) is a decision about magnitudes.
#[derive(Clone, Copy, Debug)]
pub struct Dual<const N: usize> {
    /// Real part.
    pub re: f64,
    /// Partial derivatives `d(re)/d(v_i)`.
    pub du: [f64; N],
}

impl<const N: usize> Dual<N> {
    /// A constant: no dependence on any variable.
    #[inline]
    pub const fn constant(re: f64) -> Self {
        Self { re, du: [0.0; N] }
    }

    /// Variable `i`, seeded with unit derivative.
    ///
    /// # Panics
    /// If `i >= N`.
    #[inline]
    pub fn variable(re: f64, i: usize) -> Self {
        assert!(i < N, "variable index {i} out of range for Dual<{N}>");
        let mut du = [0.0; N];
        du[i] = 1.0;
        Self { re, du }
    }

    /// Apply a scalar function given its value and its derivative at `self.re`.
    ///
    /// This is the single place the chain rule lives; every transcendental below is
    /// one line on top of it.
    #[inline]
    pub fn chain(self, val: f64, deriv: f64) -> Self {
        let mut du = [0.0; N];
        for i in 0..N {
            du[i] = deriv * self.du[i];
        }
        Self { re: val, du }
    }

    /// The gradient with respect to all `N` variables.
    #[inline]
    pub fn grad(&self) -> &[f64; N] {
        &self.du
    }
}

impl<const N: usize> PartialEq for Dual<N> {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.re == other.re
    }
}

impl<const N: usize> PartialOrd for Dual<N> {
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
        self.re.partial_cmp(&other.re)
    }
}

impl<const N: usize> Add for Dual<N> {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        let mut du = [0.0; N];
        for i in 0..N {
            du[i] = self.du[i] + rhs.du[i];
        }
        Self {
            re: self.re + rhs.re,
            du,
        }
    }
}

impl<const N: usize> Sub for Dual<N> {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        let mut du = [0.0; N];
        for i in 0..N {
            du[i] = self.du[i] - rhs.du[i];
        }
        Self {
            re: self.re - rhs.re,
            du,
        }
    }
}

impl<const N: usize> Mul for Dual<N> {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Self) -> Self {
        let mut du = [0.0; N];
        for i in 0..N {
            du[i] = self.re * rhs.du[i] + self.du[i] * rhs.re;
        }
        Self {
            re: self.re * rhs.re,
            du,
        }
    }
}

impl<const N: usize> Div for Dual<N> {
    type Output = Self;
    #[inline]
    fn div(self, rhs: Self) -> Self {
        let inv = 1.0 / rhs.re;
        let inv2 = inv * inv;
        let mut du = [0.0; N];
        for i in 0..N {
            du[i] = (self.du[i] * rhs.re - self.re * rhs.du[i]) * inv2;
        }
        Self {
            re: self.re * inv,
            du,
        }
    }
}

impl<const N: usize> Neg for Dual<N> {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        let mut du = [0.0; N];
        for i in 0..N {
            du[i] = -self.du[i];
        }
        Self { re: -self.re, du }
    }
}

macro_rules! assign_op {
    ($trait:ident, $method:ident, $op:tt) => {
        impl<const N: usize> $trait for Dual<N> {
            #[inline]
            fn $method(&mut self, rhs: Self) {
                *self = *self $op rhs;
            }
        }
    };
}
assign_op!(AddAssign, add_assign, +);
assign_op!(SubAssign, sub_assign, -);
assign_op!(MulAssign, mul_assign, *);
assign_op!(DivAssign, div_assign, /);

impl<const N: usize> Scalar for Dual<N> {
    #[inline]
    fn from_f64(v: f64) -> Self {
        Self::constant(v)
    }
    #[inline]
    fn value(self) -> f64 {
        self.re
    }
    #[inline]
    fn sqrt(self) -> Self {
        let s = self.re.sqrt();
        self.chain(s, 0.5 / s)
    }
    #[inline]
    fn abs(self) -> Self {
        if self.re < 0.0 {
            -self
        } else {
            self
        }
    }
    #[inline]
    fn sin(self) -> Self {
        self.chain(self.re.sin(), self.re.cos())
    }
    #[inline]
    fn cos(self) -> Self {
        self.chain(self.re.cos(), -self.re.sin())
    }
    #[inline]
    fn tan(self) -> Self {
        let t = self.re.tan();
        self.chain(t, 1.0 + t * t)
    }
    #[inline]
    fn asin(self) -> Self {
        self.chain(self.re.asin(), 1.0 / (1.0 - self.re * self.re).sqrt())
    }
    #[inline]
    fn acos(self) -> Self {
        self.chain(self.re.acos(), -1.0 / (1.0 - self.re * self.re).sqrt())
    }
    #[inline]
    fn atan2(self, x: Self) -> Self {
        let denom = self.re * self.re + x.re * x.re;
        let mut du = [0.0; N];
        for i in 0..N {
            du[i] = (x.re * self.du[i] - self.re * x.du[i]) / denom;
        }
        Self {
            re: self.re.atan2(x.re),
            du,
        }
    }
    #[inline]
    fn powi(self, n: i32) -> Self {
        if n == 0 {
            return Self::constant(1.0);
        }
        self.chain(self.re.powi(n), n as f64 * self.re.powi(n - 1))
    }
}
