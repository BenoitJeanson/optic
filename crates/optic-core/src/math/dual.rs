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

#[cfg(test)]
mod tests {
    use super::*;

    /// Derivatives of `f` at `x`, compared against a hand-written derivative.
    fn check(x0: f64, f: impl Fn(Dual<1>) -> Dual<1>, value: f64, deriv: f64) {
        let y = f(Dual::<1>::variable(x0, 0));
        assert!((y.re - value).abs() < 1e-12, "value {} vs {value}", y.re);
        assert!(
            (y.grad()[0] - deriv).abs() < 1e-11,
            "derivative {} vs {deriv}",
            y.grad()[0]
        );
    }

    #[test]
    fn a_constant_has_no_derivative() {
        let c = Dual::<3>::constant(4.0);
        assert_eq!(c.re, 4.0);
        assert_eq!(c.du, [0.0; 3]);
    }

    #[test]
    fn a_variable_is_seeded_in_exactly_one_slot() {
        let v = Dual::<3>::variable(2.5, 1);
        assert_eq!(v.re, 2.5);
        assert_eq!(v.du, [0.0, 1.0, 0.0]);
    }

    #[test]
    #[should_panic(expected = "out of range")]
    fn seeding_past_the_end_is_rejected() {
        let _ = Dual::<2>::variable(1.0, 2);
    }

    #[test]
    fn sum_and_difference_are_linear() {
        let (x, y) = (Dual::<2>::variable(3.0, 0), Dual::<2>::variable(5.0, 1));
        assert_eq!((x + y).du, [1.0, 1.0]);
        assert_eq!((x - y).du, [1.0, -1.0]);
        assert_eq!((x + y).re, 8.0);
        assert_eq!((-x).du, [-1.0, 0.0]);
    }

    #[test]
    fn product_rule() {
        let (x, y) = (Dual::<2>::variable(3.0, 0), Dual::<2>::variable(5.0, 1));
        let p = x * y;
        assert_eq!(p.re, 15.0);
        assert_eq!(p.du, [5.0, 3.0]); // d(xy)/dx = y, d(xy)/dy = x
    }

    #[test]
    fn quotient_rule() {
        let (x, y) = (Dual::<2>::variable(3.0, 0), Dual::<2>::variable(5.0, 1));
        let q = x / y;
        assert!((q.re - 0.6).abs() < 1e-15);
        assert!((q.du[0] - 1.0 / 5.0).abs() < 1e-15);
        assert!((q.du[1] + 3.0 / 25.0).abs() < 1e-15);
    }

    #[test]
    fn assignment_operators_agree_with_their_binary_forms() {
        let (x, y) = (Dual::<2>::variable(1.5, 0), Dual::<2>::variable(-0.5, 1));

        let mut t = x;
        t += y;
        assert_eq!((t.re, t.du), ((x + y).re, (x + y).du));

        let mut t = x;
        t -= y;
        assert_eq!((t.re, t.du), ((x - y).re, (x - y).du));

        let mut t = x;
        t *= y;
        assert_eq!((t.re, t.du), ((x * y).re, (x * y).du));

        let mut t = x;
        t /= y;
        assert_eq!((t.re, t.du), ((x / y).re, (x / y).du));
    }

    #[test]
    fn ordering_looks_only_at_the_real_part() {
        // Two values that are equal but depend on different variables must compare
        // equal: every branch in the tracer is a question about magnitude.
        let a = Dual::<2>::variable(1.0, 0);
        let b = Dual::<2>::variable(1.0, 1);
        assert_eq!(a, b);
        assert_eq!(a.partial_cmp(&b), Some(core::cmp::Ordering::Equal));
        assert!(Dual::<2>::constant(0.5) < a);
    }

    #[test]
    fn transcendentals_match_their_derivatives() {
        let x0 = 0.7;
        check(x0, |x| x.sqrt(), x0.sqrt(), 0.5 / x0.sqrt());
        check(x0, |x| x.sin(), x0.sin(), x0.cos());
        check(x0, |x| x.cos(), x0.cos(), -x0.sin());
        check(x0, |x| x.tan(), x0.tan(), 1.0 + x0.tan().powi(2));
        check(x0, |x| x.asin(), x0.asin(), 1.0 / (1.0 - x0 * x0).sqrt());
        check(x0, |x| x.acos(), x0.acos(), -1.0 / (1.0 - x0 * x0).sqrt());
        check(x0, |x| x.powi(3), x0.powi(3), 3.0 * x0 * x0);
        check(x0, |x| x.powi(-2), x0.powi(-2), -2.0 * x0.powi(-3));
    }

    #[test]
    fn powi_zero_is_a_constant() {
        let y = Dual::<1>::variable(3.0, 0).powi(0);
        assert_eq!(y.re, 1.0);
        assert_eq!(y.du, [0.0]);
    }

    #[test]
    fn abs_flips_the_derivative_below_zero() {
        assert_eq!(Dual::<1>::variable(-2.0, 0).abs().du, [-1.0]);
        assert_eq!(Dual::<1>::variable(2.0, 0).abs().du, [1.0]);
    }

    #[test]
    fn atan2_differentiates_in_both_arguments() {
        let (y0, x0) = (1.3, -2.1);
        let y = Dual::<2>::variable(y0, 0);
        let x = Dual::<2>::variable(x0, 1);
        let a = y.atan2(x);
        let r2 = x0 * x0 + y0 * y0;
        assert!((a.re - y0.atan2(x0)).abs() < 1e-15);
        assert!((a.du[0] - x0 / r2).abs() < 1e-14);
        assert!((a.du[1] + y0 / r2).abs() < 1e-14);
    }

    #[test]
    fn derivatives_compose_through_several_operations() {
        // f(x) = sqrt(x) sin(x) / (1 + x^2)
        let x0 = 1.7;
        let (s, c, r) = (x0.sin(), x0.cos(), x0.sqrt());
        let d = 1.0 + x0 * x0;
        check(
            x0,
            |x| x.sqrt() * x.sin() / (Dual::constant(1.0) + x.powi(2)),
            r * s / d,
            ((0.5 / r * s + r * c) * d - r * s * 2.0 * x0) / (d * d),
        );
    }

    #[test]
    fn independent_variables_stay_independent() {
        // A function of only x must report exactly zero sensitivity to y.
        let x = Dual::<2>::variable(2.0, 0);
        let y = x.sin() * x.sqrt();
        assert_eq!(y.du[1], 0.0);
        assert!(y.du[0] != 0.0);
    }
}
