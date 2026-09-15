//! Surface profiles: sag, normal, and ray intersection in the surface's local frame.
//!
//! Local frame convention: the vertex is at the origin and the axis is `+z`. Sag is a
//! function of `r^2 = x^2 + y^2` alone for every profile here, which keeps both the
//! normal and the aspheric Newton iteration to one derivative.

use crate::math::{Scalar, Vec3};

/// Why a ray failed to reach a surface.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MissReason {
    /// The ray passed outside the region where the conic is defined (beyond the
    /// sphere's equator, or outside a hyperboloid's valid radius).
    OutsideProfile,
    /// No real intersection, or only one travelling backwards.
    NoIntersection,
    /// The aspheric Newton iteration did not converge.
    NoConvergence,
    /// The ray struck outside the clear aperture.
    ClippedByAperture,
    /// Total internal reflection at a refracting surface.
    TotalInternalReflection,
}

/// The shape of a surface, in its own local frame.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(tag = "type", rename_all = "snake_case"))]
pub enum Profile<S: Scalar> {
    /// Flat. The degenerate case is worth its own variant: it is the most common
    /// surface in real systems and its intersection is exact and branch-free.
    Plane,
    /// Sphere or conic of revolution: `curvature` is `1/R` (zero for flat), `conic` is
    /// the Schwarzschild constant `k` (0 sphere, -1 paraboloid, `k < -1` hyperboloid).
    Conic { curvature: S, conic: S },
    /// A conic plus an even polynomial. `coeffs[i]` multiplies `r^(2i+4)`, so the list
    /// reads as the familiar A4, A6, A8, ... of a lens prescription.
    EvenAsphere {
        curvature: S,
        conic: S,
        coeffs: Vec<S>,
    },
}

impl<S: Scalar> Profile<S> {
    /// A spherical surface of the given radius. Radius 0 or infinite means flat.
    pub fn sphere(radius: f64) -> Self {
        if radius == 0.0 || !radius.is_finite() {
            Profile::Plane
        } else {
            Profile::Conic {
                curvature: S::from_f64(1.0 / radius),
                conic: S::zero(),
            }
        }
    }

    fn curvature(&self) -> S {
        match self {
            Profile::Plane => S::zero(),
            Profile::Conic { curvature, .. } | Profile::EvenAsphere { curvature, .. } => *curvature,
        }
    }

    fn conic_constant(&self) -> S {
        match self {
            Profile::Plane => S::zero(),
            Profile::Conic { conic, .. } | Profile::EvenAsphere { conic, .. } => *conic,
        }
    }

    /// Sag `z` at squared radial distance `r2`.
    pub fn sag(&self, r2: S) -> Result<S, MissReason> {
        match self {
            Profile::Plane => Ok(S::zero()),
            Profile::Conic { curvature, conic } => conic_sag(*curvature, *conic, r2),
            Profile::EvenAsphere {
                curvature,
                conic,
                coeffs,
            } => {
                let mut z = conic_sag(*curvature, *conic, r2)?;
                let mut p = r2 * r2; // r^4
                for a in coeffs {
                    z += *a * p;
                    p *= r2;
                }
                Ok(z)
            }
        }
    }

    /// `d(sag) / d(r2)`. Halved, this is the slope factor the normal needs.
    fn dsag_dr2(&self, r2: S) -> Result<S, MissReason> {
        match self {
            Profile::Plane => Ok(S::zero()),
            Profile::Conic { curvature, conic } => conic_dsag(*curvature, *conic, r2),
            Profile::EvenAsphere {
                curvature,
                conic,
                coeffs,
            } => {
                let mut d = conic_dsag(*curvature, *conic, r2)?;
                let mut p = r2; // r^(2i+2), so that d/d(r2) of r^(2i+4) is (i+2) r^(2i+2)
                for (i, a) in coeffs.iter().enumerate() {
                    d += *a * S::from_f64((i + 2) as f64) * p;
                    p *= r2;
                }
                Ok(d)
            }
        }
    }

    /// Outward unit normal at a point on the surface, pointing towards `+z`.
    pub fn normal(&self, p: Vec3<S>) -> Result<Vec3<S>, MissReason> {
        match self {
            Profile::Plane => Ok(Vec3::axis()),
            _ => {
                let r2 = p.x * p.x + p.y * p.y;
                let s = self.dsag_dr2(r2)?;
                // F(x,y,z) = z - sag(x^2+y^2); grad F = (-2x s', -2y s', 1).
                let two = S::two();
                Ok(Vec3::new(-two * p.x * s, -two * p.y * s, S::one()).normalized())
            }
        }
    }

    /// Distance along a unit-direction ray, in local coordinates, to this surface.
    ///
    /// Returns the root on the near branch of the surface: for a lens surface that is
    /// the one continuous with the vertex, never the far side of the sphere.
    pub fn intersect(&self, origin: Vec3<S>, dir: Vec3<S>) -> Result<S, MissReason> {
        let c = self.curvature();

        // Flat, or flat enough that the quadratic below degenerates.
        if matches!(self, Profile::Plane) || c.value() == 0.0 {
            if dir.z.value() == 0.0 {
                return Err(MissReason::NoIntersection);
            }
            let t = -origin.z / dir.z;
            return match self {
                Profile::EvenAsphere { .. } => self.newton(origin, dir, t),
                _ => Ok(t),
            };
        }

        let t = self.conic_root(origin, dir, c, self.conic_constant())?;
        match self {
            Profile::EvenAsphere { .. } => self.newton(origin, dir, t),
            _ => Ok(t),
        }
    }

    /// Closed-form root for `c(x^2 + y^2 + (1+k) z^2) - 2z = 0`.
    fn conic_root(&self, o: Vec3<S>, d: Vec3<S>, c: S, k: S) -> Result<S, MissReason> {
        let kk = S::one() + k;
        let a = c * (d.x * d.x + d.y * d.y + kk * d.z * d.z);
        let b = S::two() * (c * (o.x * d.x + o.y * d.y + kk * o.z * d.z) - d.z);
        let cc = c * (o.x * o.x + o.y * o.y + kk * o.z * o.z) - S::two() * o.z;

        let disc = b * b - S::from_f64(4.0) * a * cc;
        if disc.value() < 0.0 {
            return Err(MissReason::OutsideProfile);
        }
        // t = 2C / (-B + sqrt(disc)) is the root that stays finite as `a -> 0` and
        // selects the branch through the vertex for a forward-going ray. The naive
        // (-B +/- sqrt)/2a form loses precision for near-flat surfaces.
        let denom = -b + disc.sqrt();
        if denom.value().abs() < 1e-300 {
            return Err(MissReason::NoIntersection);
        }
        Ok(S::two() * cc / denom)
    }

    /// Refine an intersection onto the full aspheric profile.
    ///
    /// Derivatives are recovered by the implicit function theorem rather than by
    /// differentiating through the iteration: once the value has converged, `t` is reset
    /// to a constant and a single Newton step in dual arithmetic reproduces
    /// `dt/dp = -(dF/dp)/(dF/dt)` exactly. Iteration count therefore never perturbs the
    /// Jacobian the optimiser sees.
    fn newton(&self, o: Vec3<S>, d: Vec3<S>, seed: S) -> Result<S, MissReason> {
        const MAX_ITERS: usize = 32;
        const TOL: f64 = 1e-12;

        let mut t = seed;
        let mut converged = false;
        for _ in 0..MAX_ITERS {
            let (f, dfdt) = self.residual(o, d, t)?;
            if f.value().abs() < TOL {
                converged = true;
                break;
            }
            if dfdt.value().abs() < 1e-300 {
                return Err(MissReason::NoConvergence);
            }
            t -= f / dfdt;
        }
        if !converged {
            let (f, _) = self.residual(o, d, t)?;
            if f.value().abs() >= TOL {
                return Err(MissReason::NoConvergence);
            }
        }

        let t0 = S::from_f64(t.value());
        let (f, dfdt) = self.residual(o, d, t0)?;
        Ok(t0 - f / dfdt)
    }

    /// `F(t) = p_z(t) - sag(r^2(t))` and its derivative in `t`.
    fn residual(&self, o: Vec3<S>, d: Vec3<S>, t: S) -> Result<(S, S), MissReason> {
        let p = o + d * t;
        let r2 = p.x * p.x + p.y * p.y;
        let f = p.z - self.sag(r2)?;
        let dr2_dt = S::two() * (p.x * d.x + p.y * d.y);
        let dfdt = d.z - self.dsag_dr2(r2)? * dr2_dt;
        Ok((f, dfdt))
    }
}

fn conic_sag<S: Scalar>(c: S, k: S, r2: S) -> Result<S, MissReason> {
    if c.value() == 0.0 {
        return Ok(S::zero());
    }
    let w2 = S::one() - (S::one() + k) * c * c * r2;
    if w2.value() <= 0.0 {
        return Err(MissReason::OutsideProfile);
    }
    Ok(c * r2 / (S::one() + w2.sqrt()))
}

/// `d(sag)/d(r2) = c / (2 sqrt(1 - (1+k) c^2 r^2))`, from `dz/dr = c r / sqrt(...)`.
fn conic_dsag<S: Scalar>(c: S, k: S, r2: S) -> Result<S, MissReason> {
    if c.value() == 0.0 {
        return Ok(S::zero());
    }
    let w2 = S::one() - (S::one() + k) * c * c * r2;
    if w2.value() <= 0.0 {
        return Err(MissReason::OutsideProfile);
    }
    Ok(c / (S::two() * w2.sqrt()))
}

impl<S: Scalar> Profile<S> {
    /// Curvature at the vertex, which is all the paraxial trace sees.
    ///
    /// With aspheric terms starting at `r^4`, the polynomial contributes nothing to
    /// second order, so this is just the base curvature.
    pub fn paraxial_curvature(&self) -> S {
        self.curvature()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const R: f64 = 50.0;

    fn sphere() -> Profile<f64> {
        Profile::sphere(R)
    }

    #[test]
    fn a_degenerate_radius_means_flat() {
        assert!(matches!(Profile::<f64>::sphere(0.0), Profile::Plane));
        assert!(matches!(
            Profile::<f64>::sphere(f64::INFINITY),
            Profile::Plane
        ));
        assert!(matches!(Profile::<f64>::sphere(f64::NAN), Profile::Plane));
        assert!(matches!(
            Profile::<f64>::sphere(10.0),
            Profile::Conic { .. }
        ));
    }

    #[test]
    fn a_plane_is_flat_everywhere_and_faces_the_axis() {
        let p = Profile::<f64>::Plane;
        assert_eq!(p.sag(0.0).unwrap(), 0.0);
        assert_eq!(p.sag(100.0).unwrap(), 0.0);
        assert_eq!(p.paraxial_curvature(), 0.0);
        assert_eq!(
            p.normal(Vec3::new(3.0, 4.0, 0.0)).unwrap().value(),
            [0.0, 0.0, 1.0]
        );
    }

    #[test]
    fn spherical_sag_matches_the_closed_form() {
        // z = R - sqrt(R^2 - r^2)
        for r in [0.0, 1.0, 10.0, 25.0, 49.0] {
            let expected = R - (R * R - r * r).sqrt();
            let got = sphere().sag(r * r).unwrap();
            assert!((got - expected).abs() < 1e-12, "r={r}: {got} vs {expected}");
        }
    }

    #[test]
    fn a_paraboloid_has_exactly_quadratic_sag() {
        // k = -1 removes the higher-order terms of the conic entirely.
        let p: Profile<f64> = Profile::Conic {
            curvature: 1.0 / R,
            conic: -1.0,
        };
        for r in [0.5, 5.0, 40.0, 500.0] {
            let expected = r * r / (2.0 * R);
            assert!((p.sag(r * r).unwrap() - expected).abs() < 1e-10, "r={r}");
        }
    }

    #[test]
    fn a_sphere_is_undefined_beyond_its_equator() {
        assert_eq!(sphere().sag(R * R * 1.01), Err(MissReason::OutsideProfile));
        assert!(sphere().sag(R * R * 0.99).is_ok());
    }

    #[test]
    fn a_spherical_normal_points_at_the_centre_of_curvature() {
        // Vertex at the origin puts the centre at (0, 0, R).
        let centre = Vec3::new(0.0, 0.0, R);
        for (x, y) in [(0.0, 0.0), (5.0, 0.0), (-8.0, 12.0)] {
            let z = sphere().sag(x * x + y * y).unwrap();
            let p = Vec3::new(x, y, z);
            let n = sphere().normal(p).unwrap();
            assert!((n.norm() - 1.0).abs() < 1e-14);
            // The normal and the vector towards the centre must be parallel.
            assert!(n.cross(centre - p).norm() < 1e-12, "at ({x}, {y})");
        }
    }

    #[test]
    fn an_axial_ray_meets_the_vertex() {
        let t = sphere()
            .intersect(Vec3::new(0.0, 0.0, -10.0), Vec3::axis())
            .unwrap();
        assert!((t - 10.0).abs() < 1e-13);
    }

    #[test]
    fn intersections_land_on_the_sphere_itself() {
        let centre = Vec3::new(0.0, 0.0, R);
        for (oy, dy) in [(0.0, 0.05), (10.0, 0.0), (-20.0, 0.3)] {
            let o = Vec3::new(2.0, oy, -30.0);
            let d = Vec3::new(0.01, dy, 1.0).normalized();
            let t = sphere().intersect(o, d).unwrap();
            let p = o + d * t;
            assert!(
                ((p - centre).norm() - R).abs() < 1e-11,
                "off-sphere at oy={oy}"
            );
            // The near branch: the vertex side, not the far side of the sphere.
            assert!(p.z < R, "took the far root at oy={oy}");
        }
    }

    #[test]
    fn a_ray_parallel_to_a_plane_never_meets_it() {
        let r =
            Profile::<f64>::Plane.intersect(Vec3::new(0.0, 1.0, -5.0), Vec3::new(0.0, 1.0, 0.0));
        assert_eq!(r, Err(MissReason::NoIntersection));
    }

    #[test]
    fn a_ray_that_misses_the_sphere_entirely_is_rejected() {
        let r = sphere().intersect(Vec3::new(0.0, 500.0, -30.0), Vec3::axis());
        assert_eq!(r, Err(MissReason::OutsideProfile));
    }

    #[test]
    fn an_asphere_with_no_polynomial_terms_is_its_base_conic() {
        let conic: Profile<f64> = Profile::Conic {
            curvature: 1.0 / 30.0,
            conic: -0.5,
        };
        let asphere: Profile<f64> = Profile::EvenAsphere {
            curvature: 1.0 / 30.0,
            conic: -0.5,
            coeffs: vec![],
        };
        let o = Vec3::new(1.0, -2.0, -12.0);
        let d = Vec3::new(0.03, 0.02, 1.0).normalized();
        assert!((conic.intersect(o, d).unwrap() - asphere.intersect(o, d).unwrap()).abs() < 1e-12);
        assert!((conic.sag(9.0).unwrap() - asphere.sag(9.0).unwrap()).abs() < 1e-15);
    }

    #[test]
    fn aspheric_coefficients_read_as_a4_a6_a8() {
        // coeffs[i] multiplies r^(2i+4), so a flat base plus a single term is pure r^4.
        let p: Profile<f64> = Profile::EvenAsphere {
            curvature: 0.0,
            conic: 0.0,
            coeffs: vec![1e-4],
        };
        let r: f64 = 3.0;
        assert!((p.sag(r * r).unwrap() - 1e-4 * r.powi(4)).abs() < 1e-15);

        let q: Profile<f64> = Profile::EvenAsphere {
            curvature: 0.0,
            conic: 0.0,
            coeffs: vec![0.0, 2e-6],
        };
        assert!((q.sag(r * r).unwrap() - 2e-6 * r.powi(6)).abs() < 1e-18);
    }

    #[test]
    fn a_flat_based_asphere_still_intersects() {
        // Zero curvature takes the plane branch for the seed; Newton must still refine
        // onto the polynomial.
        let p: Profile<f64> = Profile::EvenAsphere {
            curvature: 0.0,
            conic: 0.0,
            coeffs: vec![1e-4],
        };
        let o = Vec3::new(0.0, 2.0, -5.0);
        let d = Vec3::new(0.0, 0.1, 1.0).normalized();
        let t = p.intersect(o, d).unwrap();
        let hit = o + d * t;
        let residual = hit.z - p.sag(hit.x * hit.x + hit.y * hit.y).unwrap();
        assert!(residual.abs() < 1e-12, "off-surface by {residual}");
    }

    #[test]
    fn paraxial_curvature_ignores_the_polynomial() {
        let p: Profile<f64> = Profile::EvenAsphere {
            curvature: 0.02,
            conic: -3.0,
            coeffs: vec![1e-3, 1e-5],
        };
        assert_eq!(p.paraxial_curvature(), 0.02);
        assert_eq!(Profile::<f64>::Plane.paraxial_curvature(), 0.0);
    }
}
