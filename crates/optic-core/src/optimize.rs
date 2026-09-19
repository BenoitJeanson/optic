//! Damped least-squares optimisation of a design.
//!
//! Lens design is a least-squares problem with an awkward shape: tens of variables and
//! often thousands of residuals, because every traced ray at every field and wavelength
//! contributes one. Levenberg--Marquardt is the standard answer, and it needs the full
//! `m x n` Jacobian rather than a single gradient.
//!
//! That shape is why the kernel is generic over [`Scalar`] and why the derivatives are
//! forward-mode. One trace with [`Dual<N>`] in place of `f64` yields every residual
//! *and* its derivative with respect to N variables at once, exactly — no finite
//! differences, no step-size to choose, and no separate adjoint pass to keep in step
//! with the tracer.
//!
//! `N` is a compile-time constant while the number of variables is not, so the Jacobian
//! is built in column blocks of [`BLOCK`]: a design with 14 variables costs two traces
//! per residual set rather than 14 finite-difference traces, and each column is exact.

use crate::math::{Dual, Scalar};
use crate::paraxial::Paraxial;
use crate::solve;
use crate::surface::Profile;
use crate::system::{Object, Surface, System};
use crate::trace::{launch, trace};

/// Variables differentiated per trace. The Jacobian is assembled in blocks this wide.
///
/// Larger blocks mean fewer traces but a wider dual number in every arithmetic
/// operation, most of whose lanes are idle on the last block. Eight is a reasonable
/// middle for designs with a handful to a few dozen variables.
pub const BLOCK: usize = 8;

/// A design parameter the optimiser may change.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(tag = "kind", rename_all = "snake_case"))]
pub enum Variable {
    /// Curvature of a surface — the reciprocal of its radius.
    ///
    /// Curvature, not radius, because a flat surface has curvature zero, an ordinary
    /// number the optimiser can step through in either direction, while its radius is
    /// infinite and cannot be stepped at all. A design that wants to pass through flat
    /// on its way from convex to concave is completely routine, so the variable has to
    /// admit it. Editors still show the radius.
    Curvature { surface: usize },
    /// Axial distance from this surface to the next.
    Thickness { surface: usize },
    /// Conic constant of a surface.
    Conic { surface: usize },
}

impl Variable {
    /// Current value of this parameter.
    pub fn get(&self, sys: &System<f64>) -> f64 {
        match *self {
            Variable::Curvature { surface } => sys.surfaces[surface].profile.curvature(),
            Variable::Thickness { surface } => sys.surfaces[surface].thickness,
            Variable::Conic { surface } => match sys.surfaces[surface].profile {
                Profile::Conic { conic, .. } | Profile::EvenAsphere { conic, .. } => conic,
                Profile::Plane => 0.0,
            },
        }
    }

    /// Write a new value into the design.
    pub fn set(&self, sys: &mut System<f64>, v: f64) {
        let s = &mut sys.surfaces[self.surface()];
        match *self {
            Variable::Curvature { .. } => match &mut s.profile {
                Profile::Conic { curvature, .. } | Profile::EvenAsphere { curvature, .. } => {
                    *curvature = v
                }
                // A plane is a conic of zero curvature; giving it one makes it a sphere.
                Profile::Plane => {
                    s.profile = Profile::Conic {
                        curvature: v,
                        conic: 0.0,
                    }
                }
            },
            Variable::Thickness { .. } => s.thickness = v,
            Variable::Conic { .. } => match &mut s.profile {
                Profile::Conic { conic, .. } | Profile::EvenAsphere { conic, .. } => *conic = v,
                Profile::Plane => {
                    s.profile = Profile::Conic {
                        curvature: 0.0,
                        conic: v,
                    }
                }
            },
        }
    }

    fn surface(&self) -> usize {
        match *self {
            Variable::Curvature { surface }
            | Variable::Thickness { surface }
            | Variable::Conic { surface } => surface,
        }
    }
}

/// One term of the merit function.
///
/// An operand produces residuals that the optimiser drives toward zero. Keeping them as
/// residuals rather than pre-summing into a single cost is what lets Levenberg--Marquardt
/// see the structure of the problem: the Gauss-Newton step needs each ray separately.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(tag = "op", rename_all = "snake_case"))]
pub enum Operand {
    /// Effective focal length at `wavelength`, driven to `target`.
    Efl {
        wavelength: usize,
        target: f64,
        weight: f64,
    },
    /// Transverse distance from the chief ray at the image plane, driven to zero.
    ///
    /// This is the operand a spot-size merit function is built from. It contributes two
    /// residuals, x and y: a spot is two-dimensional, and collapsing it to a radius
    /// first would hide from the optimiser which way to move.
    TransverseError {
        field: usize,
        wavelength: usize,
        /// Normalised pupil coordinates, each in `[-1, 1]`.
        px: f64,
        py: f64,
        weight: f64,
    },
    /// Penalise a thickness that leaves `[min, max]`, and nothing inside it.
    ///
    /// A one-sided residual: zero while the constraint holds, so an inactive bound adds
    /// nothing to the merit function and nothing to the Jacobian.
    ThicknessBounds {
        surface: usize,
        min: f64,
        max: f64,
        weight: f64,
    },
}

impl Operand {
    /// Append this operand's residuals to `out`.
    fn residuals<S: Scalar>(&self, sys: &System<S>, out: &mut Vec<S>) {
        match *self {
            Operand::Efl {
                wavelength,
                target,
                weight,
            } => {
                let Some(wl) = sys.wavelengths.get(wavelength) else {
                    return;
                };
                let efl = Paraxial::compute(sys, wl.um).efl;
                out.push((efl - S::from_f64(target)) * S::from_f64(weight));
            }

            Operand::TransverseError {
                field,
                wavelength,
                px,
                py,
                weight,
            } => {
                let (Some(&f), Some(wl)) = (sys.fields.get(field), sys.wavelengths.get(wavelength))
                else {
                    return;
                };
                let par = Paraxial::compute(sys, wl.um);
                // Measured against the chief ray rather than the axis, so the operand
                // asks for a tight spot and stays silent about where the spot sits.
                // Distortion is a separate question and belongs to a separate operand.
                let Some(reference) =
                    trace(sys, wl.um, launch(sys, &par, f.unvignetted(), 0.0, 0.0)).image_point()
                else {
                    return;
                };
                let Some(p) = trace(sys, wl.um, launch(sys, &par, f, px, py)).image_point() else {
                    // A blocked ray contributes nothing. It cannot contribute a large
                    // residual either: that would reward the optimiser for clipping the
                    // rays it cannot correct.
                    return;
                };
                let w = S::from_f64(weight);
                out.push((p.x - reference.x) * w);
                out.push((p.y - reference.y) * w);
            }

            Operand::ThicknessBounds {
                surface,
                min,
                max,
                weight,
            } => {
                let Some(s) = sys.surfaces.get(surface) else {
                    return;
                };
                let t = s.thickness;
                let w = S::from_f64(weight);
                let below = S::from_f64(min) - t;
                let above = t - S::from_f64(max);
                out.push(if below.value() > 0.0 {
                    below * w
                } else {
                    S::zero()
                });
                out.push(if above.value() > 0.0 {
                    above * w
                } else {
                    S::zero()
                });
            }
        }
    }
}

/// A set of operands, evaluated together.
#[derive(Clone, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Merit {
    pub operands: Vec<Operand>,
}

impl Merit {
    pub fn new(operands: Vec<Operand>) -> Self {
        Self { operands }
    }

    /// Every residual, in operand order.
    pub fn residuals<S: Scalar>(&self, sys: &System<S>) -> Vec<S> {
        let mut out = Vec::new();
        for op in &self.operands {
            op.residuals(sys, &mut out);
        }
        out
    }

    /// Root-mean-square residual: the number a report shows.
    ///
    /// The optimiser itself works with the sum of squares, but that grows with the
    /// number of operands, so it says nothing on its own. The RMS is comparable between
    /// runs and between designs.
    pub fn rms(&self, sys: &System<f64>) -> f64 {
        let r = self.residuals(sys);
        if r.is_empty() {
            return 0.0;
        }
        (r.iter().map(|v| v * v).sum::<f64>() / r.len() as f64).sqrt()
    }

    /// Spot-size operands over every field, wavelength and a hexapolar pupil sample,
    /// together with a focal length to hold.
    ///
    /// This is the merit function a designer writes first, and the one the tests use.
    pub fn rms_spot(sys: &System<f64>, rings: usize, efl_target: Option<f64>) -> Self {
        let mut operands = Vec::new();
        if let Some(target) = efl_target {
            // Weighted well above the ray errors: a focal length is a specification,
            // not a preference, and millimetres of focal error would otherwise be traded
            // away against micrometres of blur.
            operands.push(Operand::Efl {
                wavelength: sys.wavelengths.len() / 2,
                target,
                weight: 100.0,
            });
        }
        for field in 0..sys.fields.len() {
            for wavelength in 0..sys.wavelengths.len() {
                for (px, py) in hexapolar(rings) {
                    operands.push(Operand::TransverseError {
                        field,
                        wavelength,
                        px,
                        py,
                        weight: 1.0,
                    });
                }
            }
        }
        Self { operands }
    }
}

/// Pupil sample points on `rings` rings, uniform in area.
fn hexapolar(rings: usize) -> Vec<(f64, f64)> {
    let mut pts = vec![(0.0, 0.0)];
    for r in 1..=rings {
        let radius = r as f64 / rings as f64;
        let count = 6 * r;
        for k in 0..count {
            let theta = std::f64::consts::TAU * k as f64 / count as f64;
            pts.push((radius * theta.cos(), radius * theta.sin()));
        }
    }
    pts
}

/// Rebuild `sys` with `vars[offset .. offset + N]` seeded as differentiation variables.
fn lift<const N: usize>(sys: &System<f64>, vars: &[Variable], offset: usize) -> System<Dual<N>> {
    type D<const N: usize> = Dual<N>;

    // Which lane, if any, each variable occupies in this block.
    let lane = |v: Variable| -> Option<usize> {
        vars.iter()
            .skip(offset)
            .take(N)
            .position(|candidate| *candidate == v)
    };

    let surfaces: Vec<Surface<D<N>>> = sys
        .surfaces
        .iter()
        .enumerate()
        .map(|(i, s)| {
            let seed = |v: Variable, value: f64| match lane(v) {
                Some(k) => D::<N>::variable(value, k),
                None => D::<N>::from_f64(value),
            };
            let curvature = seed(Variable::Curvature { surface: i }, s.profile.curvature());
            let conic_value = match s.profile {
                Profile::Conic { conic, .. } | Profile::EvenAsphere { conic, .. } => conic,
                Profile::Plane => 0.0,
            };
            let conic = seed(Variable::Conic { surface: i }, conic_value);

            let profile = match &s.profile {
                Profile::EvenAsphere { coeffs, .. } => Profile::EvenAsphere {
                    curvature,
                    conic,
                    coeffs: coeffs.iter().map(|c| D::<N>::from_f64(*c)).collect(),
                },
                // A plane is carried as a conic so that giving it curvature during
                // optimisation is a change of value rather than a change of type.
                _ => Profile::Conic { curvature, conic },
            };

            Surface {
                profile,
                thickness: seed(Variable::Thickness { surface: i }, s.thickness),
                thickness_solve: s.thickness_solve.lift(),
                material: s.material.clone(),
                semi_diameter: s.semi_diameter.map(D::<N>::from_f64),
                is_stop: s.is_stop,
                label: s.label.clone(),
            }
        })
        .collect();

    let object = match sys.object {
        Object::Infinity => Object::Infinity,
        Object::Finite { distance } => Object::Finite {
            distance: D::<N>::from_f64(distance),
        },
    };

    let mut out = System::new(object, surfaces);
    out.title = sys.title.clone();
    out.aperture = sys.aperture;
    out.fields = sys.fields.clone();
    out.wavelengths = sys.wavelengths.clone();
    out.object_medium = sys.object_medium.clone();
    out
}

/// A copy of `sys` with its thickness solves satisfied.
///
/// Residuals are always measured on a resolved design, because a solve is part of the
/// prescription rather than a post-processing step: a system whose image plane has not
/// been refocused is not the design anyone is asking about.
pub fn resolved(sys: &System<f64>) -> System<f64> {
    let mut out = sys.clone();
    let _ = solve::resolve(&mut out, primary(sys));
    out
}

/// Residuals and the Jacobian of those residuals with respect to `vars`.
///
/// Returns `(r, jacobian)` where `jacobian[i][j]` is `d r_i / d vars[j]`.
///
/// Solves are re-resolved on the lifted system, so the derivatives flow *through* them:
/// the reported sensitivity of a spot to a curvature includes the refocus that an
/// autofocus solve performs in response. Holding the image plane still instead would
/// linearise a design nobody is proposing to build.
pub fn jacobian(sys: &System<f64>, vars: &[Variable], merit: &Merit) -> (Vec<f64>, Vec<Vec<f64>>) {
    let r = merit.residuals(&resolved(sys));
    let m = r.len();
    let mut jac = vec![vec![0.0; vars.len()]; m];

    let mut offset = 0;
    while offset < vars.len() {
        let mut lifted = lift::<BLOCK>(sys, vars, offset);
        let _ = solve::resolve(&mut lifted, primary(sys));
        let block = merit.residuals(&lifted);
        // A residual that vanished under duals means a ray that traced differently, which
        // it must not: the two systems hold the same numbers.
        debug_assert_eq!(
            block.len(),
            m,
            "the lifted system produced a different residual set"
        );

        for (i, value) in block.iter().enumerate() {
            let grad = value.grad();
            for (k, g) in grad.iter().enumerate() {
                if offset + k < vars.len() {
                    jac[i][offset + k] = *g;
                }
            }
        }
        offset += BLOCK;
    }
    (r, jac)
}

/// How an optimisation run ended.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(tag = "status", rename_all = "snake_case"))]
pub enum Stop {
    /// The merit function stopped improving.
    Converged,
    /// The step became too small to matter.
    StepTooSmall,
    /// The iteration limit was reached with the merit still falling.
    IterationLimit,
    /// No variable, no operand, or a system that will not trace.
    NothingToDo,
}

/// What an optimisation run did.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Report {
    pub stop: Stop,
    pub iterations: usize,
    /// RMS residual before and after.
    pub before: f64,
    pub after: f64,
}

/// Settings for [`optimise`].
#[derive(Clone, Copy, Debug)]
pub struct Settings {
    pub max_iterations: usize,
    /// Stop when one iteration improves the sum of squares by less than this fraction.
    pub tolerance: f64,
    /// Initial Levenberg--Marquardt damping.
    pub damping: f64,
}

impl Default for Settings {
    fn default() -> Self {
        Self {
            max_iterations: 60,
            tolerance: 1e-10,
            damping: 1e-3,
        }
    }
}

/// Improve `sys` in place by damped least squares.
///
/// The damping is what makes this usable rather than merely correct. A pure Gauss-Newton
/// step is excellent near a minimum and wild far from one, which is the normal state of
/// a design being optimised; raising the damping on a rejected step slides continuously
/// toward gradient descent, and lowering it on an accepted one restores the fast
/// convergence once the step lands.
pub fn optimise(
    sys: &mut System<f64>,
    vars: &[Variable],
    merit: &Merit,
    settings: Settings,
) -> Report {
    *sys = resolved(sys);
    let before = merit.rms(sys);
    if vars.is_empty() || merit.operands.is_empty() {
        return Report {
            stop: Stop::NothingToDo,
            iterations: 0,
            before,
            after: before,
        };
    }

    let n = vars.len();
    let mut lambda = settings.damping;
    let mut cost = sum_of_squares(&merit.residuals(sys));
    let mut stop = Stop::IterationLimit;
    let mut iterations = 0;

    for _ in 0..settings.max_iterations {
        iterations += 1;
        let (r, jac) = jacobian(sys, vars, merit);
        if r.is_empty() {
            stop = Stop::NothingToDo;
            break;
        }

        // Normal equations: (JtJ + lambda diag(JtJ)) dx = -Jt r.
        let mut jtj = vec![vec![0.0; n]; n];
        let mut jtr = vec![0.0; n];
        for (i, row) in jac.iter().enumerate() {
            for a in 0..n {
                jtr[a] += row[a] * r[i];
                for b in a..n {
                    jtj[a][b] += row[a] * row[b];
                }
            }
        }
        // Mirror the upper triangle into the lower one. Indices are the subject here --
        // the point is precisely that element (a, b) equals element (b, a) -- so an
        // iterator form would obscure the only thing the loop says.
        #[allow(clippy::needless_range_loop)]
        for a in 0..n {
            for b in 0..a {
                jtj[a][b] = jtj[b][a];
            }
        }

        let mut accepted = false;
        // Try progressively heavier damping before giving up on this iteration.
        for _ in 0..12 {
            let mut a = jtj.clone();
            for k in 0..n {
                // Scaling the damping by the diagonal keeps it meaningful when variables
                // have wildly different units, which curvatures and thicknesses do.
                let d = if jtj[k][k] > 0.0 { jtj[k][k] } else { 1.0 };
                a[k][k] += lambda * d;
            }
            let rhs: Vec<f64> = jtr.iter().map(|v| -v).collect();
            let Some(step) = solve_spd(a, rhs) else {
                lambda *= 10.0;
                continue;
            };

            let mut candidate = sys.clone();
            for (v, dx) in vars.iter().zip(&step) {
                v.set(&mut candidate, v.get(sys) + dx);
            }
            // Solves are part of the design, so they are re-applied before the candidate
            // is judged: otherwise a step is scored against a system that has not
            // refocused, and every step looks worse than it is.
            let _ = solve::resolve(&mut candidate, primary(sys));

            let trial = sum_of_squares(&merit.residuals(&candidate));
            if trial.is_finite() && trial < cost {
                let improvement = (cost - trial) / cost.max(f64::MIN_POSITIVE);
                *sys = candidate;
                cost = trial;
                lambda = (lambda * 0.1).max(1e-12);
                accepted = true;
                if improvement < settings.tolerance {
                    stop = Stop::Converged;
                }
                break;
            }
            lambda *= 10.0;
        }

        if !accepted {
            stop = Stop::StepTooSmall;
            break;
        }
        if stop == Stop::Converged {
            break;
        }
    }

    Report {
        stop,
        iterations,
        before,
        after: merit.rms(sys),
    }
}

fn primary(sys: &System<f64>) -> f64 {
    sys.wavelengths
        .get(sys.wavelengths.len() / 2)
        .map(|w| w.um)
        .unwrap_or(crate::material::lines::D)
}

fn sum_of_squares(r: &[f64]) -> f64 {
    r.iter().map(|v| v * v).sum()
}

/// Solve `a x = b` for a symmetric positive-definite `a`, by Cholesky.
///
/// Returns `None` when `a` is not positive definite, which is the caller's signal to
/// damp harder rather than to take a meaningless step.
// Cholesky reads across two rows at once (`a[i][j]` against `a[k][j]`) while writing a
// third, which no iterator form expresses without index arithmetic of its own.
#[allow(clippy::needless_range_loop)]
fn solve_spd(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    for k in 0..n {
        let mut d = a[k][k];
        for j in 0..k {
            d -= a[k][j] * a[k][j];
        }
        if d <= 0.0 || !d.is_finite() {
            return None;
        }
        let d = d.sqrt();
        a[k][k] = d;
        for i in k + 1..n {
            let mut s = a[i][k];
            for j in 0..k {
                s -= a[i][j] * a[k][j];
            }
            a[i][k] = s / d;
        }
    }
    // Forward substitution, then back substitution.
    for i in 0..n {
        let mut s = b[i];
        for j in 0..i {
            s -= a[i][j] * b[j];
        }
        b[i] = s / a[i][i];
    }
    for i in (0..n).rev() {
        let mut s = b[i];
        for j in i + 1..n {
            s -= a[j][i] * b[j];
        }
        b[i] = s / a[i][i];
    }
    if b.iter().all(|v| v.is_finite()) {
        Some(b)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::material::{lines, Material};
    use crate::system::{Aperture, Field, Surface, Wavelength};

    fn triplet() -> System<f64> {
        crate::samples::smith_triplet_moderate::<f64>()
    }

    fn curvatures(sys: &System<f64>) -> Vec<Variable> {
        (0..sys.surfaces.len() - 1)
            .filter(|i| sys.surfaces[*i].profile.curvature() != 0.0)
            .map(|surface| Variable::Curvature { surface })
            .collect()
    }

    #[test]
    fn the_jacobian_matches_central_differences() {
        // The one test that decides whether anything else here means anything. Finite
        // differences are what the analytic Jacobian exists to replace, so they are the
        // only independent check available: agreement to this precision cannot happen by
        // accident across a whole matrix.
        let sys = triplet();
        // Curvatures and the unsolved thicknesses together: more variables than one
        // block holds, so the column-blocking loop runs more than once, and of mixed
        // kinds, so a lane is not silently reading the wrong parameter.
        let mut vars = curvatures(&sys);
        vars.extend(
            (0..sys.surfaces.len() - 2)
                .filter(|i| sys.surfaces[*i].thickness_solve.is_fixed())
                .map(|surface| Variable::Thickness { surface }),
        );
        assert!(vars.len() > BLOCK, "only {} variables", vars.len());

        let merit = Merit::new(vec![
            Operand::Efl {
                wavelength: 1,
                target: 100.0,
                weight: 1.0,
            },
            Operand::TransverseError {
                field: 2,
                wavelength: 1,
                px: 0.0,
                py: 0.7,
                weight: 1.0,
            },
        ]);

        let (_, jac) = jacobian(&sys, &vars, &merit);

        for (j, v) in vars.iter().enumerate() {
            let h = 1e-7;
            let mut up = sys.clone();
            v.set(&mut up, v.get(&sys) + h);
            let mut down = sys.clone();
            v.set(&mut down, v.get(&sys) - h);

            // Resolved on both sides, like the analytic path: the autofocus solve moves
            // when a curvature does, and a derivative that ignored that would be the
            // derivative of a different system.
            let ru = merit.residuals(&resolved(&up));
            let rd = merit.residuals(&resolved(&down));
            for i in 0..ru.len() {
                let fd = (ru[i] - rd[i]) / (2.0 * h);
                let exact = jac[i][j];
                let scale = exact.abs().max(fd.abs()).max(1.0);
                assert!(
                    (fd - exact).abs() / scale < 1e-5,
                    "d r{i} / d v{j}: {exact} exact against {fd} by differences"
                );
            }
        }
    }

    #[test]
    fn a_single_variable_hits_an_exact_focal_length_target() {
        // One variable, one operand, and a residual that can reach zero: there is a right
        // answer and the optimiser must find it, not merely improve on the start.
        let mut sys = triplet();
        let merit = Merit::new(vec![Operand::Efl {
            wavelength: 1,
            target: 95.0,
            weight: 1.0,
        }]);
        let vars = [Variable::Curvature { surface: 0 }];

        let report = optimise(&mut sys, &vars, &merit, Settings::default());
        let efl = Paraxial::compute(&sys, sys.wavelengths[1].um).efl;

        assert!(
            (efl - 95.0).abs() < 1e-9,
            "{efl} after {} iterations, {report:?}",
            report.iterations
        );
        assert!(report.after < report.before);
    }

    #[test]
    fn optimising_a_mirror_finds_the_parabola() {
        // A concave mirror images a point at infinity without spherical aberration when
        // its conic constant is exactly -1. The answer is a closed form the optimiser
        // knows nothing about, so arriving at it is evidence the machinery works rather
        // than evidence it is self-consistent.
        let mut sys = System::new(
            Object::Infinity,
            vec![
                Surface::new(-200.0, -100.0, Material::Mirror).stop(),
                Surface::plane(0.0, Material::Vacuum).labelled("image"),
            ],
        )
        .with_aperture(Aperture::EntrancePupilDiameter(40.0))
        .with_fields(vec![Field::angle(0.0)])
        .with_wavelengths(vec![Wavelength::new(lines::D)]);

        let merit = Merit::new(
            [0.3, 0.6, 0.9]
                .iter()
                .map(|py| Operand::TransverseError {
                    field: 0,
                    wavelength: 0,
                    px: 0.0,
                    py: *py,
                    weight: 1.0,
                })
                .collect(),
        );
        let vars = [Variable::Conic { surface: 0 }];

        let report = optimise(&mut sys, &vars, &merit, Settings::default());
        let conic = Variable::Conic { surface: 0 }.get(&sys);

        assert!(
            (conic + 1.0).abs() < 1e-6,
            "conic {conic} after {} iterations, {report:?}",
            report.iterations
        );
    }

    #[test]
    fn optimisation_never_leaves_a_design_worse_than_it_found_it() {
        // Every step is accepted only if it lowers the sum of squares, so this holds
        // whatever the starting point -- including one the optimiser cannot improve.
        let mut sys = triplet();
        let merit = Merit::rms_spot(&sys, 2, Some(100.0));
        let vars = curvatures(&sys);

        let report = optimise(&mut sys, &vars, &merit, Settings::default());
        assert!(
            report.after <= report.before,
            "{} became {}",
            report.before,
            report.after
        );
    }

    #[test]
    fn a_badly_spoiled_design_is_recovered() {
        // Bend the published triplet 5% out of shape -- an RMS residual of 32 against
        // its own 2 -- and ask for it back.
        //
        // The optimiser reaches a *lower* residual than the published design itself.
        // That is not a claim to have improved on a 1919 patent: it says this merit
        // function is narrower than the problem its designer actually solved, who was
        // also balancing Petzval curvature, secondary spectrum, the glasses available to
        // him and the cost of making the thing. What the test establishes is the part
        // that is being tested -- that the machinery descends a long way from a bad
        // start, and lands on a lens of the same family rather than on something
        // numerically clever and physically absurd.
        let published = triplet();
        let merit = Merit::rms_spot(&published, 2, Some(100.0));
        let vars = curvatures(&published);

        let mut sys = published.clone();
        for v in &vars {
            v.set(&mut sys, v.get(&published) * 1.05);
        }
        let spoiled = merit.rms(&sys);
        assert!(
            spoiled > merit.rms(&published) * 10.0,
            "the perturbation barely moved it: {spoiled}"
        );

        let report = optimise(&mut sys, &vars, &merit, Settings::default());
        assert!(
            report.after < merit.rms(&published),
            "recovered to {}, published sits at {}, {report:?}",
            report.after,
            merit.rms(&published)
        );

        // Same lens family: no surface turned itself inside out to get there.
        for v in &vars {
            assert_eq!(
                v.get(&sys).signum(),
                v.get(&published).signum(),
                "{v:?} changed sign: {} became {}",
                v.get(&published),
                v.get(&sys)
            );
        }
    }

    #[test]
    fn a_bound_is_silent_until_it_is_crossed() {
        let sys = triplet();
        let inside = Merit::new(vec![Operand::ThicknessBounds {
            surface: 1,
            min: 0.0,
            max: 100.0,
            weight: 1.0,
        }]);
        assert_eq!(inside.rms(&sys), 0.0);

        let outside = Merit::new(vec![Operand::ThicknessBounds {
            surface: 1,
            min: 50.0,
            max: 100.0,
            weight: 1.0,
        }]);
        assert!(outside.rms(&sys) > 0.0);
    }

    #[test]
    fn nothing_to_optimise_is_reported_rather_than_iterated() {
        let mut sys = triplet();
        let merit = Merit::rms_spot(&sys, 1, None);
        let report = optimise(&mut sys, &[], &merit, Settings::default());
        assert_eq!(report.stop, Stop::NothingToDo);
        assert_eq!(report.iterations, 0);
    }
}
