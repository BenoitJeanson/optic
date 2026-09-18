//! Solves: values the system computes for itself rather than the designer typing them.
//!
//! A solve turns a thickness into a *constraint* — "put the image where the marginal ray
//! crosses the axis", "keep this spacing equal to that one" — so the design stays
//! self-consistent as other parameters move. It is the difference between refocusing by
//! hand after every edit and never being out of focus at all.
//!
//! Solves are evaluated before any ray is traced, so everything downstream (analysis,
//! optimisation, tolerancing) sees a system that already satisfies them.
//!
//! Because the arithmetic is done in the generic scalar type, derivatives flow through a
//! solve exactly as they do through anything else: the optimiser sees the true
//! sensitivity of a solved thickness to the variables that drive it, with no special
//! handling.

use crate::math::Scalar;
use crate::system::System;

/// How a surface's thickness is determined.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(tag = "type", rename_all = "snake_case"))]
pub enum ThicknessSolve<S: Scalar> {
    /// The designer's value, used as written.
    Fixed,
    /// Set the thickness so the paraxial marginal ray reaches `height` at the next
    /// surface. With `height = 0` on the last thickness this is autofocus: the image
    /// plane sits where the axial ray crosses the axis.
    MarginalRayHeight { height: S },
    /// Set the thickness so the paraxial chief ray reaches `height` at the next surface.
    /// Placing a stop at the chief ray's axis crossing is the usual use.
    ChiefRayHeight { height: S },
    /// Copy another surface's thickness: `scale * thickness[from] + offset`.
    Pickup { from: usize, scale: S, offset: S },
}

// Written out rather than derived: `#[derive(Default)]` would add an `S: Default`
// bound, and `Dual<N>` deliberately does not implement `Default` -- a derivative with
// no variable attached to it is a bug waiting to happen, not a sensible zero value.
#[allow(clippy::derivable_impls)]
impl<S: Scalar> Default for ThicknessSolve<S> {
    fn default() -> Self {
        ThicknessSolve::Fixed
    }
}

impl<S: Scalar> ThicknessSolve<S> {
    pub fn is_fixed(&self) -> bool {
        matches!(self, ThicknessSolve::Fixed)
    }

    /// A short label for the editor, in the style of a prescription table.
    pub fn code(&self) -> &'static str {
        match self {
            ThicknessSolve::Fixed => "",
            ThicknessSolve::MarginalRayHeight { .. } => "M",
            ThicknessSolve::ChiefRayHeight { .. } => "C",
            ThicknessSolve::Pickup { .. } => "P",
        }
    }
}

/// Outcome of evaluating a system's solves.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SolveReport {
    /// No surface carries a solve.
    None,
    /// Every solve settled, in this many passes.
    Converged { passes: usize },
    /// Ran out of passes. The thicknesses are the best available, not a fixed point.
    NotConverged { passes: usize },
    /// A solve could not be evaluated; the thickness was left as it was.
    Degenerate {
        surface: usize,
        reason: &'static str,
    },
}

impl SolveReport {
    pub fn is_ok(&self) -> bool {
        !matches!(
            self,
            SolveReport::NotConverged { .. } | SolveReport::Degenerate { .. }
        )
    }
}

/// Evaluate every solve in `sys`, rewriting the thicknesses they control.
///
/// Solves can depend on each other through the paraxial rays they are defined against —
/// a solve before the stop moves the entrance pupil, which moves the rays, which moves
/// the solve. Rather than order them into a dependency graph, which is wrong anyway when
/// the dependency is genuinely circular, this iterates to a fixed point. Systems in
/// practice settle in two or three passes; the cap is there for the ones that cannot.
pub fn resolve<S: Scalar>(sys: &mut System<S>, wavelength: f64) -> SolveReport {
    if sys.surfaces.iter().all(|s| s.thickness_solve.is_fixed()) {
        return SolveReport::None;
    }

    const MAX_PASSES: usize = 12;
    const TOL: f64 = 1e-12;

    for pass in 1..=MAX_PASSES {
        let before: Vec<f64> = sys.surfaces.iter().map(|s| s.thickness.value()).collect();

        // Recomputed each pass: the pupil these rays are referred to moves as solves
        // change the geometry ahead of it.
        let par = crate::paraxial::Paraxial::compute(sys, wavelength);
        let widest = sys
            .fields
            .iter()
            .copied()
            .max_by(|a, b| a.radius().total_cmp(&b.radius()))
            .unwrap_or(crate::system::Field::angle(0.0));
        let chief = crate::paraxial::chief_ray(sys, wavelength, widest, par.ep_z);

        // Indexed rather than iterated: the body writes back into `sys.surfaces`, which
        // it is simultaneously reading a solve out of.
        #[allow(clippy::needless_range_loop)]
        for i in 0..sys.surfaces.len() {
            let solve = sys.surfaces[i].thickness_solve;
            let solved = match solve {
                ThicknessSolve::Fixed => continue,
                ThicknessSolve::MarginalRayHeight { height } => {
                    let state = par.marginal[i];
                    if state.u.value().abs() < 1e-300 {
                        return SolveReport::Degenerate {
                            surface: i,
                            reason: "the marginal ray is parallel to the axis here",
                        };
                    }
                    (height - state.y) / state.u
                }
                ThicknessSolve::ChiefRayHeight { height } => {
                    let state = chief[i];
                    if state.u.value().abs() < 1e-300 {
                        return SolveReport::Degenerate {
                            surface: i,
                            reason: "the chief ray is parallel to the axis here",
                        };
                    }
                    (height - state.y) / state.u
                }
                ThicknessSolve::Pickup {
                    from,
                    scale,
                    offset,
                } => {
                    if from >= sys.surfaces.len() {
                        return SolveReport::Degenerate {
                            surface: i,
                            reason: "the pickup refers to a surface that does not exist",
                        };
                    }
                    if from == i {
                        return SolveReport::Degenerate {
                            surface: i,
                            reason: "the pickup refers to itself",
                        };
                    }
                    sys.surfaces[from].thickness * scale + offset
                }
            };
            sys.surfaces[i].thickness = solved;
        }

        let settled = sys
            .surfaces
            .iter()
            .zip(&before)
            .all(|(s, b)| (s.thickness.value() - b).abs() <= TOL * b.abs().max(1.0));
        if settled {
            return SolveReport::Converged { passes: pass };
        }
    }
    SolveReport::NotConverged { passes: MAX_PASSES }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::material::lines;
    use crate::math::Dual;
    use crate::samples;
    use crate::surface::Profile;

    fn cooke_autofocused<S: Scalar>() -> System<S> {
        let mut sys = samples::cooke_triplet::<S>();
        let last = sys.image_index() - 1;
        sys.surfaces[last].thickness_solve =
            ThicknessSolve::MarginalRayHeight { height: S::zero() };
        sys
    }

    #[test]
    fn a_system_without_solves_reports_nothing_to_do() {
        let mut sys = samples::cooke_triplet::<f64>();
        // The sample ships with an autofocus solve, which is the point of it.
        assert!(!sys.surfaces[sys.image_index() - 1]
            .thickness_solve
            .is_fixed());

        for s in &mut sys.surfaces {
            s.thickness_solve = ThicknessSolve::Fixed;
        }
        assert_eq!(resolve(&mut sys, lines::D), SolveReport::None);
    }

    #[test]
    fn autofocus_agrees_with_computing_the_back_focal_distance() {
        // Two independent routes to the same plane: the solve marches the real marginal
        // ray, while `focus` uses the first-order back focal distance. On an infinite
        // conjugate they must land in the same place.
        let manual = samples::cooke_triplet::<f64>();
        let expected = manual.surfaces[manual.image_index() - 1].thickness;

        let mut solved = cooke_autofocused::<f64>();
        let last = solved.image_index() - 1;
        solved.surfaces[last].thickness = 1.0; // deliberately wrong
        let report = resolve(&mut solved, lines::D);

        assert!(report.is_ok(), "{report:?}");
        let got = solved.surfaces[last].thickness;
        assert!((got - expected).abs() < 1e-10, "{got} vs {expected}");
    }

    #[test]
    fn a_solve_keeps_the_image_in_focus_when_the_design_changes() {
        // This is the whole point of a solve: edit a radius and stay focused.
        let mut sys = cooke_autofocused::<f64>();
        if let Profile::Conic { curvature, .. } = &mut sys.surfaces[0].profile {
            *curvature *= 1.15;
        }
        resolve(&mut sys, lines::D);

        let par = crate::Paraxial::compute(&sys, lines::D);
        let marginal = par.marginal[sys.image_index()];
        assert!(
            marginal.y.abs() < 1e-10,
            "the marginal ray misses the axis at the image by {}",
            marginal.y
        );
    }

    #[test]
    fn solves_converge_quickly() {
        let mut sys = cooke_autofocused::<f64>();
        let last = sys.image_index() - 1;
        sys.surfaces[last].thickness = 500.0;
        match resolve(&mut sys, lines::D) {
            SolveReport::Converged { passes } => assert!(passes <= 3, "took {passes} passes"),
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn a_chief_ray_solve_lands_on_the_axis_crossing() {
        let mut sys = samples::cooke_triplet::<f64>();
        // Make the gap before the stop place the stop at the chief ray crossing.
        sys.surfaces[3].thickness_solve = ThicknessSolve::ChiefRayHeight { height: 0.0 };
        assert!(resolve(&mut sys, lines::D).is_ok());

        let par = crate::Paraxial::compute(&sys, lines::D);
        let chief =
            crate::paraxial::chief_ray(&sys, lines::D, crate::system::Field::angle(20.0), par.ep_z);
        assert!(
            chief[4].y.abs() < 1e-9,
            "chief ray height at surface 5: {}",
            chief[4].y
        );
    }

    #[test]
    fn a_pickup_tracks_the_surface_it_copies() {
        let mut sys = samples::cooke_triplet::<f64>();
        sys.surfaces[3].thickness_solve = ThicknessSolve::Pickup {
            from: 1,
            scale: 0.5,
            offset: 0.25,
        };
        assert!(resolve(&mut sys, lines::D).is_ok());
        assert!(
            (sys.surfaces[3].thickness - (sys.surfaces[1].thickness * 0.5 + 0.25)).abs() < 1e-12
        );

        // And it follows when the source moves.
        sys.surfaces[1].thickness = 10.0;
        assert!(resolve(&mut sys, lines::D).is_ok());
        assert!((sys.surfaces[3].thickness - 5.25).abs() < 1e-12);
    }

    #[test]
    fn nonsense_solves_are_reported_rather_than_applied() {
        let mut sys = samples::cooke_triplet::<f64>();
        sys.surfaces[2].thickness_solve = ThicknessSolve::Pickup {
            from: 2,
            scale: 1.0,
            offset: 0.0,
        };
        assert!(matches!(
            resolve(&mut sys, lines::D),
            SolveReport::Degenerate { surface: 2, .. }
        ));

        let mut sys = samples::cooke_triplet::<f64>();
        sys.surfaces[2].thickness_solve = ThicknessSolve::Pickup {
            from: 99,
            scale: 1.0,
            offset: 0.0,
        };
        assert!(matches!(
            resolve(&mut sys, lines::D),
            SolveReport::Degenerate { .. }
        ));
    }

    #[test]
    fn derivatives_flow_through_a_solved_thickness() {
        // A solved thickness is a function of the design, so the optimiser must see its
        // sensitivity. Compare the dual-number derivative against a central difference.
        let build = |dc: Dual<1>| {
            let mut sys = cooke_autofocused::<Dual<1>>();
            if let Profile::Conic { curvature, .. } = &mut sys.surfaces[0].profile {
                *curvature += dc;
            }
            resolve(&mut sys, lines::D);
            sys.surfaces[sys.image_index() - 1].thickness
        };
        let ad = build(Dual::<1>::variable(0.0, 0)).grad()[0];

        let value_at = |dc: f64| {
            let mut sys = cooke_autofocused::<f64>();
            if let Profile::Conic { curvature, .. } = &mut sys.surfaces[0].profile {
                *curvature += dc;
            }
            resolve(&mut sys, lines::D);
            sys.surfaces[sys.image_index() - 1].thickness
        };
        let h = 1e-7;
        let fd = (value_at(h) - value_at(-h)) / (2.0 * h);

        assert!(ad.abs() > 1.0, "the derivative is suspiciously small: {ad}");
        assert!(
            (ad - fd).abs() / fd.abs() < 1e-5,
            "d(thickness)/d(curvature): AD {ad} vs FD {fd}"
        );
    }

    #[test]
    fn solve_codes_are_stable_labels_for_the_editor() {
        assert_eq!(ThicknessSolve::<f64>::Fixed.code(), "");
        assert_eq!(
            ThicknessSolve::MarginalRayHeight { height: 0.0 }.code(),
            "M"
        );
        assert_eq!(ThicknessSolve::ChiefRayHeight { height: 0.0 }.code(), "C");
        assert_eq!(
            ThicknessSolve::Pickup {
                from: 0,
                scale: 1.0,
                offset: 0.0
            }
            .code(),
            "P"
        );
        assert!(ThicknessSolve::<f64>::default().is_fixed());
    }
}
