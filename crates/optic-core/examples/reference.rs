//! First-order check of published triplet prescriptions.

use optic_core::surface::Profile;
use optic_core::{
    glass_code, lines, Aperture, Dual, Field, Material, Object, Paraxial, Scalar, Surface, System,
};

/// Eleven variables: six radii and five thicknesses, the quantities the book rounds.
type D = Dual<11>;

fn glass(code: &str) -> Material {
    glass_code(code).expect("valid glass code")
}

fn build(rows: &[(f64, f64, &str)], epd: f64, field: f64) -> System<f64> {
    let mut surfaces: Vec<Surface<f64>> = rows
        .iter()
        .map(|(r, t, g)| {
            Surface::new(
                *r,
                *t,
                if g.is_empty() {
                    Material::Vacuum
                } else {
                    glass(g)
                },
            )
        })
        .collect();
    // The book does not print a back focal distance, so the image plane is solved for.
    let last = surfaces.len() - 1;
    surfaces[last] = surfaces[last].clone().autofocus();
    surfaces.push(Surface::plane(0.0, Material::Vacuum));
    // Stop in the wide air space between the flint and the rear crown, the usual place.
    surfaces[3].is_stop = true;

    let mut sys = System::new(Object::Infinity, surfaces)
        .with_aperture(Aperture::EntrancePupilDiameter(epd))
        .with_fields(vec![Field::angle(0.0), Field::angle(field)]);
    optic_core::resolve(&mut sys, lines::D);
    sys
}

fn main() {
    // Figure 12.13: English Patent 155,640 (1919). Published focal length 100 units.
    let a = build(
        &[
            (40.1, 6.0, "613585"),
            (-537.0, 10.0, ""),
            (-47.0, 1.0, "621362"),
            (40.0, 10.8, ""),
            (234.5, 6.0, "613585"),
            (-37.9, 0.0, ""),
        ],
        30.0,
        20.0,
    );

    // Figure 12.14: German Patent 287,089 (1913). Published focal length 100 units.
    let b = build(
        &[
            (16.8, 3.5, "611591"),
            (-116.9, 1.0, ""),
            (-56.3, 0.5, "549458"),
            (15.4, 10.3, ""),
            (f64::INFINITY, 2.1, "611591"),
            (-61.3, 0.0, ""),
        ],
        20.0,
        30.0,
    );

    for (name, sys) in [("Fig 12.13 (EP 155,640)", a), ("Fig 12.14 (DE 287,089)", b)] {
        let par = Paraxial::compute(&sys, lines::D);
        let error = 100.0 * (par.efl - 100.0) / 100.0;
        println!("{name}");
        println!(
            "    EFL {:>9.4}  (published 100)   error {:+.3}%    BFD {:.4}   f/{:.2}",
            par.efl, error, par.bfd, par.fno
        );
        budget(&sys);
        distortion(&sys);
        println!();
    }
}

/// Real chief-ray height against the paraxial prediction, at the published field angles.
///
/// The book plots distortion for both designs on a +/-1% scale, so the magnitude is
/// checkable even though the curve cannot be read off precisely. The stop position is
/// not printed, and distortion is sensitive to it, so this is a sanity range rather
/// than a measurement.
fn distortion(sys: &System<f64>) {
    let par = Paraxial::compute(sys, lines::D);
    print!("    distortion:");
    for field in sys.fields.iter().filter(|f| f.radius() > 0.0) {
        let paraxial = par.image_height(sys, lines::D, *field);
        let real = optic_core::trace(
            sys,
            lines::D,
            optic_core::launch(sys, &par, *field, 0.0, 0.0),
        )
        .image_point()
        .map(|p| (p.x * p.x + p.y * p.y).sqrt());
        match real {
            Some(r) => print!(
                "  {:.0} deg: {:+.3}%",
                field.radius(),
                100.0 * (r - paraxial) / paraxial
            ),
            None => print!("  {:.0} deg: blocked", field.radius()),
        }
    }
    println!();
}

/// How much of the focal-length disagreement the published rounding can account for.
///
/// The prescriptions quote radii and thicknesses to one decimal place, so each carries
/// up to +/-0.05 of unstated value. Differentiating the focal length with respect to all
/// eleven of them at once -- one trace, exact derivatives -- turns that into a number.
///
/// Worst case assumes every rounding conspires; the RMS figure treats them as
/// independent, which is what they are.
fn budget(sys: &System<f64>) {
    let dual = lift(sys);
    let efl = Paraxial::compute(&dual, lines::D).efl;
    let grad = efl.grad();

    const ROUNDING: f64 = 0.05;
    let worst: f64 = grad.iter().map(|g| g.abs() * ROUNDING).sum();
    let rms: f64 = grad
        .iter()
        .map(|g| (g * ROUNDING).powi(2))
        .sum::<f64>()
        .sqrt();

    println!(
        "    focal length uncertainty from the published rounding: +/-{rms:.3} worst case          +/-{worst:.3}"
    );
    let dominant = grad
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.abs().total_cmp(&b.1.abs()))
        .unwrap();
    println!(
        "    most sensitive quantity: {} (d(EFL)/d = {:.1} per unit)",
        if dominant.0 < 6 {
            format!("radius {}", dominant.0 + 1)
        } else {
            format!("thickness {}", dominant.0 - 5)
        },
        dominant.1
    );
}

/// Rebuild the system with every radius and thickness seeded as a variable.
fn lift(sys: &System<f64>) -> System<D> {
    let mut surfaces = Vec::with_capacity(sys.surfaces.len());
    for (i, s) in sys.surfaces.iter().enumerate() {
        let profile = match s.profile {
            Profile::Plane => Profile::Plane,
            Profile::Conic { curvature, conic } => {
                // Seed the radius, not the curvature: the book rounds radii.
                let radius = D::variable(1.0 / curvature, i.min(5));
                Profile::Conic {
                    curvature: D::from_f64(1.0) / radius,
                    conic: D::from_f64(conic),
                }
            }
            _ => unreachable!("reference designs are all-spherical"),
        };
        surfaces.push(Surface {
            profile,
            thickness: if i < 5 {
                D::variable(s.thickness, 6 + i)
            } else {
                D::from_f64(s.thickness)
            },
            thickness_solve: Default::default(),
            material: s.material.clone(),
            semi_diameter: s.semi_diameter.map(D::from_f64),
            is_stop: s.is_stop,
            label: s.label.clone(),
        });
    }
    let mut out = System::new(Object::Infinity, surfaces);
    out.aperture = sys.aperture;
    out.fields = sys.fields.clone();
    out.wavelengths = sys.wavelengths.clone();
    out
}
