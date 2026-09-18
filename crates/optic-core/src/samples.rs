//! Reference systems used by the test suite and by the examples.
//!
//! These are *working prescriptions*, not citations. Golden numbers derived from them
//! are regression baselines that pin today's behaviour; validating them against
//! published designs is tracked separately in `tests/prescriptions`.

use crate::material::{catalog, lines, Material};
use crate::math::Scalar;
use crate::paraxial::Paraxial;
use crate::system::{Aperture, Field, Object, Surface, System, Wavelength};

/// An equiconvex N-BK7 singlet, 100 mm radii, 5 mm thick, at f/10.
///
/// Its focal length is known in closed form from the thick-lens equation, which is what
/// makes it the first thing the paraxial solver is checked against.
pub fn singlet<S: Scalar>() -> System<S> {
    let mut sys = System::new(
        Object::Infinity,
        vec![
            Surface::new(100.0, 5.0, catalog::N_BK7)
                .stop()
                .labelled("front"),
            Surface::new(-100.0, 97.0, Material::Vacuum)
                .autofocus()
                .labelled("back"),
            Surface::plane(0.0, Material::Vacuum).labelled("image"),
        ],
    )
    .titled("Equiconvex N-BK7 singlet")
    .with_aperture(Aperture::EntrancePupilDiameter(10.0))
    .with_fields(vec![Field::angle(0.0), Field::angle(2.0)]);

    focus(&mut sys, lines::D);
    sys
}

/// A Cooke triplet of nominally 50 mm focal length at f/5.
///
/// The classic three-element anastigmat: crown, flint, crown with the stop between the
/// flint and the rear crown. It exercises everything the kernel has -- several media,
/// an interior stop, meaningful field aberration -- while staying small enough to reason
/// about by hand.
pub fn cooke_triplet<S: Scalar>() -> System<S> {
    let mut sys = System::new(
        Object::Infinity,
        vec![
            Surface::new(22.01359, 3.25861, catalog::N_SK16)
                .with_semi_diameter(9.5)
                .labelled("crown 1 front"),
            Surface::new(-435.76044, 6.00762, Material::Vacuum)
                .with_semi_diameter(9.5)
                .labelled("crown 1 back"),
            Surface::new(-22.21328, 1.0, catalog::F2)
                .with_semi_diameter(4.5)
                .labelled("flint front"),
            Surface::new(20.29192, 2.0, Material::Vacuum)
                .with_semi_diameter(4.5)
                .labelled("flint back"),
            Surface::plane(2.75043, Material::Vacuum)
                .with_semi_diameter(4.0)
                .stop()
                .labelled("stop"),
            Surface::new(79.6836, 2.95208, catalog::N_SK16)
                .with_semi_diameter(6.5)
                .labelled("crown 2 front"),
            Surface::new(-18.39533, 42.0, Material::Vacuum)
                .with_semi_diameter(6.5)
                .autofocus()
                .labelled("crown 2 back"),
            Surface::plane(0.0, Material::Vacuum).labelled("image"),
        ],
    )
    .titled("Cooke triplet, 50 mm f/5")
    .with_aperture(Aperture::EntrancePupilDiameter(10.0))
    .with_fields(vec![
        Field::angle(0.0),
        Field::angle(14.0),
        Field::angle(20.0),
    ])
    .with_wavelengths(vec![
        Wavelength::new(lines::F),
        Wavelength::new(lines::D),
        Wavelength::new(lines::C),
    ]);

    focus(&mut sys, lines::D);
    sys
}

/// A Cooke triplet of moderate aperture and coverage, at f/3.3 over +/-20 degrees.
///
/// From English Patent 155,640 (1919), as tabulated by W. Smith, *Modern Optical
/// Engineering*, figure 12.13. Published focal length 100 units.
///
/// This is a **reference design**, not an invention of ours: its purpose is to be
/// checked against published data. Two caveats come with it. The glasses are given only
/// as six-digit codes, so they are modelled from index and Abbe number rather than a
/// real catalogue entry. And the source prints no stop position, so the stop is placed
/// in the wide air space ahead of the rear crown, which is where a triplet's stop
/// conventionally sits; distortion is sensitive to that choice.
pub fn smith_triplet_moderate<S: Scalar>() -> System<S> {
    reference_triplet(
        "Cooke triplet f/3.3, EP 155,640 (1919)",
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
    )
}

/// A Cooke triplet of small aperture and wide coverage, at f/5 over +/-30 degrees.
///
/// From German Patent 287,089 (1913), as tabulated by W. Smith, *Modern Optical
/// Engineering*, figure 12.14. Published focal length 100 units. The same caveats as
/// [`smith_triplet_moderate`] apply.
pub fn smith_triplet_wide<S: Scalar>() -> System<S> {
    reference_triplet(
        "Cooke triplet f/5 wide field, DE 287,089 (1913)",
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
    )
}

/// Build a published triplet from its prescription rows.
fn reference_triplet<S: Scalar>(
    title: &str,
    rows: &[(f64, f64, &str)],
    epd: f64,
    field: f64,
) -> System<S> {
    let mut surfaces: Vec<Surface<S>> = rows
        .iter()
        .map(|(radius, thickness, code)| {
            let material = if code.is_empty() {
                Material::Vacuum
            } else {
                crate::material::glass_code(code).expect("a valid six-digit glass code")
            };
            Surface::new(*radius, *thickness, material)
        })
        .collect();

    // The sources print no back focal distance, so it is solved rather than invented.
    let last = surfaces.len() - 1;
    surfaces[last] = surfaces[last].clone().autofocus();
    surfaces.push(Surface::plane(0.0, Material::Vacuum).labelled("image"));
    surfaces[3].is_stop = true;

    let mut sys = System::new(Object::Infinity, surfaces)
        .titled(title)
        .with_aperture(Aperture::EntrancePupilDiameter(epd))
        .with_fields(vec![
            Field::angle(0.0),
            Field::angle(field * 0.7),
            Field::angle(field),
        ])
        .with_wavelengths(vec![
            Wavelength::new(lines::F),
            Wavelength::new(lines::D),
            Wavelength::new(lines::C),
        ]);
    crate::solve::resolve(&mut sys, lines::D);
    sys
}

/// Move the image plane to the paraxial focus at `wl`.
///
/// A prescription is normally written with a back focal distance already in it; doing it
/// this way means a sample cannot drift out of focus when a radius is edited.
pub fn focus<S: Scalar>(sys: &mut System<S>, wl: f64) {
    let par = Paraxial::compute(sys, wl);
    let last = sys.image_index();
    sys.surfaces[last - 1].thickness = par.bfd;
}
