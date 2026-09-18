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

/// Move the image plane to the paraxial focus at `wl`.
///
/// A prescription is normally written with a back focal distance already in it; doing it
/// this way means a sample cannot drift out of focus when a radius is edited.
pub fn focus<S: Scalar>(sys: &mut System<S>, wl: f64) {
    let par = Paraxial::compute(sys, wl);
    let last = sys.image_index();
    sys.surfaces[last - 1].thickness = par.bfd;
}
