//! Verification of `.zmx` reading and writing.
//!
//! The round trip is the backbone: a system written out and read back must be the same
//! system. That catches both directions at once, and it is the property that lets a
//! designer check our results in Zemax and trust that they compared the same lens.

use optic_core::{lines, samples, surface::Profile, system::Object, Material, Paraxial};
use optic_io::zmx;

/// A hand-written prescription for an equiconvex N-BK7 singlet, focused.
/// Its focal length is known in closed form, so this checks the reader against physics
/// rather than against our own writer.
const SINGLET: &str = r#"VERS 190513 0 000000 0
MODE SEQ
NAME Test singlet
UNIT MM X W X CM MR CPMM
ENPD 10
FTYP 0 0 2 1 0 0 0 1
XFLN 0 0
YFLN 0 2
FWGN 1 1
WAVM 1 0.5875618 1
PWAV 1
SURF 0
  TYPE STANDARD
  CURV 0.0
  DISZ INFINITY
SURF 1
  TYPE STANDARD
  CURV 0.01
  DISZ 5.0
  GLAS N-BK7 0 0 1.516800 64.1673 0 0 0 0 0 0
  DIAM 6 0 0 0 1 ""
  STOP
SURF 2
  TYPE STANDARD
  CURV -0.01
  DISZ 95.918036
  DIAM 6 0 0 0 1 ""
SURF 3
  TYPE STANDARD
  CURV 0.0
  DISZ 0.0
"#;

#[test]
fn a_hand_written_prescription_reads_to_the_right_focal_length() {
    let imported = zmx::parse(SINGLET).expect("parsed");
    let sys = &imported.system;

    assert_eq!(sys.title, "Test singlet");
    assert_eq!(
        sys.surfaces.len(),
        3,
        "two refracting surfaces and an image plane"
    );
    assert_eq!(sys.stop_index(), Some(0));
    assert!(matches!(sys.object, Object::Infinity));
    assert_eq!(sys.wavelengths.len(), 1);
    assert_eq!(sys.fields.len(), 2);

    // The thick-lens equation, computed independently in the core test suite.
    let efl = Paraxial::compute(sys, lines::D).efl;
    assert!((efl - 97.580403).abs() < 1e-4, "EFL {efl}");
    assert!(imported.warnings.is_empty(), "{:?}", imported.warnings);
}

#[test]
fn a_system_survives_a_round_trip_unchanged() {
    for original in [samples::cooke_triplet::<f64>(), samples::singlet::<f64>()] {
        let text = zmx::write(&original);
        let back = zmx::parse(&text).expect("re-read what we wrote").system;

        assert_eq!(
            back.surfaces.len(),
            original.surfaces.len(),
            "{}",
            original.title
        );
        assert_eq!(back.fields.len(), original.fields.len());
        assert_eq!(back.wavelengths.len(), original.wavelengths.len());
        assert_eq!(back.stop_index(), original.stop_index());

        for (i, (a, b)) in back
            .surfaces
            .iter()
            .zip(original.surfaces.iter())
            .enumerate()
        {
            assert!(
                (a.thickness - b.thickness).abs() < 1e-9,
                "surface {i} thickness: {} vs {}",
                a.thickness,
                b.thickness
            );
            assert!(
                (a.profile.paraxial_curvature() - b.profile.paraxial_curvature()).abs() < 1e-12,
                "surface {i} curvature"
            );
            assert!(
                (a.material.nd() - b.material.nd()).abs() < 1e-5,
                "surface {i} glass: n_d {} vs {}",
                a.material.nd(),
                b.material.nd()
            );
            assert_eq!(
                a.semi_diameter, b.semi_diameter,
                "surface {i} semi-diameter"
            );
        }

        // The number that matters must survive to full precision.
        let before = Paraxial::compute(&original, lines::D).efl;
        let after = Paraxial::compute(&back, lines::D).efl;
        assert!(
            (before - after).abs() < 1e-6,
            "{}: EFL {before} became {after}",
            original.title
        );
    }
}

#[test]
fn utf16_files_are_decoded() {
    // Older Zemax wrote UTF-16LE with a byte-order mark. Guessing wrong yields a file
    // that looks empty rather than one that fails, so this is worth pinning.
    let mut bytes = vec![0xFF, 0xFE];
    for unit in SINGLET.encode_utf16() {
        bytes.extend_from_slice(&unit.to_le_bytes());
    }
    let text = zmx::decode(&bytes);
    let sys = zmx::parse(&text).expect("parsed from UTF-16").system;
    assert_eq!(sys.title, "Test singlet");
    assert!((Paraxial::compute(&sys, lines::D).efl - 97.580403).abs() < 1e-4);
}

#[test]
fn a_utf8_byte_order_mark_does_not_corrupt_the_first_keyword() {
    let mut bytes = vec![0xEF, 0xBB, 0xBF];
    bytes.extend_from_slice(SINGLET.as_bytes());
    let sys = zmx::parse(&zmx::decode(&bytes)).expect("parsed").system;
    assert_eq!(sys.title, "Test singlet");
}

#[test]
fn an_unknown_glass_falls_back_to_its_index_and_abbe_number() {
    // A file may name a catalogue we do not ship. Zemax records n_d and V_d beside the
    // name, so the design is still usable; the importer must say so rather than
    // silently substituting air.
    let text = SINGLET.replace(
        "GLAS N-BK7 0 0 1.516800 64.1673",
        "GLAS OHARA-S-BSL7 0 0 1.516300 64.1400",
    );
    let imported = zmx::parse(&text).expect("parsed");
    match &imported.system.surfaces[0].material {
        Material::ModelGlass { nd, vd } => {
            assert!((nd - 1.5163).abs() < 1e-6);
            assert!((vd - 64.14).abs() < 1e-6);
        }
        other => panic!("expected model glass, got {other:?}"),
    }
    assert!(
        imported
            .warnings
            .iter()
            .any(|w| w.contains("OHARA-S-BSL7") && w.contains("model glass")),
        "{:?}",
        imported.warnings
    );
}

#[test]
fn an_unknown_glass_with_no_index_is_reported_not_guessed() {
    let text = SINGLET.replace("GLAS N-BK7 0 0 1.516800 64.1673", "GLAS MYSTERY 0 0 0 0");
    let imported = zmx::parse(&text).expect("parsed");
    assert_eq!(imported.system.surfaces[0].material, Material::Vacuum);
    assert!(imported.warnings.iter().any(|w| w.contains("MYSTERY")));
}

#[test]
fn mirrors_are_recognised() {
    let text = SINGLET.replace("GLAS N-BK7 0 0 1.516800 64.1673", "GLAS MIRROR 0 0 0 0");
    let sys = zmx::parse(&text).expect("parsed").system;
    assert_eq!(sys.surfaces[0].material, Material::Mirror);
}

#[test]
fn aspheric_coefficients_map_to_the_right_powers() {
    // Zemax PARM 1 is the r^2 term, PARM 2 the r^4 term, matching our ordering.
    let text = SINGLET.replace(
        "SURF 1\n  TYPE STANDARD\n  CURV 0.01",
        "SURF 1\n  TYPE EVENASPH\n  CURV 0.01\n  PARM 1 1e-5\n  PARM 3 7e-11",
    );
    let sys = zmx::parse(&text).expect("parsed").system;
    match &sys.surfaces[0].profile {
        Profile::EvenAsphere { coeffs, .. } => {
            assert_eq!(coeffs.len(), 3, "gaps must be filled with zeros");
            assert_eq!(coeffs[0], 1e-5);
            assert_eq!(coeffs[1], 0.0);
            assert_eq!(coeffs[2], 7e-11);
        }
        other => panic!("expected an asphere, got {other:?}"),
    }
}

#[test]
fn an_asphere_with_only_zero_coefficients_is_just_a_conic() {
    let text = SINGLET.replace(
        "SURF 1\n  TYPE STANDARD\n  CURV 0.01",
        "SURF 1\n  TYPE EVENASPH\n  CURV 0.01\n  PARM 1 0\n  PARM 2 0",
    );
    let sys = zmx::parse(&text).expect("parsed").system;
    assert!(matches!(sys.surfaces[0].profile, Profile::Conic { .. }));
}

#[test]
fn unsupported_surface_types_import_with_a_warning() {
    let text = SINGLET.replace("SURF 2\n  TYPE STANDARD", "SURF 2\n  TYPE COORDBRK");
    let imported = zmx::parse(&text).expect("parsed anyway");
    assert_eq!(imported.system.surfaces.len(), 3);
    assert!(
        imported.warnings.iter().any(|w| w.contains("COORDBRK")),
        "{:?}",
        imported.warnings
    );
}

#[test]
fn vignetting_factors_are_reported_because_we_do_not_honour_them() {
    let text = SINGLET.replace("PWAV 1", "VDYN 0 0.25\nVCYN 0 0.1\nPWAV 1");
    let imported = zmx::parse(&text).expect("parsed");
    assert!(
        imported.warnings.iter().any(|w| w.contains("vignetting")),
        "{:?}",
        imported.warnings
    );
}

#[test]
fn unused_wavelength_slots_are_discarded() {
    // Zemax writes 24 WAVM lines whether or not they are used. Taking them at face value
    // would run every analysis on two dozen wavelengths, most of them identical.
    let mut text = SINGLET.replace("WAVM 1 0.5875618 1\n", "");
    let mut block = String::new();
    for i in 1..=24 {
        let um = match i {
            1 => 0.4861327,
            2 => 0.5875618,
            3 => 0.6562725,
            _ => 0.55,
        };
        block.push_str(&format!("WAVM {i} {um} 1\n"));
    }
    text = text.replace("PWAV 1", &format!("{block}PWAV 2"));
    text = text.replace("FTYP 0 0 2 1 ", "FTYP 0 0 2 3 ");

    let sys = zmx::parse(&text).expect("parsed").system;
    assert_eq!(sys.wavelengths.len(), 3, "should follow the declared count");
    assert!((sys.wavelengths[1].um - 0.5875618).abs() < 1e-9);
}

#[test]
fn a_finite_object_distance_is_read() {
    let text = SINGLET.replace("  DISZ INFINITY", "  DISZ 250.0");
    let sys = zmx::parse(&text).expect("parsed").system;
    match sys.object {
        Object::Finite { distance } => assert!((distance - 250.0).abs() < 1e-9),
        other => panic!("expected a finite object, got {other:?}"),
    }
}

#[test]
fn field_counts_are_inferred_when_the_file_does_not_declare_them() {
    let text = SINGLET
        .replace("FTYP 0 0 2 1 0 0 0 1", "FTYP 0 0 0 0 0 0 0 1")
        .replace("YFLN 0 2", "YFLN 0 7 14 0 0 0")
        .replace("XFLN 0 0", "XFLN 0 0 0 0 0 0");
    let sys = zmx::parse(&text).expect("parsed").system;
    assert_eq!(
        sys.fields.len(),
        3,
        "trailing zeros are padding, not fields"
    );
    assert_eq!(sys.fields[2].radius(), 14.0);
}

#[test]
fn a_file_with_no_surfaces_is_rejected_clearly() {
    let err = zmx::parse("VERS 190513\nMODE SEQ\nNAME nothing\n").unwrap_err();
    assert!(err.contains("surfaces"), "{err}");
    assert!(zmx::parse("").is_err());
}

#[test]
fn a_missing_stop_is_reported_and_defaulted() {
    let text = SINGLET.replace("  STOP\n", "");
    let imported = zmx::parse(&text).expect("parsed");
    assert_eq!(imported.system.stop_index(), Some(0));
    assert!(imported.warnings.iter().any(|w| w.contains("stop")));
}

#[test]
fn the_declared_primary_wavelength_is_readable() {
    assert_eq!(zmx::primary_index(SINGLET), Some(0));
    assert_eq!(
        zmx::primary_index(&SINGLET.replace("PWAV 1", "PWAV 2")),
        Some(1)
    );
    assert_eq!(zmx::primary_index("NAME nothing"), None);
}

#[test]
fn what_we_write_is_accepted_by_our_own_reader_without_complaint() {
    let text = zmx::write(&samples::cooke_triplet::<f64>());
    let imported = zmx::parse(&text).expect("parsed");
    assert!(
        imported.warnings.is_empty(),
        "our own output produced warnings: {:?}",
        imported.warnings
    );
    assert!(text.starts_with("VERS "));
    assert!(text.contains("MODE SEQ"));
    assert!(text.contains("N-SK16"));
}
