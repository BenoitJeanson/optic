//! Reading and writing Zemax `.zmx` prescriptions.
//!
//! The format is line-oriented text: a keyword, then its arguments. Global keywords sit
//! at the left margin; a `SURF n` line opens a surface block whose keywords are indented.
//! That indentation is the only reliable way to tell the two apart, because several
//! keywords are valid in both places.
//!
//! Two decisions shape this importer.
//!
//! **Nothing is dropped silently.** A file that uses a feature this kernel does not model
//! — coordinate breaks, toroids, vignetting factors — imports anyway, and says exactly
//! what it could not honour. A prescription that is quietly wrong is far worse than one
//! that arrives with a list of caveats.
//!
//! **An unknown glass is not a failure.** Zemax records each glass's index and Abbe
//! number next to its name, so a file referring to a catalogue we do not ship still
//! imports, as model glass built from those two numbers. It will be very slightly wrong
//! away from the d line, and it says so.

use optic_core::material::catalog;
use optic_core::surface::Profile;
use optic_core::system::{Aperture, Field, Object, Surface, System, Wavelength};
use optic_core::Material;
use std::collections::BTreeMap;

/// A prescription read from a `.zmx` file, with anything it could not honour.
#[derive(Clone, Debug)]
pub struct Import {
    pub system: System<f64>,
    /// Features present in the file that this kernel does not model.
    pub warnings: Vec<String>,
}

/// Decode a `.zmx` file's bytes.
///
/// Zemax has written both UTF-16LE (older, with a byte-order mark) and UTF-8. Guessing
/// wrong yields a file that looks empty rather than one that fails, so the mark is
/// checked before anything else.
pub fn decode(bytes: &[u8]) -> String {
    if bytes.len() >= 2 && bytes[0] == 0xFF && bytes[1] == 0xFE {
        let units: Vec<u16> = bytes[2..]
            .chunks_exact(2)
            .map(|c| u16::from_le_bytes([c[0], c[1]]))
            .collect();
        return String::from_utf16_lossy(&units);
    }
    // A UTF-8 mark, if present, would otherwise become part of the first keyword.
    let body = bytes.strip_prefix(&[0xEF, 0xBB, 0xBF]).unwrap_or(bytes);
    String::from_utf8_lossy(body).into_owned()
}

#[derive(Default, Debug)]
struct RawSurface {
    kind: String,
    curvature: f64,
    curvature_solve: bool,
    thickness: f64,
    thickness_infinite: bool,
    glass: Option<GlassRef>,
    conic: f64,
    semi_diameter: Option<f64>,
    params: BTreeMap<usize, f64>,
    is_stop: bool,
    comment: String,
}

#[derive(Debug)]
struct GlassRef {
    name: String,
    nd: f64,
    vd: f64,
}

/// Split a line into its keyword and the rest.
fn split_keyword(line: &str) -> Option<(&str, &str)> {
    let trimmed = line.trim_start();
    if trimmed.is_empty() {
        return None;
    }
    match trimmed.find(char::is_whitespace) {
        Some(i) => Some((&trimmed[..i], trimmed[i..].trim())),
        None => Some((trimmed, "")),
    }
}

fn numbers(rest: &str) -> Vec<f64> {
    rest.split_whitespace()
        .filter_map(|t| t.parse::<f64>().ok())
        .collect()
}

fn first_number(rest: &str) -> Option<f64> {
    rest.split_whitespace().next().and_then(|t| t.parse().ok())
}

/// Parse the text of a `.zmx` file.
pub fn parse(text: &str) -> Result<Import, String> {
    let mut warnings = Vec::new();
    let mut title = String::new();
    let mut surfaces: Vec<RawSurface> = Vec::new();
    let mut in_surface = false;

    let mut aperture: Option<Aperture> = None;
    let mut field_type = 0usize;
    let mut xfln: Vec<f64> = Vec::new();
    let mut yfln: Vec<f64> = Vec::new();
    let mut declared_fields: Option<usize> = None;
    let mut declared_waves: Option<usize> = None;
    let mut waves: BTreeMap<usize, (f64, f64)> = BTreeMap::new();
    let mut primary_wave: Option<usize> = None;
    let mut vignetting_used = false;
    let mut unit_warned = false;

    for line in text.lines() {
        let indented = line.starts_with(' ') || line.starts_with('\t');
        let Some((key, rest)) = split_keyword(line) else {
            continue;
        };

        if key == "SURF" {
            surfaces.push(RawSurface::default());
            in_surface = true;
            continue;
        }

        if in_surface && indented {
            let surface = surfaces.last_mut().expect("a surface block is open");
            match key {
                "TYPE" => surface.kind = rest.split_whitespace().next().unwrap_or("").to_string(),
                "CURV" => {
                    let v = numbers(rest);
                    surface.curvature = v.first().copied().unwrap_or(0.0);
                    surface.curvature_solve = v.get(1).copied().unwrap_or(0.0) != 0.0;
                }
                "DISZ" => {
                    if rest.trim().eq_ignore_ascii_case("INFINITY") {
                        surface.thickness_infinite = true;
                    } else {
                        surface.thickness = first_number(rest).unwrap_or(0.0);
                    }
                }
                "GLAS" => {
                    let mut tokens = rest.split_whitespace();
                    let name = tokens.next().unwrap_or("").to_string();
                    let rest_numbers: Vec<f64> =
                        tokens.filter_map(|t| t.parse::<f64>().ok()).collect();
                    // name, code, pickup, then index and Abbe number.
                    surface.glass = Some(GlassRef {
                        name,
                        nd: rest_numbers.get(2).copied().unwrap_or(0.0),
                        vd: rest_numbers.get(3).copied().unwrap_or(0.0),
                    });
                }
                "CONI" => surface.conic = first_number(rest).unwrap_or(0.0),
                "DIAM" => surface.semi_diameter = first_number(rest).filter(|v| *v > 0.0),
                "PARM" => {
                    let v = numbers(rest);
                    if let (Some(i), Some(value)) = (v.first(), v.get(1)) {
                        surface.params.insert(*i as usize, *value);
                    }
                }
                "STOP" => surface.is_stop = true,
                "COMM" => surface.comment = rest.to_string(),
                _ => {}
            }
            continue;
        }

        in_surface = false;
        match key {
            "NAME" => title = rest.to_string(),
            "UNIT" => {
                let unit = rest.split_whitespace().next().unwrap_or("MM");
                if !unit.eq_ignore_ascii_case("MM") && !unit_warned {
                    warnings.push(format!(
                        "file is in {unit}, not millimetres; lengths were not converted"
                    ));
                    unit_warned = true;
                }
            }
            "ENPD" => {
                if let Some(v) = first_number(rest) {
                    aperture = Some(Aperture::EntrancePupilDiameter(v));
                }
            }
            // ENPD wins if the file carries one, so these only apply when nothing has
            // set the aperture yet.
            "FNUM" if aperture.is_none() => {
                if let Some(v) = first_number(rest).filter(|v| *v > 0.0) {
                    aperture = Some(Aperture::ImageSpaceFNumber(v));
                }
            }
            "OBNA" if aperture.is_none() => {
                if let Some(v) = first_number(rest).filter(|v| *v > 0.0) {
                    aperture = Some(Aperture::ObjectSpaceNA(v));
                }
            }
            "FTYP" => {
                let v = numbers(rest);
                field_type = v.first().copied().unwrap_or(0.0) as usize;
                declared_fields = v
                    .get(2)
                    .map(|n| *n as usize)
                    .filter(|n| (1..=12).contains(n));
                declared_waves = v
                    .get(3)
                    .map(|n| *n as usize)
                    .filter(|n| (1..=24).contains(n));
            }
            "XFLN" => xfln = numbers(rest),
            "YFLN" => yfln = numbers(rest),
            "WAVM" => {
                let v = numbers(rest);
                if let (Some(i), Some(um)) = (v.first(), v.get(1)) {
                    waves.insert(*i as usize, (*um, v.get(2).copied().unwrap_or(1.0)));
                }
            }
            "PWAV" => primary_wave = first_number(rest).map(|v| v as usize),
            // All-zero vignetting factors mean the feature is unused, which is the
            // common case; only a non-zero one changes what the prescription means.
            "VDXN" | "VDYN" | "VCXN" | "VCYN" | "VANN"
                if numbers(rest).iter().any(|v| *v != 0.0) =>
            {
                vignetting_used = true;
            }
            _ => {}
        }
    }

    if surfaces.len() < 2 {
        return Err("no surfaces found: this does not look like a .zmx prescription".into());
    }

    // Surface 0 is the object; the last is the image plane.
    let object_raw = surfaces.remove(0);
    let object = if object_raw.thickness_infinite {
        Object::Infinity
    } else {
        Object::Finite {
            distance: object_raw.thickness.abs(),
        }
    };

    let mut built = Vec::with_capacity(surfaces.len());
    for (i, raw) in surfaces.iter().enumerate() {
        let number = i + 1;
        let kind = raw.kind.to_ascii_uppercase();
        let profile = match kind.as_str() {
            "" | "STANDARD" => conic_profile(raw.curvature, raw.conic),
            "EVENASPH" => {
                let highest = raw.params.keys().copied().max().unwrap_or(0);
                let coeffs: Vec<f64> = (1..=highest)
                    .map(|k| raw.params.get(&k).copied().unwrap_or(0.0))
                    .collect();
                if coeffs.iter().all(|c| *c == 0.0) {
                    conic_profile(raw.curvature, raw.conic)
                } else {
                    Profile::EvenAsphere {
                        curvature: raw.curvature,
                        conic: raw.conic,
                        coeffs,
                    }
                }
            }
            other => {
                warnings.push(format!(
                    "surface {number}: type {other} is not modelled yet; treated as a \
                     {} surface",
                    if raw.curvature == 0.0 {
                        "flat"
                    } else {
                        "conic"
                    }
                ));
                conic_profile(raw.curvature, raw.conic)
            }
        };

        if raw.curvature_solve {
            warnings.push(format!(
                "surface {number}: a curvature solve was dropped; the radius is fixed at \
                 its current value"
            ));
        }

        let material = match &raw.glass {
            None => Material::Vacuum,
            Some(g) if g.name.eq_ignore_ascii_case("MIRROR") => Material::Mirror,
            Some(g) => {
                match catalog::by_name(&g.name).or_else(|| optic_core::glass_code(&g.name)) {
                    Some(m) => m,
                    None if g.nd > 1.0 && g.vd > 0.0 => {
                        warnings.push(format!(
                            "surface {number}: {} is not in the built-in catalogue; using a \
                         model glass from its n_d {:.5} and V_d {:.2}",
                            g.name, g.nd, g.vd
                        ));
                        Material::ModelGlass { nd: g.nd, vd: g.vd }
                    }
                    None => {
                        warnings.push(format!(
                            "surface {number}: glass {} is unknown and the file carries no \
                         index for it; treated as air",
                            g.name
                        ));
                        Material::Vacuum
                    }
                }
            }
        };

        built.push(Surface {
            profile,
            thickness: if raw.thickness_infinite {
                0.0
            } else {
                raw.thickness
            },
            thickness_solve: Default::default(),
            material,
            semi_diameter: raw.semi_diameter,
            is_stop: raw.is_stop,
            label: raw.comment.clone(),
        });
    }

    // Fields.
    let count = declared_fields.unwrap_or_else(|| infer_count(&yfln, &xfln));
    let fields: Vec<Field> = (0..count.max(1))
        .map(|i| {
            let x = xfln.get(i).copied().unwrap_or(0.0);
            let y = yfln.get(i).copied().unwrap_or(0.0);
            match field_type {
                0 => Field::Angle { x, y },
                _ => Field::Height { x, y },
            }
        })
        .collect();
    if field_type > 1 {
        warnings.push(
            "fields are given as image heights, which this kernel does not solve for; \
             they were read as object heights"
                .into(),
        );
    }

    // Wavelengths. Zemax writes 24 slots whether or not they are used, so the declared
    // count matters: without it the analysis would run on two dozen identical lines.
    let wave_count = declared_waves.unwrap_or_else(|| waves.len().clamp(1, 3));
    let mut wavelengths: Vec<Wavelength> = (1..=wave_count)
        .filter_map(|i| waves.get(&i))
        .map(|(um, weight)| Wavelength {
            um: *um,
            weight: *weight,
        })
        .collect();
    if wavelengths.is_empty() {
        wavelengths.push(Wavelength::new(optic_core::lines::D));
        warnings.push("no wavelengths found; assuming the d line".into());
    }

    if vignetting_used {
        warnings.push(
            "this prescription uses vignetting factors, which are not implemented yet; \
             rays are clipped on the clear apertures instead, so the outer field will \
             not match Zemax"
                .into(),
        );
    }

    let mut system = System::new(object, built);
    system.title = if title.is_empty() {
        "Imported from .zmx".into()
    } else {
        title
    };
    system.fields = fields;
    system.wavelengths = wavelengths;
    system.aperture = aperture.unwrap_or_else(|| {
        warnings.push("no aperture found in the file; assuming a 10 mm entrance pupil".into());
        Aperture::EntrancePupilDiameter(10.0)
    });

    if system.stop_index().is_none() {
        warnings.push("no stop surface was marked; using the first surface".into());
        system.surfaces[0].is_stop = true;
    }

    if let Some(p) = primary_wave {
        if p == 0 || p > system.wavelengths.len() {
            warnings.push(format!(
                "primary wavelength {p} is out of range; using the middle one"
            ));
        }
    }

    Ok(Import { system, warnings })
}

fn conic_profile(curvature: f64, conic: f64) -> Profile<f64> {
    if curvature == 0.0 {
        Profile::Plane
    } else {
        Profile::Conic { curvature, conic }
    }
}

/// How many fields a file defines, when it does not say.
///
/// The first entry always counts, even when it is the axis; beyond that, trailing zeros
/// are padding rather than field points.
fn infer_count(y: &[f64], x: &[f64]) -> usize {
    let last_used = (0..y.len().max(x.len())).rev().find(|i| {
        y.get(*i).copied().unwrap_or(0.0) != 0.0 || x.get(*i).copied().unwrap_or(0.0) != 0.0
    });
    last_used.map(|i| i + 1).unwrap_or(1)
}

/// Index of the primary wavelength declared by the file, if any (zero-based).
pub fn primary_index(text: &str) -> Option<usize> {
    for line in text.lines() {
        if let Some(("PWAV", rest)) = split_keyword(line) {
            if !(line.starts_with(' ') || line.starts_with('\t')) {
                return first_number(rest).map(|v| (v as usize).saturating_sub(1));
            }
        }
    }
    None
}

/// Write a system as a `.zmx` prescription.
///
/// The point of writing this format is not archival — we have our own — but verification.
/// A designer can open the result in Zemax and check our numbers against a tool they
/// already trust, which is the only argument for correctness that carries weight.
pub fn write(sys: &System<f64>) -> String {
    let mut out = String::new();
    let line = |out: &mut String, s: &str| {
        out.push_str(s);
        out.push('\n');
    };

    line(&mut out, "VERS 190513 0 000000 0");
    line(&mut out, "MODE SEQ");
    line(&mut out, &format!("NAME {}", sys.title));
    line(&mut out, "UNIT MM X W X CM MR CPMM");
    line(&mut out, "GCAT SCHOTT");

    match sys.aperture {
        Aperture::EntrancePupilDiameter(d) => line(&mut out, &format!("ENPD {d}")),
        Aperture::ImageSpaceFNumber(f) => line(&mut out, &format!("FNUM {f} 0")),
        Aperture::ObjectSpaceNA(na) => line(&mut out, &format!("OBNA {na} 0")),
        Aperture::StopDiameter(d) => {
            // No direct equivalent; the stop's own DIAM carries it below.
            line(&mut out, &format!("ENPD {d}"));
        }
    }

    let angular = matches!(sys.fields.first(), Some(Field::Angle { .. }));
    line(
        &mut out,
        &format!(
            "FTYP {} 0 {} {} 0 0 0 1",
            if angular { 0 } else { 1 },
            sys.fields.len(),
            sys.wavelengths.len()
        ),
    );
    let (xs, ys): (Vec<f64>, Vec<f64>) = sys
        .fields
        .iter()
        .map(|f| match *f {
            Field::Angle { x, y } | Field::Height { x, y } => (x, y),
        })
        .unzip();
    line(&mut out, &format!("XFLN {}", join(&xs)));
    line(&mut out, &format!("YFLN {}", join(&ys)));
    line(
        &mut out,
        &format!("FWGN {}", join(&vec![1.0; sys.fields.len().max(1)])),
    );
    for (i, w) in sys.wavelengths.iter().enumerate() {
        line(&mut out, &format!("WAVM {} {} {}", i + 1, w.um, w.weight));
    }
    line(&mut out, &format!("PWAV {}", sys.wavelengths.len() / 2 + 1));

    // Surface 0 is the object.
    line(&mut out, "SURF 0");
    line(&mut out, "  TYPE STANDARD");
    line(&mut out, "  CURV 0.0");
    match sys.object {
        Object::Infinity => line(&mut out, "  DISZ INFINITY"),
        Object::Finite { distance } => line(&mut out, &format!("  DISZ {distance}")),
    }

    for (i, s) in sys.surfaces.iter().enumerate() {
        line(&mut out, &format!("SURF {}", i + 1));
        match &s.profile {
            Profile::EvenAsphere {
                curvature,
                conic,
                coeffs,
            } => {
                line(&mut out, "  TYPE EVENASPH");
                line(&mut out, &format!("  CURV {curvature}"));
                if *conic != 0.0 {
                    line(&mut out, &format!("  CONI {conic}"));
                }
                for (k, c) in coeffs.iter().enumerate() {
                    line(&mut out, &format!("  PARM {} {c}", k + 1));
                }
            }
            Profile::Conic { curvature, conic } => {
                line(&mut out, "  TYPE STANDARD");
                line(&mut out, &format!("  CURV {curvature}"));
                if *conic != 0.0 {
                    line(&mut out, &format!("  CONI {conic}"));
                }
            }
            Profile::Plane => {
                line(&mut out, "  TYPE STANDARD");
                line(&mut out, "  CURV 0.0");
            }
        }
        // The image plane's thickness is meaningless; Zemax expects 0.
        let thickness = if i + 1 == sys.surfaces.len() {
            0.0
        } else {
            s.thickness
        };
        line(&mut out, &format!("  DISZ {thickness}"));

        if let Some(g) = glass_record(&s.material) {
            line(&mut out, &format!("  GLAS {g}"));
        }
        if let Some(sd) = s.semi_diameter {
            line(&mut out, &format!("  DIAM {sd} 0 0 0 1 \"\""));
        }
        if s.is_stop {
            line(&mut out, "  STOP");
        }
        if !s.label.is_empty() {
            line(&mut out, &format!("  COMM {}", s.label));
        }
    }
    out
}

fn join(values: &[f64]) -> String {
    values
        .iter()
        .map(|v| v.to_string())
        .collect::<Vec<_>>()
        .join(" ")
}

/// The `GLAS` argument list for a medium, or `None` for air.
fn glass_record(m: &Material) -> Option<String> {
    match m {
        Material::Vacuum | Material::Air => None,
        Material::Mirror => Some("MIRROR 0 0 0 0 0 0".into()),
        other => {
            let name = catalog::PUBLISHED
                .iter()
                .find(|(_, g, ..)| *g == other)
                .map(|(n, ..)| (*n).to_string())
                .unwrap_or_else(|| format!("MODEL_{:.4}_{:.2}", other.nd(), other.abbe()));
            Some(format!(
                "{name} 0 0 {:.6} {:.4} 0 0 0 0 0 0",
                other.nd(),
                other.abbe()
            ))
        }
    }
}
