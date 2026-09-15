//! The JSON analysis API.
//!
//! This is deliberately ordinary Rust with no WebAssembly in sight, so it can be tested
//! natively rather than only through a browser. The C ABI in `lib.rs` is a thin shim
//! over [`analyze_json`].
//!
//! The request describes a system the way a prescription table does — one row per
//! surface, the last row being the image plane — and the response carries everything the
//! page needs to draw: first-order data, closed polygons for the glass, ray polylines
//! and spot diagrams.

use optic_core::{
    catalog, launch, samples,
    surface::Profile,
    system::{Aperture, Field, Object, Surface, Wavelength},
    trace, Material, Paraxial, System,
};
use serde::{Deserialize, Serialize};

/// One row of the prescription, as the editor presents it.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct SurfaceSpec {
    /// Radius of curvature. Use `null` or a non-finite value for a flat surface.
    #[serde(default)]
    pub radius: Option<f64>,
    pub thickness: f64,
    /// Glass name, or an empty string / "AIR" for no glass, or "MIRROR".
    #[serde(default)]
    pub glass: String,
    #[serde(default)]
    pub semi_diameter: Option<f64>,
    #[serde(default)]
    pub conic: f64,
    #[serde(default)]
    pub stop: bool,
    #[serde(default)]
    pub label: String,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct SystemSpec {
    #[serde(default)]
    pub title: String,
    /// Entrance pupil diameter, in millimetres.
    pub entrance_pupil_diameter: f64,
    /// Wavelengths in micrometres. The middle one is treated as primary.
    pub wavelengths: Vec<f64>,
    /// Field angles in degrees, for an object at infinity.
    pub fields: Vec<f64>,
    /// Surfaces in order; the last one is the image plane.
    pub surfaces: Vec<SurfaceSpec>,
}

#[derive(Clone, Debug, Deserialize)]
pub struct Request {
    pub system: SystemSpec,
    /// Rays per field in the layout drawing.
    #[serde(default = "default_fan")]
    pub rays_per_fan: usize,
    /// Grid resolution across the pupil for spot diagrams.
    #[serde(default = "default_grid")]
    pub spot_grid: usize,
    /// Move the image plane to the paraxial focus before analysing.
    #[serde(default)]
    pub refocus: bool,
}

fn default_fan() -> usize {
    11
}
fn default_grid() -> usize {
    25
}

#[derive(Clone, Debug, Serialize)]
pub struct FirstOrder {
    pub efl: f64,
    pub bfd: f64,
    pub epd: f64,
    pub fno: f64,
    pub entrance_pupil_z: f64,
    pub total_track: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct Polyline {
    pub points: Vec<[f64; 2]>,
}

#[derive(Clone, Debug, Serialize)]
pub struct RayPath {
    pub field: f64,
    pub points: Vec<[f64; 2]>,
    /// False when the ray was stopped before the image plane.
    pub complete: bool,
}

#[derive(Clone, Debug, Serialize)]
pub struct Layout {
    /// Closed outlines of the glass elements.
    pub elements: Vec<Polyline>,
    /// Every surface profile, including the ones in air.
    pub surfaces: Vec<Polyline>,
    pub rays: Vec<RayPath>,
    pub image_z: f64,
    pub bounds: [f64; 4],
}

#[derive(Clone, Debug, Serialize)]
pub struct SpotPoint {
    pub x: f64,
    pub y: f64,
    /// Index into the wavelength list, for colouring.
    pub w: usize,
}

#[derive(Clone, Debug, Serialize)]
pub struct Spot {
    pub field: f64,
    /// RMS spot radius about the centroid, micrometres.
    pub rms: f64,
    /// Radius enclosing every ray, micrometres.
    pub geometric: f64,
    /// Airy disc radius for the working f-number, micrometres.
    pub airy: f64,
    /// Paraxial image height, millimetres.
    pub image_height: f64,
    pub points: Vec<SpotPoint>,
    /// Fraction of launched rays that reached the image plane.
    pub throughput: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct Analysis {
    pub ok: bool,
    pub title: String,
    pub first_order: FirstOrder,
    pub layout: Layout,
    pub spots: Vec<Spot>,
    pub warnings: Vec<String>,
    /// The prescription actually analysed, which differs from the request if refocused.
    pub system: SystemSpec,
}

/// Resolve a glass name to a medium.
pub fn material_for(name: &str) -> Result<Material, String> {
    let key = name.trim();
    if key.is_empty() || key.eq_ignore_ascii_case("air") || key.eq_ignore_ascii_case("vacuum") {
        return Ok(Material::Vacuum);
    }
    if key.eq_ignore_ascii_case("mirror") {
        return Ok(Material::Mirror);
    }
    catalog::by_name(key).ok_or_else(|| format!("unknown glass \"{key}\""))
}

/// Every glass the demo offers, including the two pseudo-media.
pub fn glass_names() -> Vec<String> {
    let mut v = vec!["AIR".to_string(), "MIRROR".to_string()];
    v.extend(catalog::PUBLISHED.iter().map(|(n, ..)| n.to_string()));
    v
}

impl SystemSpec {
    /// Turn a prescription into a traceable system.
    pub fn build(&self) -> Result<System<f64>, String> {
        if self.surfaces.len() < 2 {
            return Err("a system needs at least one surface and an image plane".into());
        }
        if self.wavelengths.is_empty() {
            return Err("at least one wavelength is required".into());
        }
        if self.entrance_pupil_diameter.is_nan() || self.entrance_pupil_diameter <= 0.0 {
            return Err("the entrance pupil diameter must be positive".into());
        }

        let mut surfaces = Vec::with_capacity(self.surfaces.len());
        for (i, s) in self.surfaces.iter().enumerate() {
            let radius = s.radius.unwrap_or(f64::INFINITY);
            let profile = if s.conic == 0.0 {
                Profile::sphere(radius)
            } else if radius.is_finite() && radius != 0.0 {
                Profile::Conic {
                    curvature: 1.0 / radius,
                    conic: s.conic,
                }
            } else {
                Profile::Plane
            };
            let mut surf = Surface {
                profile,
                thickness: s.thickness,
                material: material_for(&s.glass).map_err(|e| format!("surface {}: {e}", i + 1))?,
                semi_diameter: s.semi_diameter.filter(|v| *v > 0.0),
                is_stop: s.stop,
                label: s.label.clone(),
            };
            // The image plane never carries glass or a thickness of its own.
            if i + 1 == self.surfaces.len() {
                surf.material = Material::Vacuum;
                surf.thickness = 0.0;
            }
            surfaces.push(surf);
        }

        let mut sys = System::new(Object::Infinity, surfaces)
            .titled(&self.title)
            .with_aperture(Aperture::EntrancePupilDiameter(
                self.entrance_pupil_diameter,
            ))
            .with_fields(self.fields.iter().map(|d| Field::angle(*d)).collect())
            .with_wavelengths(
                self.wavelengths
                    .iter()
                    .map(|w| Wavelength::new(*w))
                    .collect(),
            );

        if sys.stop_index().is_none() {
            // Without a stop the pupil is undefined; the first surface is the usual
            // fallback and matches what a prescription with no STO row means.
            sys.surfaces[0].is_stop = true;
        }
        if sys.fields.is_empty() {
            sys.fields = vec![Field::angle(0.0)];
        }
        Ok(sys)
    }
}

/// The wavelength treated as primary: the middle of the list.
fn primary(wavelengths: &[f64]) -> f64 {
    wavelengths[wavelengths.len() / 2]
}

/// Largest radius at which a profile still has a defined sag.
fn usable_semi_diameter(profile: &Profile<f64>, wanted: f64) -> f64 {
    let mut r = wanted;
    for _ in 0..64 {
        if r <= 0.0 || profile.sag(r * r).is_ok() {
            return r;
        }
        r *= 0.97;
    }
    0.0
}

fn sample_profile(profile: &Profile<f64>, vertex_z: f64, sd: f64, n: usize) -> Vec<[f64; 2]> {
    let sd = usable_semi_diameter(profile, sd);
    (0..=n)
        .map(|i| {
            let y = -sd + 2.0 * sd * i as f64 / n as f64;
            let z = vertex_z + profile.sag(y * y).unwrap_or(0.0);
            [z, y]
        })
        .collect()
}

/// Run the full analysis.
pub fn analyze(request: &Request) -> Result<Analysis, String> {
    let mut sys = request.system.build()?;
    let wl = primary(&request.system.wavelengths);
    let mut warnings = Vec::new();

    if request.refocus {
        samples::focus(&mut sys, wl);
    }

    let par = Paraxial::compute(&sys, wl);
    if !par.efl.is_finite() {
        return Err("the system has no focal length: check radii and glasses".into());
    }

    let vertices = sys.vertices();
    let image_z = vertices[sys.image_index()];
    let fields: Vec<f64> = sys.fields.iter().map(|f| f.radius()).collect();

    // Semi-diameters that were not given are taken from where the rays actually land.
    let mut traced_sd = vec![0.0f64; sys.surfaces.len()];
    let fan = request.rays_per_fan.clamp(3, 101);
    for field in &sys.fields {
        for i in 0..fan {
            let py = -1.0 + 2.0 * i as f64 / (fan - 1) as f64;
            let t = trace(&sys, wl, launch(&sys, &par, *field, 0.0, py));
            for h in &t.hits {
                traced_sd[h.surface] = traced_sd[h.surface].max(h.local.radius());
            }
        }
    }
    let semi: Vec<f64> = sys
        .surfaces
        .iter()
        .enumerate()
        .map(|(i, s)| match s.semi_diameter {
            Some(v) => v,
            None => (traced_sd[i] * 1.05).max(1e-3),
        })
        .collect();

    // Glass elements: a closed outline spanning each gap that contains glass. Both faces
    // are drawn to the same height so the outline closes into a lens shape.
    let mut elements = Vec::new();
    for i in 0..sys.surfaces.len() - 1 {
        let m = sys.medium_after(i);
        if matches!(m, Material::Vacuum | Material::Air | Material::Mirror) {
            continue;
        }
        let edge = semi[i].max(semi[i + 1]);
        let mut pts = sample_profile(&sys.surfaces[i].profile, vertices[i], edge, 48);
        let mut back = sample_profile(&sys.surfaces[i + 1].profile, vertices[i + 1], edge, 48);
        back.reverse();
        pts.extend(back);
        elements.push(Polyline { points: pts });
    }

    let surfaces: Vec<Polyline> = sys
        .surfaces
        .iter()
        .enumerate()
        .map(|(i, s)| Polyline {
            points: sample_profile(&s.profile, vertices[i], semi[i], 48),
        })
        .collect();

    // Ray paths, started from a plane comfortably in front of the first surface.
    let track = image_z.max(1.0);
    let z_start = -0.15 * track;
    let mut rays = Vec::new();
    for (fi, field) in sys.fields.iter().enumerate() {
        for i in 0..fan {
            let py = -1.0 + 2.0 * i as f64 / (fan - 1) as f64;
            let start = launch(&sys, &par, *field, 0.0, py);
            // The launch plane is the entrance pupil, which may sit inside the lens, so
            // walk the ray back to a plane in front of everything before drawing it.
            let back = (z_start - start.pos.z) / start.dir.z;
            let head = start.pos + start.dir * back;

            let traced = trace(&sys, wl, start);
            let mut points = vec![[head.z, head.y]];
            points.extend(traced.hits.iter().map(|h| [h.global.z, h.global.y]));
            rays.push(RayPath {
                field: fields[fi],
                points,
                complete: traced.is_complete(),
            });
        }
    }

    // Spot diagrams, over the full wavelength set.
    let grid = request.spot_grid.clamp(3, 81);
    let mut spots = Vec::new();
    for (fi, field) in sys.fields.iter().enumerate() {
        let mut pts: Vec<(f64, f64, usize)> = Vec::new();
        let mut launched = 0usize;
        for (wi, w) in request.system.wavelengths.iter().enumerate() {
            let par_w = Paraxial::compute(&sys, *w);
            for a in 0..grid {
                for b in 0..grid {
                    let px = -1.0 + 2.0 * (a as f64 + 0.5) / grid as f64;
                    let py = -1.0 + 2.0 * (b as f64 + 0.5) / grid as f64;
                    if px * px + py * py > 1.0 {
                        continue;
                    }
                    launched += 1;
                    if let Some(p) =
                        trace(&sys, *w, launch(&sys, &par_w, *field, px, py)).image_point()
                    {
                        pts.push((p.x, p.y, wi));
                    }
                }
            }
        }

        if pts.is_empty() {
            warnings.push(format!(
                "no ray reached the image plane at {:.1} degrees",
                fields[fi]
            ));
            spots.push(Spot {
                field: fields[fi],
                rms: f64::NAN,
                geometric: f64::NAN,
                airy: 1.22 * wl * par.fno,
                image_height: par.image_height(&sys, wl, *field),
                points: Vec::new(),
                throughput: 0.0,
            });
            continue;
        }

        let k = pts.len() as f64;
        let cx = pts.iter().map(|p| p.0).sum::<f64>() / k;
        let cy = pts.iter().map(|p| p.1).sum::<f64>() / k;
        let mut sum2 = 0.0;
        let mut worst: f64 = 0.0;
        let points: Vec<SpotPoint> = pts
            .iter()
            .map(|(x, y, w)| {
                let (dx, dy) = ((x - cx) * 1000.0, (y - cy) * 1000.0);
                sum2 += dx * dx + dy * dy;
                worst = worst.max((dx * dx + dy * dy).sqrt());
                SpotPoint {
                    x: dx,
                    y: dy,
                    w: *w,
                }
            })
            .collect();

        spots.push(Spot {
            field: fields[fi],
            rms: (sum2 / k).sqrt(),
            geometric: worst,
            airy: 1.22 * wl * par.fno,
            image_height: par.image_height(&sys, wl, *field),
            points,
            throughput: pts.len() as f64 / launched.max(1) as f64,
        });
    }

    for s in &spots {
        if s.throughput > 0.0 && s.throughput < 0.98 {
            warnings.push(format!(
                "{:.0}% of rays are vignetted at {:.1} degrees",
                100.0 * (1.0 - s.throughput),
                s.field
            ));
        }
    }

    let y_max = semi
        .iter()
        .cloned()
        .fold(0.0f64, f64::max)
        .max(
            spots
                .iter()
                .map(|s| s.image_height.abs())
                .fold(0.0, f64::max),
        )
        .max(1e-3);

    // Report back the prescription actually analysed, so a refocus is visible in the editor.
    let mut echoed = request.system.clone();
    for (spec, surf) in echoed.surfaces.iter_mut().zip(sys.surfaces.iter()) {
        spec.thickness = surf.thickness;
    }

    Ok(Analysis {
        ok: true,
        title: sys.title.clone(),
        first_order: FirstOrder {
            efl: par.efl,
            bfd: par.bfd,
            epd: par.epd,
            fno: par.fno,
            entrance_pupil_z: par.ep_z,
            total_track: image_z,
        },
        layout: Layout {
            elements,
            surfaces,
            rays,
            image_z,
            bounds: [z_start, -y_max * 1.15, image_z, y_max * 1.15],
        },
        spots,
        warnings,
        system: echoed,
    })
}

/// Parse, analyse and serialise. Errors come back as `{"ok": false, "error": ...}`.
pub fn analyze_json(request: &str) -> String {
    let parsed: Result<Request, _> = serde_json::from_str(request);
    let result = match parsed {
        Err(e) => Err(format!("could not read the request: {e}")),
        Ok(req) => analyze(&req),
    };
    match result {
        Ok(a) => serde_json::to_string(&a).unwrap_or_else(|e| error_json(&e.to_string())),
        Err(e) => error_json(&e),
    }
}

fn error_json(message: &str) -> String {
    serde_json::json!({ "ok": false, "error": message }).to_string()
}

/// The built-in systems offered by the demo, as prescriptions.
pub fn presets_json() -> String {
    serde_json::json!({
        "glasses": glass_names(),
        "presets": [spec_of(&samples::cooke_triplet::<f64>()), spec_of(&samples::singlet::<f64>())],
    })
    .to_string()
}

/// Describe a built system as the prescription the editor edits.
fn spec_of(sys: &System<f64>) -> SystemSpec {
    let surfaces = sys
        .surfaces
        .iter()
        .map(|s| {
            let (radius, conic) = match &s.profile {
                Profile::Plane => (None, 0.0),
                Profile::Conic { curvature, conic } => (curvature_to_radius(*curvature), *conic),
                Profile::EvenAsphere {
                    curvature, conic, ..
                } => (curvature_to_radius(*curvature), *conic),
            };
            SurfaceSpec {
                radius,
                thickness: s.thickness,
                glass: glass_name_of(&s.material),
                semi_diameter: s.semi_diameter,
                conic,
                stop: s.is_stop,
                label: s.label.clone(),
            }
        })
        .collect();

    SystemSpec {
        title: sys.title.clone(),
        entrance_pupil_diameter: match sys.aperture {
            Aperture::EntrancePupilDiameter(d) => d,
            _ => 10.0,
        },
        wavelengths: sys.wavelengths.iter().map(|w| w.um).collect(),
        fields: sys.fields.iter().map(|f| f.radius()).collect(),
        surfaces,
    }
}

fn curvature_to_radius(c: f64) -> Option<f64> {
    if c == 0.0 {
        None
    } else {
        Some(1.0 / c)
    }
}

fn glass_name_of(m: &Material) -> String {
    match m {
        Material::Vacuum | Material::Air => "AIR".into(),
        Material::Mirror => "MIRROR".into(),
        other => catalog::PUBLISHED
            .iter()
            .find(|(_, g, ..)| *g == other)
            .map(|(n, ..)| n.to_string())
            .unwrap_or_else(|| "AIR".into()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn presets() -> serde_json::Value {
        serde_json::from_str(&presets_json()).unwrap()
    }

    fn cooke() -> SystemSpec {
        serde_json::from_value(presets()["presets"][0].clone()).unwrap()
    }

    fn run(spec: SystemSpec) -> Analysis {
        analyze(&Request {
            system: spec,
            rays_per_fan: 9,
            spot_grid: 15,
            refocus: false,
        })
        .expect("analysis succeeded")
    }

    #[test]
    fn a_preset_survives_the_round_trip_through_the_editor() {
        // The editor reads a prescription, may write it back unchanged, and must get the
        // same system. Anything lost here would silently alter a design on reload.
        let spec = cooke();
        let rebuilt = spec_of(&spec.build().unwrap());

        assert_eq!(rebuilt.surfaces.len(), spec.surfaces.len());
        assert_eq!(
            rebuilt.entrance_pupil_diameter,
            spec.entrance_pupil_diameter
        );
        assert_eq!(rebuilt.wavelengths, spec.wavelengths);
        assert_eq!(rebuilt.fields, spec.fields);
        for (a, b) in rebuilt.surfaces.iter().zip(spec.surfaces.iter()) {
            assert_eq!(a.glass, b.glass);
            assert_eq!(a.thickness, b.thickness);
            assert_eq!(a.stop, b.stop);
            assert_eq!(a.conic, b.conic);
            match (a.radius, b.radius) {
                (None, None) => {}
                (Some(x), Some(y)) => assert!((x - y).abs() < 1e-9, "{x} vs {y}"),
                _ => panic!("radius changed shape: {:?} vs {:?}", a.radius, b.radius),
            }
        }
    }

    #[test]
    fn the_cooke_preset_analyses_to_its_known_first_order_data() {
        let a = run(cooke());
        assert!(a.ok);
        assert!((a.first_order.efl - 50.021332).abs() < 1e-4);
        assert!((a.first_order.fno - 5.002133).abs() < 1e-4);
        assert!((a.first_order.epd - 10.0).abs() < 1e-9);
        assert!(a.first_order.total_track > 0.0);
    }

    #[test]
    fn the_layout_has_one_outline_per_glass_element() {
        let a = run(cooke());
        assert_eq!(a.layout.elements.len(), 3, "a triplet has three elements");
        assert_eq!(a.layout.surfaces.len(), cooke().surfaces.len());
        for e in &a.layout.elements {
            assert!(e.points.len() > 8);
            assert!(e
                .points
                .iter()
                .all(|p| p[0].is_finite() && p[1].is_finite()));
        }
    }

    #[test]
    fn rays_are_drawn_from_in_front_of_the_lens_to_wherever_they_stop() {
        // Some rays really are blocked: the preset's clear apertures vignette the outer
        // field, which is exactly what the drawing should show. What must hold is that
        // every path starts in front of the system, advances monotonically, and either
        // reaches the image plane or stops at the surface that blocked it.
        let a = run(cooke());
        let spec = cooke();
        assert!(!a.layout.rays.is_empty());

        for r in &a.layout.rays {
            assert!(r.points[0][0] < 0.0, "the ray does not start in front");
            assert!(
                r.points.windows(2).all(|w| w[1][0] >= w[0][0] - 1e-9),
                "a ray path doubles back"
            );
            if r.complete {
                assert_eq!(r.points.len(), spec.surfaces.len() + 1);
                let last = r.points.last().unwrap();
                assert!((last[0] - a.layout.image_z).abs() < 1e-9);
            } else {
                assert!(r.points.len() < spec.surfaces.len() + 1);
            }
        }

        // On axis nothing should be clipped at all.
        assert!(
            a.layout
                .rays
                .iter()
                .filter(|r| r.field == 0.0)
                .all(|r| r.complete),
            "an axial ray was blocked"
        );
        let complete = a.layout.rays.iter().filter(|r| r.complete).count();
        assert!(
            complete * 5 > a.layout.rays.len() * 4,
            "only {complete} of {} rays got through",
            a.layout.rays.len()
        );
    }

    #[test]
    fn spot_diagrams_are_reported_per_field_with_plausible_sizes() {
        let a = run(cooke());
        assert_eq!(a.spots.len(), cooke().fields.len());
        for s in &a.spots {
            assert!(!s.points.is_empty());
            assert!(s.rms > 1.0 && s.rms < 100.0, "RMS {} um", s.rms);
            assert!(s.geometric >= s.rms, "geometric radius below RMS");
            assert!(s.airy > 0.1 && s.airy < 20.0, "Airy radius {} um", s.airy);
            assert!(s.throughput > 0.9);
            // Points are centroid-relative, so they must straddle zero.
            assert!(s.points.iter().any(|p| p.x <= 0.0));
            assert!(s.points.iter().any(|p| p.x >= 0.0));
        }
        // Image height grows with field angle.
        let heights: Vec<f64> = a.spots.iter().map(|s| s.image_height).collect();
        assert!(heights.windows(2).all(|w| w[1] > w[0]));
    }

    #[test]
    fn spot_points_carry_their_wavelength() {
        let a = run(cooke());
        let used: std::collections::BTreeSet<usize> =
            a.spots[0].points.iter().map(|p| p.w).collect();
        assert_eq!(used.len(), cooke().wavelengths.len());
    }

    #[test]
    fn refocusing_reports_the_thickness_it_changed() {
        let mut spec = cooke();
        let last = spec.surfaces.len() - 2;
        spec.surfaces[last].thickness += 3.0; // throw the image plane out of focus

        let blurred = run(spec.clone());
        let refocused = analyze(&Request {
            system: spec,
            rays_per_fan: 9,
            spot_grid: 15,
            refocus: true,
        })
        .unwrap();

        assert!(
            refocused.spots[0].rms < blurred.spots[0].rms * 0.5,
            "refocusing did not sharpen the axial spot: {} -> {}",
            blurred.spots[0].rms,
            refocused.spots[0].rms
        );
        // The echoed prescription must show the new thickness, so the editor can display it.
        assert!(
            (refocused.system.surfaces[last].thickness - 42.436702).abs() < 1e-3,
            "echoed thickness {}",
            refocused.system.surfaces[last].thickness
        );
    }

    #[test]
    fn glass_names_are_resolved_leniently_but_not_silently() {
        assert_eq!(material_for("").unwrap(), Material::Vacuum);
        assert_eq!(material_for("air").unwrap(), Material::Vacuum);
        assert_eq!(material_for("  AIR  ").unwrap(), Material::Vacuum);
        assert_eq!(material_for("Mirror").unwrap(), Material::Mirror);
        assert_eq!(material_for("n-bk7").unwrap(), catalog::N_BK7);

        let err = material_for("UNOBTAINIUM").unwrap_err();
        assert!(err.contains("UNOBTAINIUM"), "unhelpful error: {err}");
    }

    #[test]
    fn the_glass_list_offers_the_pseudo_media_first() {
        let names = glass_names();
        assert_eq!(names[0], "AIR");
        assert_eq!(names[1], "MIRROR");
        assert!(names.contains(&"N-SK16".to_string()));
    }

    #[test]
    fn malformed_systems_are_rejected_with_a_readable_reason() {
        let mut too_short = cooke();
        too_short.surfaces.truncate(1);
        assert!(too_short.build().unwrap_err().contains("image plane"));

        let mut no_light = cooke();
        no_light.wavelengths.clear();
        assert!(no_light.build().unwrap_err().contains("wavelength"));

        let mut shut = cooke();
        shut.entrance_pupil_diameter = 0.0;
        assert!(shut.build().unwrap_err().contains("entrance pupil"));

        let mut bad_glass = cooke();
        bad_glass.surfaces[0].glass = "PLASTIC".into();
        let err = bad_glass.build().unwrap_err();
        assert!(
            err.contains("surface 1") && err.contains("PLASTIC"),
            "{err}"
        );
    }

    #[test]
    fn a_system_with_no_stop_falls_back_to_the_first_surface() {
        let mut spec = cooke();
        for s in &mut spec.surfaces {
            s.stop = false;
        }
        let sys = spec.build().unwrap();
        assert_eq!(sys.stop_index(), Some(0));
    }

    #[test]
    fn a_null_radius_means_a_flat_surface() {
        let mut spec = cooke();
        spec.surfaces[0].radius = None;
        let sys = spec.build().unwrap();
        assert!(matches!(sys.surfaces[0].profile, Profile::Plane));
    }

    #[test]
    fn a_conic_constant_is_carried_through() {
        let mut spec = cooke();
        spec.surfaces[0].conic = -0.75;
        let sys = spec.build().unwrap();
        match sys.surfaces[0].profile {
            Profile::Conic { conic, .. } => assert_eq!(conic, -0.75),
            ref other => panic!("expected a conic, got {other:?}"),
        }
        // And it must survive the trip back to the editor.
        assert_eq!(spec_of(&sys).surfaces[0].conic, -0.75);
    }

    #[test]
    fn vignetting_is_reported_rather_than_hidden() {
        let mut spec = cooke();
        spec.surfaces[4].semi_diameter = Some(1.2); // choke the stop
        let a = run(spec);
        assert!(
            a.warnings.iter().any(|w| w.contains("vignetted")),
            "no vignetting warning: {:?}",
            a.warnings
        );
        assert!(a.spots.iter().any(|s| s.throughput < 0.98));
    }

    #[test]
    fn the_singlet_preset_also_analyses() {
        let spec: SystemSpec = serde_json::from_value(presets()["presets"][1].clone()).unwrap();
        let a = run(spec);
        assert!(a.ok);
        assert!((a.first_order.efl - 97.5804).abs() < 1e-3);
        assert_eq!(a.layout.elements.len(), 1);
    }

    #[test]
    fn bad_json_produces_a_response_rather_than_a_panic() {
        for bad in ["", "{", "null", r#"{"system": 3}"#, "[]"] {
            let out: serde_json::Value = serde_json::from_str(&analyze_json(bad)).unwrap();
            assert_eq!(out["ok"], false, "input {bad:?}");
            assert!(out["error"].is_string());
        }
    }

    #[test]
    fn a_mirror_system_analyses_end_to_end() {
        let spec = SystemSpec {
            title: "Concave mirror".into(),
            entrance_pupil_diameter: 20.0,
            wavelengths: vec![0.5875618],
            fields: vec![0.0],
            surfaces: vec![
                SurfaceSpec {
                    radius: Some(-200.0),
                    thickness: -100.0,
                    glass: "MIRROR".into(),
                    semi_diameter: Some(15.0),
                    conic: -1.0,
                    stop: true,
                    label: "primary".into(),
                },
                SurfaceSpec {
                    radius: None,
                    thickness: 0.0,
                    glass: "AIR".into(),
                    semi_diameter: None,
                    conic: 0.0,
                    stop: false,
                    label: "image".into(),
                },
            ],
        };
        let a = run(spec);
        assert!(
            (a.first_order.efl + 100.0).abs() < 1e-6,
            "EFL {}",
            a.first_order.efl
        );
        // A paraboloid is free of spherical aberration on axis.
        assert!(a.spots[0].rms < 1e-6, "RMS {} um", a.spots[0].rms);
        assert!(
            a.layout.elements.is_empty(),
            "a mirror is not a glass element"
        );
    }
}
