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
    solve::{self, ThicknessSolve},
    surface::Profile,
    system::{Aperture, Field, Object, Surface, Vignette, Wavelength},
    trace, Material, Paraxial, System,
};
use serde::{Deserialize, Serialize};

/// How a thickness is determined, in the form the editor exchanges.
#[derive(Clone, Copy, Debug, Default, Deserialize, Serialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum SolveSpec {
    /// Use the typed value.
    #[default]
    Fixed,
    /// Put the next surface where the paraxial marginal ray reaches `height`.
    /// With height 0 on the last thickness this is autofocus.
    MarginalRayHeight {
        #[serde(default)]
        height: f64,
    },
    /// Put the next surface where the paraxial chief ray reaches `height`.
    ChiefRayHeight {
        #[serde(default)]
        height: f64,
    },
    /// Copy surface `from`'s thickness: `scale * t + offset`.
    Pickup {
        from: usize,
        #[serde(default = "one")]
        scale: f64,
        #[serde(default)]
        offset: f64,
    },
}

fn one() -> f64 {
    1.0
}

impl SolveSpec {
    fn to_core(self) -> ThicknessSolve<f64> {
        match self {
            SolveSpec::Fixed => ThicknessSolve::Fixed,
            SolveSpec::MarginalRayHeight { height } => ThicknessSolve::MarginalRayHeight { height },
            SolveSpec::ChiefRayHeight { height } => ThicknessSolve::ChiefRayHeight { height },
            SolveSpec::Pickup {
                from,
                scale,
                offset,
            } => ThicknessSolve::Pickup {
                from,
                scale,
                offset,
            },
        }
    }

    fn from_core(s: ThicknessSolve<f64>) -> Self {
        match s {
            ThicknessSolve::Fixed => SolveSpec::Fixed,
            ThicknessSolve::MarginalRayHeight { height } => SolveSpec::MarginalRayHeight { height },
            ThicknessSolve::ChiefRayHeight { height } => SolveSpec::ChiefRayHeight { height },
            ThicknessSolve::Pickup {
                from,
                scale,
                offset,
            } => SolveSpec::Pickup {
                from,
                scale,
                offset,
            },
        }
    }
}

/// A field point, with both transverse components named as a prescription names them.
///
/// The five vignetting factors say how much of the pupil this field actually uses. All
/// zero, the default, means the whole pupil, so a spec written before they existed still
/// means exactly what it meant.
#[derive(Clone, Copy, Debug, Default, Deserialize, Serialize, PartialEq)]
pub struct FieldSpec {
    #[serde(default)]
    pub x: f64,
    #[serde(default)]
    pub y: f64,
    /// Pupil decentre, in normalised pupil coordinates.
    #[serde(default)]
    pub vdx: f64,
    #[serde(default)]
    pub vdy: f64,
    /// Pupil compression: 0 keeps the full width, 0.5 keeps half of it.
    #[serde(default)]
    pub vcx: f64,
    #[serde(default)]
    pub vcy: f64,
    /// Rotation of the vignetted pupil, in degrees.
    #[serde(default)]
    pub van: f64,
}

impl FieldSpec {
    pub fn radius(&self) -> f64 {
        (self.x * self.x + self.y * self.y).sqrt()
    }

    fn vignette(&self) -> Vignette {
        Vignette {
            dx: self.vdx,
            dy: self.vdy,
            cx: self.vcx,
            cy: self.vcy,
            angle: self.van,
        }
    }
}

/// How the pupil is sampled when building a spot diagram.
#[derive(Clone, Copy, Debug, Default, Deserialize, Serialize, PartialEq)]
#[serde(rename_all = "snake_case")]
pub enum PupilPattern {
    /// A square lattice clipped to the pupil. Even coverage, visible rows and columns.
    #[default]
    Square,
    /// Rings of increasing population. The usual choice for reading a spot's shape,
    /// because the sampling density is uniform in area rather than in x and y.
    Hexapolar,
    /// A square lattice with each point jittered inside its own cell. Breaks up the
    /// lattice artefacts that make a square grid look structured, and converges on the
    /// RMS faster. The jitter is derived from the point's index, not a random number
    /// generator, so the same system always produces the same diagram.
    Dithered,
}

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
    pub solve: SolveSpec,
    #[serde(default)]
    pub label: String,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct SystemSpec {
    #[serde(default)]
    pub title: String,
    /// Entrance pupil diameter, in millimetres.
    pub entrance_pupil_diameter: f64,
    /// Wavelengths in micrometres.
    pub wavelengths: Vec<f64>,
    /// Index into `wavelengths` of the primary line: the one first-order data, the
    /// layout and the Airy radius are referred to. Defaults to the middle entry.
    #[serde(default)]
    pub primary_wavelength: Option<usize>,
    /// Field points in degrees, for an object at infinity.
    pub fields: Vec<FieldSpec>,
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
    /// How to sample the pupil for spot diagrams.
    #[serde(default)]
    pub pupil_pattern: PupilPattern,
    /// Shift the image plane by this much, in millimetres, *after* solves are applied.
    ///
    /// Keeping defocus here rather than in the prescription means exploring focus does
    /// not fight an autofocus solve, and does not edit the design to ask a question
    /// about it.
    #[serde(default)]
    pub defocus: f64,
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
    /// Index into the system's field list, for colouring.
    pub field_index: usize,
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

/// Real chief-ray height against the paraxial prediction, at one field and wavelength.
///
/// Distortion is a chief-ray property, so it is untouched by aperture, vignetting and
/// pupil sampling. That makes it the cleanest single number for comparing two tools:
/// if it disagrees, the prescriptions differ, not the settings.
#[derive(Clone, Debug, Serialize)]
pub struct DistortionPoint {
    pub field_index: usize,
    pub wavelength: f64,
    pub percent: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct Spot {
    pub field: FieldSpec,
    pub field_index: usize,
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
    /// Distortion at every field and wavelength.
    pub distortion: Vec<DistortionPoint>,
    /// The wavelength first-order data refers to, in micrometres.
    pub primary_wavelength: f64,
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
    if let Some(m) = catalog::by_name(key) {
        return Ok(m);
    }
    // Literature quotes designs by six-digit code when the original glass is obsolete.
    optic_core::glass_code(key).ok_or_else(|| {
        format!("unknown glass \"{key}\" (try a catalogue name or a code like 613585)")
    })
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
                thickness_solve: s.solve.to_core(),
                semi_diameter: s.semi_diameter.filter(|v| *v > 0.0),
                is_stop: s.stop,
                label: s.label.clone(),
            };
            // The image plane never carries glass or a thickness of its own.
            if i + 1 == self.surfaces.len() {
                surf.material = Material::Vacuum;
                surf.thickness = 0.0;
                surf.thickness_solve = ThicknessSolve::Fixed;
            }
            surfaces.push(surf);
        }

        let mut sys = System::new(Object::Infinity, surfaces)
            .titled(&self.title)
            .with_aperture(Aperture::EntrancePupilDiameter(
                self.entrance_pupil_diameter,
            ))
            .with_fields(
                self.fields
                    .iter()
                    .map(|f| Field::angle_xy(f.x, f.y).vignetted(f.vignette()))
                    .collect(),
            )
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

impl SystemSpec {
    /// Index of the primary wavelength, clamped into range.
    pub fn primary_index(&self) -> usize {
        self.primary_wavelength
            .filter(|i| *i < self.wavelengths.len())
            .unwrap_or(self.wavelengths.len() / 2)
    }

    /// The wavelength first-order data and the Airy radius are referred to.
    pub fn primary(&self) -> f64 {
        self.wavelengths[self.primary_index()]
    }
}

/// Normalised pupil coordinates for one sampling pattern.
///
/// `n` sets the density: for a square or dithered lattice it is the side of the grid,
/// for hexapolar it fixes the ring count so the total lands near the same number.
pub fn pupil_points(pattern: PupilPattern, n: usize) -> Vec<(f64, f64)> {
    let n = n.clamp(3, 81);
    match pattern {
        PupilPattern::Square | PupilPattern::Dithered => {
            let dithered = pattern == PupilPattern::Dithered;
            let mut out = Vec::with_capacity(n * n);
            for a in 0..n {
                for b in 0..n {
                    let (jx, jy) = if dithered {
                        jitter((a * n + b) as u64)
                    } else {
                        (0.0, 0.0)
                    };
                    let px = -1.0 + 2.0 * (a as f64 + 0.5 + jx) / n as f64;
                    let py = -1.0 + 2.0 * (b as f64 + 0.5 + jy) / n as f64;
                    if px * px + py * py <= 1.0 {
                        out.push((px, py));
                    }
                }
            }
            out
        }
        PupilPattern::Hexapolar => {
            let rings = (n / 2).max(1);
            let mut out = vec![(0.0, 0.0)];
            for ring in 1..=rings {
                let r = ring as f64 / rings as f64;
                let count = 6 * ring;
                for k in 0..count {
                    let theta = std::f64::consts::TAU * k as f64 / count as f64;
                    out.push((r * theta.cos(), r * theta.sin()));
                }
            }
            out
        }
    }
}

/// Deterministic sub-cell offset in `[-0.5, 0.5)^2`.
///
/// A hash of the point's index rather than a random number generator, so a dithered
/// diagram is reproducible: the same system always yields the same picture, which
/// matters when a spot diagram is evidence in a design review.
fn jitter(i: u64) -> (f64, f64) {
    let mut z = i.wrapping_add(1).wrapping_mul(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^= z >> 31;
    let a = (z & 0xFFFF_FFFF) as f64 / u32::MAX as f64 - 0.5;
    let b = ((z >> 32) & 0xFFFF_FFFF) as f64 / u32::MAX as f64 - 0.5;
    (a, b)
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
    let wl = request.system.primary();
    let mut warnings = Vec::new();

    if request.refocus {
        samples::focus(&mut sys, wl);
    }

    // Solves run before anything is traced, so every analysis below sees a system that
    // already satisfies them.
    match solve::resolve(&mut sys, wl) {
        optic_core::SolveReport::Degenerate { surface, reason } => {
            warnings.push(format!("surface {}: solve ignored, {reason}", surface + 1));
        }
        optic_core::SolveReport::NotConverged { passes } => {
            warnings.push(format!("solves did not settle after {passes} passes"));
        }
        _ => {}
    }

    // Snapshot what the design says before defocus perturbs it. Defocus is a question
    // asked *of* the design, so the editor must keep showing the design's own numbers.
    let design_thicknesses: Vec<f64> = sys.surfaces.iter().map(|s| s.thickness).collect();

    if request.defocus != 0.0 {
        let last = sys.image_index() - 1;
        sys.surfaces[last].thickness += request.defocus;
    }

    let par = Paraxial::compute(&sys, wl);
    if !par.efl.is_finite() {
        return Err("the system has no focal length: check radii and glasses".into());
    }

    let vertices = sys.vertices();
    let image_z = vertices[sys.image_index()];
    let field_specs = request.system.fields.clone();

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
                field_index: fi,
                points,
                complete: traced.is_complete(),
            });
        }
    }

    // Spot diagrams over the whole wavelength set, and distortion per wavelength.
    let pupil = pupil_points(request.pupil_pattern, request.spot_grid);
    let mut spots = Vec::new();
    let mut distortion = Vec::new();

    for (fi, field) in sys.fields.iter().enumerate() {
        let mut pts: Vec<(f64, f64, usize)> = Vec::new();
        let mut launched = 0usize;

        for (wi, w) in request.system.wavelengths.iter().enumerate() {
            let par_w = Paraxial::compute(&sys, *w);

            // Distortion: the real chief ray against the paraxial prediction. Both are
            // measured along the paraxial image point's own direction rather than by
            // magnitude, so a chief ray that lands on the far side of the axis reads as
            // the gross distortion it is, and so that this stays right once decentres
            // and tilts can push the real intercept out of the meridional plane.
            let want = par_w.image_point(&sys, *w, *field);
            let reference = (want[0] * want[0] + want[1] * want[1]).sqrt();
            if reference > 1e-9 {
                // Deliberately the *unvignetted* chief ray. Vignetting says which rays
                // we choose to sample, and distortion is a property of the lens: if a
                // decentred pupil could move it, the figure would describe our sampling
                // rather than the design.
                let chief = launch(&sys, &par_w, field.unvignetted(), 0.0, 0.0);
                if let Some(p) = trace(&sys, *w, chief).image_point() {
                    let real = (p.x * want[0] + p.y * want[1]) / reference;
                    distortion.push(DistortionPoint {
                        field_index: fi,
                        wavelength: *w,
                        percent: 100.0 * (real - reference) / reference,
                    });
                }
            }

            for (px, py) in &pupil {
                launched += 1;
                if let Some(p) =
                    trace(&sys, *w, launch(&sys, &par_w, *field, *px, *py)).image_point()
                {
                    pts.push((p.x, p.y, wi));
                }
            }
        }

        let spec = field_specs.get(fi).copied().unwrap_or_default();
        if pts.is_empty() {
            warnings.push(format!(
                "no ray reached the image plane at field {:.1}, {:.1}",
                spec.x, spec.y
            ));
            spots.push(Spot {
                field: spec,
                field_index: fi,
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
            field: spec,
            field_index: fi,
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
                "{:.0}% of rays are vignetted at field {:.1}, {:.1}",
                100.0 * (1.0 - s.throughput),
                s.field.x,
                s.field.y
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

    // Report back the prescription actually analysed, so solved thicknesses and a
    // refocus are visible in the editor rather than silently applied.
    let mut echoed = request.system.clone();
    for ((spec, surf), thickness) in echoed
        .surfaces
        .iter_mut()
        .zip(sys.surfaces.iter())
        .zip(design_thicknesses.iter())
    {
        spec.thickness = *thickness;
        spec.solve = SolveSpec::from_core(surf.thickness_solve);
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
        distortion,
        primary_wavelength: wl,
        warnings,
        system: echoed,
    })
}

/// Read a `.zmx` file's raw bytes into a prescription the editor can hold.
///
/// Takes bytes rather than a string because the encoding is part of the problem: Zemax
/// has written both UTF-16LE and UTF-8, and decoding is the reader's job.
pub fn import_zmx_json(bytes: &[u8]) -> String {
    let text = optic_io::zmx::decode(bytes);
    match optic_io::zmx::parse(&text) {
        Ok(imported) => {
            let mut spec = spec_of(&imported.system);
            if let Some(p) = optic_io::zmx::primary_index(&text) {
                if p < spec.wavelengths.len() {
                    spec.primary_wavelength = Some(p);
                }
            }
            serde_json::json!({
                "ok": true,
                "system": spec,
                "warnings": imported.warnings,
            })
            .to_string()
        }
        Err(e) => error_json(&e),
    }
}

/// Write a prescription out as `.zmx`, so it can be checked in Zemax.
pub fn export_zmx_json(request: &str) -> String {
    let parsed: Result<SystemSpec, _> = serde_json::from_str(request);
    match parsed
        .map_err(|e| e.to_string())
        .and_then(|spec| spec.build())
    {
        Ok(sys) => {
            serde_json::json!({ "ok": true, "text": optic_io::zmx::write(&sys) }).to_string()
        }
        Err(e) => error_json(&e),
    }
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
        "presets": [
            spec_of(&samples::cooke_triplet::<f64>()),
            spec_of(&samples::singlet::<f64>()),
            spec_of(&samples::smith_triplet_moderate::<f64>()),
            spec_of(&samples::smith_triplet_wide::<f64>()),
        ],
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
                solve: SolveSpec::from_core(s.thickness_solve),
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
        primary_wavelength: Some(sys.wavelengths.len() / 2),
        fields: sys
            .fields
            .iter()
            .map(|f| {
                let v = f.vignette();
                let (x, y) = match *f {
                    Field::Angle { x, y, .. } | Field::Height { x, y, .. } => (x, y),
                };
                FieldSpec {
                    x,
                    y,
                    vdx: v.dx,
                    vdy: v.dy,
                    vcx: v.cx,
                    vcy: v.cy,
                    van: v.angle,
                }
            })
            .collect(),
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

/// Name a medium so that [`material_for`] reads back the same medium.
///
/// A catalogue name is preferred because it carries the full dispersion. Failing that,
/// a six-digit code preserves n_d and V_d, which is what a design quoted from the
/// literature has anyway. Falling back to "AIR" would silently strip a lens of its
/// glass and leave a prescription with no focal length, so anything unnameable is
/// reported as such and refused on the way back in.
fn glass_name_of(m: &Material) -> String {
    match m {
        Material::Vacuum | Material::Air => "AIR".into(),
        Material::Mirror => "MIRROR".into(),
        other => catalog::PUBLISHED
            .iter()
            .find(|(_, g, ..)| *g == other)
            .map(|(n, ..)| n.to_string())
            .or_else(|| other.code())
            .unwrap_or_else(|| "UNNAMEABLE".into()),
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

    #[test]
    fn every_preset_survives_the_journey_to_the_browser() {
        // A preset reaches the demo as JSON, so a medium the projection cannot name is a
        // medium the browser never receives. Naming an unknown glass "AIR" once turned
        // both published triplets into stacks of air: they loaded with no focal length
        // at all, while the Rust samples they were built from were perfectly correct.
        // Checking the samples is therefore not enough; the projection needs its own test.
        let sources = [
            samples::cooke_triplet::<f64>(),
            samples::singlet::<f64>(),
            samples::smith_triplet_moderate::<f64>(),
            samples::smith_triplet_wide::<f64>(),
        ];
        let presets = presets();
        let presets = presets["presets"].as_array().unwrap();
        assert_eq!(presets.len(), sources.len(), "a preset was added untested");

        for (preset, source) in presets.iter().zip(&sources) {
            let spec: SystemSpec = serde_json::from_value(preset.clone()).unwrap();
            let title = spec.title.clone();

            for s in &spec.surfaces {
                assert!(
                    material_for(&s.glass).is_ok(),
                    "{title}: glass {:?} cannot be read back",
                    s.glass
                );
            }

            // Against the sample itself, not against the projection's own output, so a
            // glass that degrades into a nearby one is caught as well as one that
            // vanishes entirely.
            let primary = source.wavelengths[source.wavelengths.len() / 2].um;
            let want = Paraxial::compute(source, primary).efl;
            let got = run(spec).first_order.efl;
            assert!(
                (got - want).abs() < 1e-9,
                "{title}: focal length {got} through JSON, {want} in the sample"
            );
        }
    }

    fn run(spec: SystemSpec) -> Analysis {
        analyze(&Request {
            system: spec,
            rays_per_fan: 9,
            spot_grid: 15,
            refocus: false,
            pupil_pattern: PupilPattern::Square,
            defocus: 0.0,
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
                .filter(|r| r.field_index == 0)
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
    fn a_solved_thickness_ignores_a_typed_value() {
        // The Cooke preset autofocuses. Typing a different image distance must have no
        // effect, because the solve owns that number -- and the editor must be told the
        // value that was actually used, not the one that was typed.
        let mut spec = cooke();
        let last = spec.surfaces.len() - 2;
        assert!(
            !matches!(spec.surfaces[last].solve, SolveSpec::Fixed),
            "the preset should ship with a solve"
        );

        let focused = run(spec.clone());
        spec.surfaces[last].thickness += 3.0;
        let edited = run(spec.clone());

        assert!(
            (edited.spots[0].rms - focused.spots[0].rms).abs() < 1e-9,
            "the typed thickness defeated the solve: {} vs {}",
            focused.spots[0].rms,
            edited.spots[0].rms
        );
        assert!(
            (edited.system.surfaces[last].thickness - 42.436702).abs() < 1e-3,
            "the echoed thickness is the typed one, not the solved one: {}",
            edited.system.surfaces[last].thickness
        );
    }

    #[test]
    fn releasing_the_solve_hands_the_thickness_back() {
        let mut spec = cooke();
        let last = spec.surfaces.len() - 2;
        spec.surfaces[last].solve = SolveSpec::Fixed;
        spec.surfaces[last].thickness = 45.0;

        let a = run(spec);
        assert!((a.system.surfaces[last].thickness - 45.0).abs() < 1e-12);
        assert!(a.spots[0].rms > 100.0, "45 mm should be badly defocused");
    }

    #[test]
    fn defocus_blurs_without_editing_the_prescription() {
        // Exploring focus must not alter the design, and must not fight the solve.
        let spec = cooke();
        let sharp = run(spec.clone());
        let shifted = analyze(&Request {
            system: spec.clone(),
            rays_per_fan: 9,
            spot_grid: 15,
            refocus: false,
            pupil_pattern: PupilPattern::Square,
            defocus: 1.5,
        })
        .unwrap();

        assert!(
            shifted.spots[0].rms > sharp.spots[0].rms * 2.0,
            "1.5 mm of defocus barely changed the spot: {} -> {}",
            sharp.spots[0].rms,
            shifted.spots[0].rms
        );
        let last = spec.surfaces.len() - 2;
        assert_eq!(
            shifted.system.surfaces[last].thickness, sharp.system.surfaces[last].thickness,
            "defocus edited the prescription"
        );
    }

    #[test]
    fn a_solve_holds_focus_when_the_design_changes() {
        let mut spec = cooke();
        let sharp = run(spec.clone());
        spec.surfaces[0].radius = Some(spec.surfaces[0].radius.unwrap() * 1.05);
        let moved = run(spec);
        // A 5% radius change would ruin an unsolved system; with autofocus the axial
        // spot stays in the same order of magnitude.
        assert!(
            moved.spots[0].rms < sharp.spots[0].rms * 3.0,
            "autofocus did not hold: {} -> {}",
            sharp.spots[0].rms,
            moved.spots[0].rms
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
    fn a_prescription_survives_export_and_reimport_through_the_api() {
        let original = cooke();
        let exported: serde_json::Value =
            serde_json::from_str(&export_zmx_json(&serde_json::to_string(&original).unwrap()))
                .unwrap();
        assert_eq!(exported["ok"], true);
        let text = exported["text"].as_str().expect("zmx text");

        let imported: serde_json::Value =
            serde_json::from_str(&import_zmx_json(text.as_bytes())).unwrap();
        assert_eq!(imported["ok"], true, "{imported}");
        let spec: SystemSpec = serde_json::from_value(imported["system"].clone()).unwrap();

        assert_eq!(spec.surfaces.len(), original.surfaces.len());
        assert_eq!(spec.wavelengths.len(), original.wavelengths.len());
        for (a, b) in spec.surfaces.iter().zip(original.surfaces.iter()) {
            assert_eq!(a.glass, b.glass);
            assert!((a.thickness - b.thickness).abs() < 1e-9);
        }

        // And the round-tripped system analyses to the same first-order data.
        let before = run(original).first_order.efl;
        let after = run(spec).first_order.efl;
        assert!((before - after).abs() < 1e-6, "{before} vs {after}");
    }

    #[test]
    fn importing_a_utf16_file_works_through_the_api() {
        let text = export_zmx_json(&serde_json::to_string(&cooke()).unwrap());
        let inner: serde_json::Value = serde_json::from_str(&text).unwrap();
        let zmx = inner["text"].as_str().unwrap();

        let mut bytes = vec![0xFF, 0xFE];
        for unit in zmx.encode_utf16() {
            bytes.extend_from_slice(&unit.to_le_bytes());
        }
        let imported: serde_json::Value = serde_json::from_str(&import_zmx_json(&bytes)).unwrap();
        assert_eq!(imported["ok"], true, "{imported}");
    }

    #[test]
    fn importing_something_that_is_not_a_prescription_fails_politely() {
        for junk in [&b""[..], b"hello", b"\x00\x01\x02\x03"] {
            let out: serde_json::Value = serde_json::from_str(&import_zmx_json(junk)).unwrap();
            assert_eq!(out["ok"], false, "{out}");
            assert!(out["error"].is_string());
        }
    }

    #[test]
    fn import_carries_vignetting_factors_onto_the_field() {
        let zmx = "VERS 1\nMODE SEQ\nNAME t\nENPD 10\nFTYP 0 0 1 1 0 0 0 1\nXFLN 0\nYFLN 0\n\
                   WAVM 1 0.5876 1\nVDYN 0.3\nSURF 0\n  TYPE STANDARD\n  CURV 0\n  DISZ INFINITY\n\
                   SURF 1\n  TYPE STANDARD\n  CURV 0.01\n  DISZ 5\n  GLAS N-BK7 0 0 1.5168 64.17\n  STOP\n\
                   SURF 2\n  TYPE STANDARD\n  CURV -0.01\n  DISZ 95\n\
                   SURF 3\n  TYPE STANDARD\n  CURV 0\n  DISZ 0\n";
        let out: serde_json::Value =
            serde_json::from_str(&import_zmx_json(zmx.as_bytes())).unwrap();
        assert_eq!(out["ok"], true, "{out}");
        assert_eq!(out["system"]["fields"][0]["vdy"], 0.3, "{out}");

        // They are honoured now, so nothing should be reported as unsupported.
        let warnings = out["warnings"].as_array().unwrap();
        assert!(
            !warnings
                .iter()
                .any(|w| w.as_str().unwrap().contains("vignetting")),
            "{warnings:?}"
        );
    }

    #[test]
    fn vignetting_shrinks_the_spot_without_moving_the_distortion() {
        // The two halves of the promise. Sampling less of the pupil must drop the rays
        // that carry the most aberration, so the spot shrinks. Distortion is a chief-ray
        // property, so it must not move at all -- which is what lets an outer-field
        // distortion figure be compared against Zemax even when the vignetting settings
        // are not yet known to match.
        let mut spec = cooke();
        let wide = spec.fields.len() - 1;

        let before = run(spec.clone());
        spec.fields[wide].vcy = 0.5;
        spec.fields[wide].vcx = 0.5;
        let after = run(spec);

        assert!(
            after.spots[wide].rms < before.spots[wide].rms * 0.9,
            "{} against {}",
            after.spots[wide].rms,
            before.spots[wide].rms
        );
        assert!(
            (after.spots[0].rms - before.spots[0].rms).abs() < 1e-12,
            "an unvignetted field was disturbed"
        );

        for (a, b) in after.distortion.iter().zip(&before.distortion) {
            assert!(
                (a.percent - b.percent).abs() < 1e-12,
                "distortion moved: {} against {}",
                a.percent,
                b.percent
            );
        }
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
            primary_wavelength: None,
            fields: vec![FieldSpec::default()],
            surfaces: vec![
                SurfaceSpec {
                    radius: Some(-200.0),
                    thickness: -100.0,
                    glass: "MIRROR".into(),
                    semi_diameter: Some(15.0),
                    conic: -1.0,
                    stop: true,
                    solve: SolveSpec::Fixed,
                    label: "primary".into(),
                },
                SurfaceSpec {
                    radius: None,
                    thickness: 0.0,
                    glass: "AIR".into(),
                    semi_diameter: None,
                    conic: 0.0,
                    stop: false,
                    solve: SolveSpec::Fixed,
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
