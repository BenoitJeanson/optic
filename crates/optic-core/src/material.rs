//! Optical media and their dispersion.
//!
//! Coefficients here are *physical constants*, not design variables, so they are plain
//! `f64` rather than [`Scalar`](crate::Scalar). The one case that wants differentiation
//! is model glass during glass optimisation; that will be handled by promoting
//! `ModelGlass` at the variable layer rather than by making every catalog entry generic.
//!
//! Wavelengths are micrometres throughout, matching every published dispersion formula.

/// Spectral lines, in micrometres.
pub mod lines {
    /// Helium d line, 587.5618 nm. The reference for `n_d` and Abbe `V_d`.
    pub const D: f64 = 0.5875618;
    /// Hydrogen F line, 486.1327 nm.
    pub const F: f64 = 0.4861327;
    /// Hydrogen C line, 656.2725 nm.
    pub const C: f64 = 0.6562725;
    /// Mercury e line, 546.0740 nm.
    pub const E: f64 = 0.5460740;
    /// Mercury g line, 435.8343 nm.
    pub const G: f64 = 0.4358343;
}

/// A medium between two surfaces.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(tag = "model", rename_all = "snake_case"))]
pub enum Material {
    /// Index 1 at all wavelengths.
    Vacuum,
    /// Standard air at 20 C, 1 atm (Kohlrausch/Edlen fit).
    Air,
    /// Reflective. The medium is unchanged; the ray direction is folded.
    Mirror,
    /// Non-dispersive medium of fixed index.
    Fixed { n: f64 },
    /// Two-term Cauchy fit through `n_d` and the Abbe number.
    ///
    /// Used to let the optimiser move continuously through glass space before
    /// substituting a real catalog glass. Valid in the visible only.
    ModelGlass { nd: f64, vd: f64 },
    /// `n^2 - 1 = sum_i b_i L^2 / (L^2 - c_i)`. The near-universal catalog format.
    Sellmeier1 { b: [f64; 3], c: [f64; 3] },
    /// `n^2 = a0 + a1 L^2 + a2 L^-2 + a3 L^-4 + a4 L^-6 + a5 L^-8`. Older Schott format.
    Schott { a: [f64; 6] },
}

impl Material {
    /// Refractive index at `wavelength` micrometres.
    pub fn index(&self, wavelength: f64) -> f64 {
        let l2 = wavelength * wavelength;
        match *self {
            Material::Vacuum | Material::Mirror => 1.0,
            Material::Air => air_index(wavelength),
            Material::Fixed { n } => n,
            Material::ModelGlass { nd, vd } => {
                // n = a + b / L^2, pinned by n_d and by n_F - n_C = (n_d - 1) / V_d.
                let (fd, ff, fc) = (
                    1.0 / (lines::D * lines::D),
                    1.0 / (lines::F * lines::F),
                    1.0 / (lines::C * lines::C),
                );
                let b = (nd - 1.0) / (vd * (ff - fc));
                let a = nd - b * fd;
                a + b / l2
            }
            Material::Sellmeier1 { b, c } => {
                let n2 = 1.0
                    + b[0] * l2 / (l2 - c[0])
                    + b[1] * l2 / (l2 - c[1])
                    + b[2] * l2 / (l2 - c[2]);
                n2.sqrt()
            }
            Material::Schott { a } => {
                let n2 = a[0]
                    + a[1] * l2
                    + a[2] / l2
                    + a[3] / l2.powi(2)
                    + a[4] / l2.powi(3)
                    + a[5] / l2.powi(4);
                n2.sqrt()
            }
        }
    }

    /// Index at the helium d line.
    pub fn nd(&self) -> f64 {
        self.index(lines::D)
    }

    /// Abbe number `V_d = (n_d - 1) / (n_F - n_C)`.
    ///
    /// Returns [`f64::INFINITY`] for a non-dispersive medium.
    pub fn abbe(&self) -> f64 {
        let d = self.index(lines::F) - self.index(lines::C);
        if d == 0.0 {
            f64::INFINITY
        } else {
            (self.nd() - 1.0) / d
        }
    }

    /// Whether a ray meeting this medium is reflected rather than refracted.
    pub fn is_reflective(&self) -> bool {
        matches!(self, Material::Mirror)
    }
}

/// Refractive index of standard air, Kohlrausch's form of the Edlen equation.
///
/// Catalog glass indices are quoted *relative to air*, so tracing a lens in `Air`
/// rather than `Vacuum` double-counts this. It is here for completeness and for the
/// eventual environmental model; centred visible-light designs should use [`Material::Vacuum`]
/// for the surrounding medium to match catalog conventions.
fn air_index(wavelength: f64) -> f64 {
    let s2 = 1.0 / (wavelength * wavelength);
    1.0 + (6432.8 + 2_949_810.0 / (146.0 - s2) + 25_540.0 / (41.0 - s2)) * 1e-8
}

/// A small built-in catalog.
///
/// These are Schott N-series coefficients. The test suite checks each entry reproduces
/// its published `n_d` and `V_d`, so a transcription error fails the build rather than
/// silently biasing every design.
pub mod catalog {
    use super::Material;

    macro_rules! glass {
        ($name:ident, $b:expr, $c:expr) => {
            pub const $name: Material = Material::Sellmeier1 { b: $b, c: $c };
        };
    }

    glass!(
        N_BK7,
        [1.03961212, 0.231792344, 1.01046945],
        [0.00600069867, 0.0200179144, 103.560653]
    );
    glass!(
        N_SK16,
        [1.34317774, 0.241144399, 0.994317969],
        [0.00704687339, 0.0229005, 92.7508526]
    );
    glass!(
        F2,
        [1.34533359, 0.209073176, 0.937357162],
        [0.00997743871, 0.0470450767, 111.886764]
    );
    glass!(
        N_SF11,
        [1.73759695, 0.313747346, 1.89878101],
        [0.013188707, 0.0623068142, 155.23629]
    );
    glass!(
        N_BAF10,
        [1.5851495, 0.143559385, 1.08521269],
        [0.00926681282, 0.0424489805, 105.613573]
    );

    /// Published `(n_d, V_d)` for each entry above, used by the verification test.
    pub const PUBLISHED: &[(&str, &Material, f64, f64)] = &[
        ("N-BK7", &N_BK7, 1.5168, 64.17),
        ("N-SK16", &N_SK16, 1.62041, 60.32),
        ("F2", &F2, 1.62004, 36.37),
        ("N-SF11", &N_SF11, 1.78472, 25.68),
        ("N-BAF10", &N_BAF10, 1.67003, 47.11),
    ];

    /// Look a glass up by name, case- and separator-insensitively.
    pub fn by_name(name: &str) -> Option<Material> {
        let key: String = name
            .chars()
            .filter(|c| c.is_ascii_alphanumeric())
            .map(|c| c.to_ascii_uppercase())
            .collect();
        PUBLISHED
            .iter()
            .find(|(n, ..)| {
                let k: String = n
                    .chars()
                    .filter(|c| c.is_ascii_alphanumeric())
                    .map(|c| c.to_ascii_uppercase())
                    .collect();
                k == key
            })
            .map(|(_, m, ..)| (*m).clone())
    }
}
