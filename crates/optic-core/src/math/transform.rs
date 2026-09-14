//! Rigid-body placement of a surface in the global frame.
//!
//! In a purely centred system every transform is a translation along `z`, and the
//! compiler folds most of this away. It exists from the start because coordinate
//! breaks, tilted and decentred elements are the single most invasive feature to
//! retrofit into a tracer that assumed a common axis.

use super::{Scalar, Vec3};

/// Maps local surface coordinates to global coordinates: `global = rot * local + trans`.
///
/// `rot` is stored row-major and is always orthonormal, so the inverse is the transpose.
#[derive(Clone, Copy, Debug)]
pub struct Transform<S: Scalar> {
    pub rot: [[S; 3]; 3],
    pub trans: Vec3<S>,
}

impl<S: Scalar> Transform<S> {
    pub fn identity() -> Self {
        let (o, z) = (S::one(), S::zero());
        Self {
            rot: [[o, z, z], [z, o, z], [z, z, o]],
            trans: Vec3::zero(),
        }
    }

    /// A pure shift of the vertex along the optical axis.
    pub fn along_axis(z: S) -> Self {
        Self {
            trans: Vec3::new(S::zero(), S::zero(), z),
            ..Self::identity()
        }
    }

    /// Decentre in `x`/`y` then tilt about x, y, z (radians), in that order.
    ///
    /// This matches the convention of a coordinate-break surface: decentre first,
    /// then rotate about the decentred origin.
    pub fn decenter_tilt(dx: S, dy: S, tx: S, ty: S, tz: S) -> Self {
        let rot = mul3(mul3(rot_x(tx), rot_y(ty)), rot_z(tz));
        Self {
            rot,
            trans: Vec3::new(dx, dy, S::zero()),
        }
    }

    #[inline]
    pub fn is_identity_rotation(&self) -> bool {
        self.rot[0][0].value() == 1.0
            && self.rot[1][1].value() == 1.0
            && self.rot[2][2].value() == 1.0
    }

    /// Rotate a direction from local into global coordinates.
    #[inline]
    pub fn dir_to_global(&self, v: Vec3<S>) -> Vec3<S> {
        Vec3::new(
            self.rot[0][0] * v.x + self.rot[0][1] * v.y + self.rot[0][2] * v.z,
            self.rot[1][0] * v.x + self.rot[1][1] * v.y + self.rot[1][2] * v.z,
            self.rot[2][0] * v.x + self.rot[2][1] * v.y + self.rot[2][2] * v.z,
        )
    }

    /// Rotate a direction from global into local coordinates (transpose of the above).
    #[inline]
    pub fn dir_to_local(&self, v: Vec3<S>) -> Vec3<S> {
        Vec3::new(
            self.rot[0][0] * v.x + self.rot[1][0] * v.y + self.rot[2][0] * v.z,
            self.rot[0][1] * v.x + self.rot[1][1] * v.y + self.rot[2][1] * v.z,
            self.rot[0][2] * v.x + self.rot[1][2] * v.y + self.rot[2][2] * v.z,
        )
    }

    #[inline]
    pub fn point_to_global(&self, p: Vec3<S>) -> Vec3<S> {
        self.dir_to_global(p) + self.trans
    }

    #[inline]
    pub fn point_to_local(&self, p: Vec3<S>) -> Vec3<S> {
        self.dir_to_local(p - self.trans)
    }
}

fn rot_x<S: Scalar>(a: S) -> [[S; 3]; 3] {
    let (c, s, z, o) = (a.cos(), a.sin(), S::zero(), S::one());
    [[o, z, z], [z, c, -s], [z, s, c]]
}

fn rot_y<S: Scalar>(a: S) -> [[S; 3]; 3] {
    let (c, s, z, o) = (a.cos(), a.sin(), S::zero(), S::one());
    [[c, z, s], [z, o, z], [-s, z, c]]
}

fn rot_z<S: Scalar>(a: S) -> [[S; 3]; 3] {
    let (c, s, z, o) = (a.cos(), a.sin(), S::zero(), S::one());
    [[c, -s, z], [s, c, z], [z, z, o]]
}

fn mul3<S: Scalar>(a: [[S; 3]; 3], b: [[S; 3]; 3]) -> [[S; 3]; 3] {
    let mut out = [[S::zero(); 3]; 3];
    for (i, row) in out.iter_mut().enumerate() {
        for (j, cell) in row.iter_mut().enumerate() {
            *cell = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j];
        }
    }
    out
}
