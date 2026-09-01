//! Spatial inertia and rigid-body types.
//!
//! Provides:
//! - [`SpatialInertia`] — 6×6 spatial mass matrix capturing mass, center of
//!   mass, and rotational inertia in a single entity (Featherstone convention).
//! - [`RigidBody`] — spatial inertia + reference frame + collision geometry.
//!
//! Cross-product operator structs:
//! - [`MotionCross`] — `v×` applied to a [`SpatialVelocity`].
//! - [`ForceCross`] — `v×*` applied to a [`SpatialForce`].
//!
//! # Examples
//!
//! ```
//! use tpt_mbd_core::{SpatialInertia, RigidBody, Frame};
//! use tpt_mbd_core::spatial::SpatialVelocity;
//! use tpt_math_linalg_fixed::{Matrix3, Vector3};
//!
//! let si = SpatialInertia::new(
//!     2.0,
//!     Vector3::new([0.0, 0.0, 0.0]),
//!     Matrix3::new([
//!         [1.0, 0.0, 0.0],
//!         [0.0, 2.0, 0.0],
//!         [0.0, 0.0, 3.0],
//!     ]),
//! );
//! assert_eq!(si.mass, 2.0);
//!
//! let body = RigidBody::new(si, Frame::identity(), "link0", 0);
//! assert_eq!(body.name, "link0");
//! ```

use tpt_math_geometry::Isometry3;
use tpt_math_linalg_fixed::{Matrix, Matrix3, Vector3, Vector6};
use tpt_math_numeric::Scalar;

use crate::spatial::{force_cross, motion_cross, SpatialForce, SpatialMomentum, SpatialVelocity};

type Matrix6<T> = Matrix<T, 6, 6>;

// ===========================================================================
// SpatialInertia
// ===========================================================================

/// 6×6 spatial mass matrix capturing mass, center of mass, and rotational inertia.
///
/// Encodes a rigid body's inertia in Featherstone's spatial vector convention:
/// the upper-left 3×3 block is the mass, the upper-right is the COM skew term,
/// and the lower-right is the rotational inertia about the COM.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpatialInertia<T: Copy> {
    /// Full 6×6 spatial mass matrix.
    pub matrix: Matrix6<T>,
    /// Total mass.
    pub mass: T,
    /// Center of mass position relative to the reference frame origin.
    pub com: Vector3<T>,
    /// Rotational inertia about the center of mass (3×3).
    pub inertia_com: Matrix3<T>,
}

impl<T: Copy + Scalar> SpatialInertia<T> {
    /// Construct from mass, center-of-mass offset, and rotational inertia about the COM.
    pub fn new(mass: T, com: Vector3<T>, inertia_com: Matrix3<T>) -> Self {
        let m = mass;
        let z = T::from(0.0).unwrap();
        let o = T::one();

        let cx = skew_matrix(&com);
        let c2 = com.data[0] * com.data[0] + com.data[1] * com.data[1] + com.data[2] * com.data[2];
        let c00 = com.data[0] * com.data[0];
        let c01 = com.data[0] * com.data[1];
        let c02 = com.data[0] * com.data[2];
        let c10 = com.data[1] * com.data[0];
        let c11 = com.data[1] * com.data[1];
        let c12 = com.data[1] * com.data[2];
        let c20 = com.data[2] * com.data[0];
        let c21 = com.data[2] * com.data[1];
        let c22 = com.data[2] * com.data[2];

        let i00 = inertia_com.data[0][0] + m * (o * c2 - c00);
        let i01 = inertia_com.data[0][1] + m * (z - c01);
        let i02 = inertia_com.data[0][2] + m * (z - c02);
        let i10 = inertia_com.data[1][0] + m * (z - c10);
        let i11 = inertia_com.data[1][1] + m * (o * c2 - c11);
        let i12 = inertia_com.data[1][2] + m * (z - c12);
        let i20 = inertia_com.data[2][0] + m * (z - c20);
        let i21 = inertia_com.data[2][1] + m * (z - c21);
        let i22 = inertia_com.data[2][2] + m * (o * c2 - c22);

        let matrix = Matrix6::new([
            [m, z, z, z, -m * cx.data[0][2], m * cx.data[0][1]],
            [z, m, z, m * cx.data[0][2], z, -m * cx.data[0][0]],
            [z, z, m, -m * cx.data[1][2], m * cx.data[1][1], z],
            [z, m * cx.data[0][2], -m * cx.data[1][2], i00, i01, i02],
            [-m * cx.data[0][2], z, m * cx.data[1][0], i10, i11, i12],
            [m * cx.data[0][1], -m * cx.data[1][0], z, i20, i21, i22],
        ]);

        SpatialInertia {
            matrix,
            mass: m,
            com,
            inertia_com,
        }
    }

    /// Transform this inertia to a new reference point displaced by `d` from the current COM.
    pub fn inertia_about(&self, d: &Vector3<T>) -> Self {
        let m = self.mass;
        let z = T::from(0.0).unwrap();
        let o = T::one();

        let dx = skew_matrix(d);
        let dxt = dx.transpose();
        let cx = skew_matrix(&self.com);

        let mut dcd = [[z; 3]; 3];
        let mut ddd = [[z; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    dcd[i][j] = dcd[i][j] + dxt.data[i][k] * cx.data[k][j];
                    ddd[i][j] = ddd[i][j] + dxt.data[i][k] * dx.data[k][j];
                }
            }
        }

        let d2 = d.data[0] * d.data[0] + d.data[1] * d.data[1] + d.data[2] * d.data[2];
        let d2m = dx * dx;

        let mut nm = [[z; 6]; 6];
        for i in 0..3 {
            for j in 0..3 {
                nm[i][j] = self.matrix.data[i][j] - m * d2m.data[i][j];
                nm[i][3 + j] = self.matrix.data[i][3 + j] - m * dcd[i][j];
                nm[3 + i][j] = self.matrix.data[3 + i][j] - m * dcd[j][i];
                nm[3 + i][3 + j] = self.matrix.data[3 + i][3 + j] + m * ddd[i][j];
            }
        }

        let new_com = Vector3::new([
            self.com.data[0] - d.data[0],
            self.com.data[1] - d.data[1],
            self.com.data[2] - d.data[2],
        ]);
        let mut new_ic = self.inertia_com;
        for i in 0..3 {
            for j in 0..3 {
                let diag = if i == j { o } else { z };
                new_ic.data[i][j] =
                    self.inertia_com.data[i][j] + m * (diag * d2 - d.data[i] * d.data[j]);
            }
        }

        SpatialInertia {
            matrix: Matrix6::new(nm),
            mass: m,
            com: new_com,
            inertia_com: new_ic,
        }
    }

    /// Compose (add) two inertias into an equivalent single inertia.
    pub fn compose(&self, other: &SpatialInertia<T>) -> Self {
        let mut new_matrix = self.matrix;
        for i in 0..6 {
            for j in 0..6 {
                new_matrix.data[i][j] = new_matrix.data[i][j] + other.matrix.data[i][j];
            }
        }
        let total_mass = self.mass + other.mass;
        let new_com = if total_mass == T::from(0.0).unwrap() {
            Vector3::new([z(); 3])
        } else {
            let inv: T = T::one() / total_mass;
            Vector3::new([
                (self.mass * self.com.data[0] + other.mass * other.com.data[0]) * inv,
                (self.mass * self.com.data[1] + other.mass * other.com.data[1]) * inv,
                (self.mass * self.com.data[2] + other.mass * other.com.data[2]) * inv,
            ])
        };
        let mut new_ic = self.inertia_com;
        for i in 0..3 {
            for j in 0..3 {
                new_ic.data[i][j] = new_ic.data[i][j] + other.inertia_com.data[i][j];
            }
        }
        SpatialInertia {
            matrix: new_matrix,
            mass: total_mass,
            com: new_com,
            inertia_com: new_ic,
        }
    }

    /// Transform this inertia into a new coordinate frame using the 6×6 adjoint of `x`.
    pub fn transform(&self, x: &Isometry3<T>) -> Self {
        let ad = spatial_inertia_transform(x);
        let z = T::from(0.0).unwrap();
        let mut new_matrix = Matrix6::new([[z; 6]; 6]);
        for i in 0..6 {
            for j in 0..6 {
                let mut val = z;
                for k in 0..6 {
                    for l in 0..6 {
                        val = val + ad.data[i][k] * self.matrix.data[k][l] * ad.data[j][l];
                    }
                }
                new_matrix.data[i][j] = val;
            }
        }
        let new_com = x.rotation.transform_vector(&self.com);
        SpatialInertia {
            matrix: new_matrix,
            mass: self.mass,
            com: new_com,
            inertia_com: self.inertia_com,
        }
    }

    /// Compute spatial momentum `h = I * v`.
    pub fn momentum(&self, velocity: &SpatialVelocity<T>) -> SpatialMomentum<T> {
        let d = velocity.0.data;
        let v6 = Vector6::new([
            d.data[0], d.data[1], d.data[2], d.data[3], d.data[4], d.data[5],
        ]);
        let data = self.matrix.mul_vec(&v6);
        SpatialMomentum::new(
            Vector3::new([data.data[0], data.data[1], data.data[2]]),
            Vector3::new([data.data[3], data.data[4], data.data[5]]),
        )
    }

    /// Alias for [`Self::new`].
    pub fn from_rotational(mass: T, com: Vector3<T>, inertia_com: Matrix3<T>) -> Self {
        Self::new(mass, com, inertia_com)
    }
}

#[inline]
fn z<T: Copy + Scalar>() -> T {
    T::from(0.0).unwrap()
}

// ===========================================================================
// RigidBody
// ===========================================================================

/// A rigid body defined by its spatial inertia, reference pose, and collision geometry.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidBody<T: Copy + Scalar> {
    /// Spatial inertia tensor.
    pub spatial_inertia: SpatialInertia<T>,
    /// World-frame pose (position + orientation).
    pub transform: Isometry3<T>,
    /// Body identifier.
    pub name: &'static str,
    /// Collision geometry index (for contact detection).
    pub collision_geometry: usize,
}

impl<T: Copy + Scalar> RigidBody<T> {
    /// Construct a new rigid body definition.
    pub fn new(
        spatial_inertia: SpatialInertia<T>,
        transform: Isometry3<T>,
        name: &'static str,
        collision_geometry: usize,
    ) -> Self {
        RigidBody {
            spatial_inertia,
            transform,
            name,
            collision_geometry,
        }
    }
}

// ===========================================================================
// Cross-product operator structs
// ===========================================================================

/// Motion cross-product operator `v×` (Featherstone `crm`).
///
/// Encapsulates an angular velocity vector and applies the cross-product to a
/// spatial velocity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MotionCross<T: Copy> {
    /// Angular velocity vector defining the cross-product operator.
    pub angular: Vector3<T>,
}

impl<T: Copy + Scalar> MotionCross<T> {
    /// Build from an angular velocity vector.
    pub fn new(angular: Vector3<T>) -> Self {
        MotionCross { angular }
    }
    /// Apply the motion cross-product operator `v×` to a spatial velocity.
    pub fn apply(&self, v: &SpatialVelocity<T>) -> SpatialVelocity<T> {
        motion_cross(
            &SpatialVelocity::new(self.angular, Vector3::new([z(); 3])),
            v,
        )
    }
}

/// Force cross-product operator `v×*` (Featherstone `crf`).
///
/// Encapsulates an angular (moment) vector and applies the cross-product to a
/// spatial force.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ForceCross<T: Copy> {
    /// Angular (moment) vector defining the cross-product operator.
    pub angular: Vector3<T>,
}

impl<T: Copy + Scalar> ForceCross<T> {
    /// Build from an angular vector.
    pub fn new(angular: Vector3<T>) -> Self {
        ForceCross { angular }
    }
    /// Apply the force cross-product operator `v×*` to a spatial force.
    pub fn apply(&self, f: &SpatialForce<T>) -> SpatialVelocity<T> {
        force_cross(&SpatialForce::new(self.angular, Vector3::new([z(); 3])), f)
    }
}

// ===========================================================================
// Free helpers
// ===========================================================================

/// Compute the skew-symmetric cross-product matrix for a vector `v`.
#[inline]
pub fn skew_matrix<T: Copy + Scalar>(v: &Vector3<T>) -> Matrix3<T> {
    let z = T::from(0.0).unwrap();
    Matrix3::new([[z, -v.z(), v.y()], [v.z(), z, -v.x()], [-v.y(), v.x(), z]])
}

#[inline]
fn spatial_inertia_transform<T: Copy + Scalar>(x: &Isometry3<T>) -> Matrix6<T> {
    let r = x.rotation.matrix();
    let t = &x.translation.vector;
    let tl = skew_matrix(t) * *r;
    let z = T::from(0.0).unwrap();
    let mut ad = Matrix6::new([[z; 6]; 6]);
    for i in 0..3 {
        for j in 0..3 {
            ad.data[i][j] = r.data[i][j];
            ad.data[3 + i][3 + j] = r.data[i][j];
            ad.data[3 + i][j] = tl.data[i][j];
        }
    }
    ad
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn spatial_inertia_unit_box() {
        let mass = 1.0f64;
        let com = Vector3::new([0.0f64, 0.0, 0.0]);
        let si = SpatialInertia::new(
            mass,
            com,
            Matrix3::new([
                [1.0 / 6.0, 0.0, 0.0],
                [0.0, 1.0 / 6.0, 0.0],
                [0.0, 0.0, 1.0 / 6.0],
            ]),
        );
        assert!((si.mass - 1.0f64).abs() < 1e-12f64);
        for i in 0..3 {
            assert!((si.matrix.data[i][i] - 1.0f64).abs() < 1e-12f64);
        }
    }

    #[test]
    fn spatial_inertia_momentum() {
        let si = SpatialInertia::new(
            2.0f64,
            Vector3::new([0.0f64; 3]),
            Matrix3::new([[2.0f64; 3]; 3]),
        );
        let v = SpatialVelocity::new(Vector3::new([1.0f64, 0.0, 0.0]), Vector3::new([0.0f64; 3]));
        let h = si.momentum(&v);
        assert!((h.angular().x() - 2.0f64).abs() < 1e-12f64);
        assert!(h.linear().norm() < 1e-12f64);
    }

    #[test]
    fn spatial_inertia_compose() {
        let i1 = SpatialInertia::new(
            1.0f64,
            Vector3::new([0.0f64; 3]),
            Matrix3::new([[1.0f64; 3]; 3]),
        );
        let i2 = SpatialInertia::new(
            2.0f64,
            Vector3::new([1.0f64, 0.0, 0.0]),
            Matrix3::new([[1.0f64; 3]; 3]),
        );
        assert!((i1.compose(&i2).mass - 3.0f64).abs() < 1e-12f64);
    }

    #[test]
    fn motion_cross_identity() {
        let v = SpatialVelocity::new(
            Vector3::new([1.0f64, 0.0, 0.0]),
            Vector3::new([0.0f64, 1.0, 0.0]),
        );
        let zero = SpatialVelocity::new(Vector3::new([0.0f64; 3]), Vector3::new([0.0f64; 3]));
        let r = motion_cross(&v, &zero);
        assert!(r.angular().norm() < 1e-12f64);
        assert!(r.linear().norm() < 1e-12f64);
    }

    #[test]
    fn force_cross_identity() {
        let f = SpatialForce::new(
            Vector3::new([0.0f64, 0.0, 1.0]),
            Vector3::new([1.0f64, 0.0, 0.0]),
        );
        let g = SpatialForce::new(Vector3::new([0.0f64; 3]), Vector3::new([1.0f64, 0.0, 0.0]));
        let r = force_cross(&f, &g);
        assert!((r.angular().y() - 1.0f64).abs() < 1e-12f64);
    }

    #[test]
    fn rigid_body_construction() {
        let si = SpatialInertia::new(
            1.0f64,
            Vector3::new([0.0f64; 3]),
            Matrix3::new([[1.0f64; 3]; 3]),
        );
        let body = RigidBody::new(si, Isometry3::identity(), "link0", 0);
        assert_eq!(body.name, "link0");
        assert_eq!(body.collision_geometry, 0);
    }

    #[test]
    fn skew_matrix_properties() {
        let v = Vector3::new([1.0f64, 2.0, 3.0]);
        let s = skew_matrix(&v);
        for i in 0..3 {
            for j in 0..3 {
                assert!((s.data[i][j] + s.data[j][i]).abs() < 1e-12f64);
            }
        }
        let sv = s * v;
        assert!(sv.norm() < 1e-12f64);
    }

    #[test]
    fn spatial_inertia_inertia_about() {
        let mass = 1.0f64;
        let l = 2.0f64;
        // Rod along x from -l/2 to l/2: I_com = diag(0, ML²/12, ML²/12)
        let si = SpatialInertia::new(
            mass,
            Vector3::new([0.0f64; 3]),
            Matrix3::new([
                [0.0f64, 0.0, 0.0],
                [0.0, mass * l * l / 12.0, 0.0],
                [0.0, 0.0, mass * l * l / 12.0],
            ]),
        );
        // Transfer to end: I_end_zz = I_com_zz + M·(L/2)² = ML²/12 + ML²/4 = ML²/3
        let si_end = si.inertia_about(&Vector3::new([l / 2.0f64, 0.0, 0.0]));
        let expected = mass * l * l / 3.0f64;
        assert!((si_end.inertia_com.data[2][2] - expected).abs() < 1e-10f64);
    }
}
