//! Product of Exponentials (PoE) forward kinematics.
//!
//! Implements the PoE formulation: `T = exp(ξ₁θ₁) ... exp(ξₙθₙ) M`
//!
//! Where `ξᵢ` are the screw axes for each joint and `M` is the home
//! configuration (end-effector pose when all joint variables are zero).

extern crate alloc;

use alloc::vec::Vec;

use tpt_math_geometry::{Isometry3, Rotation3, Translation};
use tpt_math_linalg_fixed::{Matrix3, Vector3};

use crate::chain::DhLink;
use crate::forward::forward_kinematics;

/// A screw axis for a single joint: `ξ = [ω; v]` where `ω` is the unit
/// angular velocity and `v` is the linear velocity of a point on the axis.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScrewAxis {
    pub angular: Vector3<f64>,
    pub linear: Vector3<f64>,
}

impl ScrewAxis {
    /// Create a new screw axis.
    pub fn new(angular: Vector3<f64>, linear: Vector3<f64>) -> Self {
        Self { angular, linear }
    }

    /// Create a revolute screw axis from axis direction and a point on the axis.
    pub fn revolute(axis: Vector3<f64>, point_on_axis: Vector3<f64>) -> Self {
        let omega = axis.normalize();
        let v = omega.cross(&point_on_axis) * -1.0;
        Self::new(omega, v)
    }

    /// Create a prismatic screw axis from the translation direction.
    pub fn prismatic(direction: Vector3<f64>) -> Self {
        let v = direction.normalize();
        Self::new(Vector3::new([0.0, 0.0, 0.0]), v)
    }
}

/// Compute the homogeneous transform for a twist `ξ` scaled by `theta`.
///
/// Uses Rodrigues' rotation formula for the rotational part and the
/// `G(θ)` matrix for the translation part.
pub fn twist_exponential(screw: &ScrewAxis, theta: f64) -> Isometry3<f64> {
    let omega = screw.angular;
    let v = screw.linear;

    let omega_norm = omega.norm();
    if omega_norm < 1e-12 {
        let t = Translation::new(v * theta);
        return Isometry3::new(t, Rotation3::identity());
    }

    let omega_hat = omega / omega_norm;
    let angle = theta * omega_norm;

    let wx = omega_hat.data[0];
    let wy = omega_hat.data[1];
    let wz = omega_hat.data[2];

    let c = angle.cos();
    let s = angle.sin();
    let one_minus_c = 1.0 - c;

    let r = Matrix3::new([
        [
            c + wx * wx * one_minus_c,
            wx * wy * one_minus_c - wz * s,
            wx * wz * one_minus_c + wy * s,
        ],
        [
            wy * wx * one_minus_c + wz * s,
            c + wy * wy * one_minus_c,
            wy * wz * one_minus_c - wx * s,
        ],
        [
            wz * wx * one_minus_c - wy * s,
            wz * wy * one_minus_c + wx * s,
            c + wz * wz * one_minus_c,
        ],
    ]);

    let vx = v.data[0];
    let vy = v.data[1];
    let vz = v.data[2];

    let g00 = angle - wx * wx * angle + wx * wx * angle * c + wx * wy * s - wx * wy * angle * s
        + wx * wz * s + wx * wz * angle * s;
    let g01 = wx * wy * angle - wx * wy * angle * c - wz * s + wz * angle * s
        + wy * wz * angle - wy * wz * angle * c - wx * s;
    let g02 = wx * wz * angle - wx * wz * angle * c + wy * s - wy * angle * s
        + wy * wz * angle - wy * wz * angle * c - wx * s;

    let g10 = wy * wx * angle - wy * wx * angle * c + wz * s - wz * angle * s
        + wx * wy * angle - wx * wy * angle * c - wz * s;
    let g11 = angle - wy * wy * angle + wy * wy * angle * c + wz * wx * s + wz * wx * angle * s
        + wy * wz * angle - wy * wz * angle * c + wx * s;
    let g12 = wy * wz * angle - wy * wz * angle * c - wx * s + wx * angle * s
        + wz * wx * angle - wz * wx * angle * c - wy * s;

    let g20 = wz * wx * angle - wz * wx * angle * c - wy * s + wy * angle * s
        + wx * wz * angle - wx * wz * angle * c + wy * s;
    let g21 = wz * wy * angle - wz * wy * angle * c + wx * s - wx * angle * s
        + wy * wz * angle - wy * wz * angle * c + wx * s;
    let g22 = angle - wz * wz * angle + wz * wz * angle * c + wx * wy * s - wx * wy * angle * s
        + wy * wx * s - wy * wx * angle * s;

    let t0 = (vx * g00 + vy * g01 + vz * g02) / omega_norm;
    let t1 = (vx * g10 + vy * g11 + vz * g12) / omega_norm;
    let t2 = (vx * g20 + vy * g21 + vz * g22) / omega_norm;

    let translation = Translation::new(Vector3::new([t0, t1, t2]));
    let rotation = Rotation3::from_matrix_unchecked(r);
    Isometry3::new(translation, rotation)
}

/// Compute forward kinematics using the Product of Exponentials formulation.
///
/// Given screw axes for each joint and the home configuration `M`, computes:
/// `T = exp(ξ₁θ₁) ... exp(ξₙθₙ) M`
///
/// # Arguments
///
/// * `screws` — screw axis for each joint
/// * `thetas` — joint variables (angles for revolute, displacements for prismatic)
/// * `home` — end-effector pose when all thetas = 0
///
/// # Returns
///
/// The end-effector pose in the base frame.
pub fn poe_forward_kinematics(screws: &[ScrewAxis], thetas: &[f64], home: Isometry3<f64>) -> Isometry3<f64> {
    let mut t = Isometry3::identity();
    for (i, screw) in screws.iter().enumerate() {
        let theta = thetas.get(i).copied().unwrap_or(0.0);
        let exp_xi = twist_exponential(screw, theta);
        t = t * exp_xi;
    }
    t * home
}

/// Compute screw axes for a serial chain from DH parameters.
///
/// Returns one screw axis per joint. For standard DH parameters, the z-axis
/// of frame `i-1` is the revolute axis, and the origin of frame `i-1` is
/// a point on that axis.
pub fn dh_to_screw_axes(links: &[DhLink]) -> Vec<ScrewAxis> {
    let mut axes = Vec::new();
    let mut z_prev = Vector3::new([0.0, 0.0, 1.0]);
    let mut o_prev = Vector3::new([0.0, 0.0, 0.0]);

    for link in links {
        let theta = link.theta;
        let d = link.d;
        let a = link.a;
        let alpha = link.alpha;

        let rz = Rotation3::from_axis_angle(&Vector3::new([0.0, 0.0, 1.0]), theta);
        let tz = Translation::new(Vector3::new([0.0, 0.0, d]));
        let tx = Translation::new(Vector3::new([a, 0.0, 0.0]));
        let rx = Rotation3::from_axis_angle(&Vector3::new([1.0, 0.0, 0.0]), alpha);

        let t_i = Isometry3::new(tz, rz) * Isometry3::new(tx, rx);

        let z_i = t_i.rotation.transform_vector(&z_prev);
        let o_i = t_i.translation.vector + t_i.rotation.transform_vector(&o_prev);

        let screw = ScrewAxis::revolute(z_i, o_i);
        axes.push(screw);

        z_prev = z_i;
        o_prev = o_i;
    }

    axes
}

/// Compute the home configuration `M` from DH parameters.
///
/// `M` is the end-effector pose when all joint variables are zero.
pub fn dh_home_configuration(links: &[DhLink]) -> Isometry3<f64> {
    forward_kinematics(links, &[])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_screw_axis_revolute() {
        let axis = Vector3::new([0.0, 0.0, 1.0]);
        let point = Vector3::new([1.0, 0.0, 0.0]);
        let screw = ScrewAxis::revolute(axis, point);
        assert!((screw.angular.norm() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_screw_axis_prismatic() {
        let dir = Vector3::new([0.0, 0.0, 1.0]);
        let screw = ScrewAxis::prismatic(dir);
        assert!((screw.linear.norm() - 1.0).abs() < 1e-10);
        assert_eq!(screw.angular, Vector3::new([0.0, 0.0, 0.0]));
    }

    #[test]
    fn test_twist_exponential_rotation_about_z() {
        let screw = ScrewAxis::revolute(Vector3::new([0.0, 0.0, 1.0]), Vector3::new([0.0, 0.0, 0.0]));
        let t = twist_exponential(&screw, core::f64::consts::PI);
        let rot = t.rotation.matrix();
        let expected = Matrix3::new([
            [-1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]);
        for i in 0..3 {
            for j in 0..3 {
                assert!((rot.data[i][j] - expected.data[i][j]).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn test_poe_matches_dh_single_joint() {
        let links = vec![DhLink::new(0.0, 0.0, 0.0, 0.0)];
        let screws = dh_to_screw_axes(&links);
        let home = dh_home_configuration(&links);
        let q = [core::f64::consts::FRAC_PI_2];
        let poe_pose = poe_forward_kinematics(&screws, &q, home);
        let dh_pose = forward_kinematics(&links, &q);
        assert!((poe_pose.translation.vector.data[0] - dh_pose.translation.vector.data[0]).abs() < 1e-9);
        assert!((poe_pose.translation.vector.data[1] - dh_pose.translation.vector.data[1]).abs() < 1e-9);
        assert!((poe_pose.translation.vector.data[2] - dh_pose.translation.vector.data[2]).abs() < 1e-9);
    }

    #[test]
    fn test_poe_matches_dh_two_joints() {
        let links = vec![DhLink::new(1.0, 0.0, 0.0, 0.0), DhLink::new(0.0, 0.0, 0.0, 0.0)];
        let screws = dh_to_screw_axes(&links);
        let home = dh_home_configuration(&links);
        let q = [0.5, -0.3];
        let poe_pose = poe_forward_kinematics(&screws, &q, home);
        let dh_pose = forward_kinematics(&links, &q);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (poe_pose.rotation.matrix().data[i][j] - dh_pose.rotation.matrix().data[i][j]).abs()
                        < 1e-9
                );
            }
        }
    }
}
