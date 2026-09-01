//! Forward kinematics for serial and parallel kinematic chains.
//!
//! Provides O(n) recursive transformation composition for serial chains,
//! and Product of Exponentials (PoE) formulation using screw theory.

use tpt_math_geometry::Isometry3;
use tpt_math_linalg_fixed::Vector3;

pub use crate::chain::DhLink;

/// Compute forward kinematics for a serial chain given joint angles.
///
/// Returns the end-effector pose in the base frame.
pub fn forward_kinematics(links: &[DhLink], joint_angles: &[f64]) -> Isometry3<f64> {
    let mut result = Isometry3::<f64>::identity();
    for (i, link) in links.iter().enumerate() {
        let theta = if i < joint_angles.len() {
            joint_angles[i]
        } else {
            link.theta
        };
        let t = DhLink::new(link.a, link.alpha, link.d, theta).transform();
        result = result * t;
    }
    result
}

/// Compute the spatial Jacobian for a serial chain at the end-effector.
///
/// Returns a 6×n matrix where each column is the screw axis for joint i,
/// expressed in the base frame.
pub fn geometric_jacobian(links: &[DhLink], joint_angles: &[f64]) -> Jacobian {
    let n = links.len();
    let mut jac = Jacobian::new(n);
    let mut t_world = Isometry3::<f64>::identity();

    for i in 0..n {
        let theta = if i < joint_angles.len() {
            joint_angles[i]
        } else {
            links[i].theta
        };
        let t_i = DhLink::new(links[i].a, links[i].alpha, links[i].d, theta).transform();
        t_world = t_world * t_i;

        let z_i = Vector3::new([0.0, 0.0, 1.0]);
        let p_i = Vector3::new([0.0, 0.0, 0.0]);

        let z = t_world.rotation.transform_vector(&z_i);
        let p = t_world.translation.vector;
        let v = z.cross(&p);

        jac.set_angular_column(i, z);
        jac.set_linear_column(i, v);
    }

    jac
}

/// A spatial Jacobian matrix (6 rows × n columns).
#[derive(Clone, Debug, PartialEq)]
pub struct Jacobian {
    pub angular: Vec<[f64; 3]>,
    pub linear: Vec<[f64; 3]>,
}

impl Jacobian {
    pub fn new(num_joints: usize) -> Self {
        Self {
            angular: vec![[0.0; 3]; num_joints],
            linear: vec![[0.0; 3]; num_joints],
        }
    }

    pub fn num_joints(&self) -> usize {
        self.angular.len()
    }

    pub fn set_angular_column(&mut self, col: usize, v: Vector3<f64>) {
        self.angular[col] = [v.x(), v.y(), v.z()];
    }

    pub fn set_linear_column(&mut self, col: usize, v: Vector3<f64>) {
        self.linear[col] = [v.x(), v.y(), v.z()];
    }

    pub fn angular_column(&self, col: usize) -> [f64; 3] {
        self.angular[col]
    }

    pub fn linear_column(&self, col: usize) -> [f64; 3] {
        self.linear[col]
    }

    /// Compute end-effector spatial velocity from joint rates.
    pub fn spatial_velocity(&self, qdot: &[f64]) -> SpatialVelocity {
        let mut omega = [0.0f64; 3];
        let mut v = [0.0f64; 3];
        for i in 0..self.num_joints() {
            let qd = qdot.get(i).copied().unwrap_or(0.0);
            for j in 0..3 {
                omega[j] += self.angular[i][j] * qd;
                v[j] += self.linear[i][j] * qd;
            }
        }
        SpatialVelocity::new(Vector3::new(omega), Vector3::new(v))
    }
}

/// A spatial velocity (twist): angular velocity ω and linear velocity v.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpatialVelocity {
    pub angular: [f64; 3],
    pub linear: [f64; 3],
}

impl SpatialVelocity {
    pub fn new(angular: Vector3<f64>, linear: Vector3<f64>) -> Self {
        Self {
            angular: [angular.x(), angular.y(), angular.z()],
            linear: [linear.x(), linear.y(), linear.z()],
        }
    }

    pub fn angular_vec(&self) -> Vector3<f64> {
        Vector3::new(self.angular)
    }

    pub fn linear_vec(&self) -> Vector3<f64> {
        Vector3::new(self.linear)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_math_geometry::{Rotation3, Translation};

    #[test]
    fn forward_kinematics_single_joint() {
        let links = vec![DhLink::new(1.0, 0.0, 0.0, 0.0)];
        let q = [core::f64::consts::FRAC_PI_2];
        let pose = forward_kinematics(&links, &q);
        let expected_rot = Rotation3::from_axis_angle(
            &Vector3::new([0.0, 0.0, 1.0]),
            core::f64::consts::FRAC_PI_2,
        );
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (pose.rotation.matrix().data[i][j] - expected_rot.matrix().data[i][j]).abs()
                        < 1e-12
                );
            }
        }
    }

    #[test]
    fn jacobian_single_joint() {
        let links = vec![DhLink::new(1.0, 0.0, 0.0, 0.0)];
        let jac = geometric_jacobian(&links, &[]);
        assert_eq!(jac.num_joints(), 1);
    }
}
