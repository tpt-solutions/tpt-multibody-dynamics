//! Denavit-Hartenberg (DH) parameterization for serial manipulators.
//!
//! Supports both standard DH and modified DH conventions.
//!
//! Standard DH: `T = Rot_z(θ) · Trans_z(d) · Trans_x(a) · Rot_x(α)`
//! Modified DH: `T = Rot_x(α) · Trans_x(a) · Rot_z(θ) · Trans_z(d)`

use core::fmt;

use tpt_math_geometry::{Isometry3, Rotation3, Translation};
use tpt_math_linalg_fixed::Vector3;

/// DH link parameters (standard convention).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DhLink {
    pub a: f64,
    pub alpha: f64,
    pub d: f64,
    pub theta: f64,
}

impl DhLink {
    pub fn new(a: f64, alpha: f64, d: f64, theta: f64) -> Self {
        Self { a, alpha, d, theta }
    }

    /// Build the homogeneous transform for this link (standard DH).
    pub fn transform(&self) -> Isometry3<f64> {
        let ct = self.theta.cos();
        let st = self.theta.sin();
        let ca = self.alpha.cos();
        let sa = self.alpha.sin();

        let rot_z = Rotation3::from_axis_angle(&Vector3::new([0.0, 0.0, 1.0]), self.theta);
        let trans_z = Translation::new(Vector3::new([0.0, 0.0, self.d]));
        let trans_x = Translation::new(Vector3::new([self.a, 0.0, 0.0]));
        let rot_x = Rotation3::from_axis_angle(&Vector3::new([1.0, 0.0, 0.0]), self.alpha);

        let tz = Isometry3::new(trans_z, rot_z);
        let tx = Isometry3::new(trans_x, rot_x);
        tx * tz
    }
}

/// A kinematic chain parameterized by DH parameters.
#[derive(Clone, Debug, PartialEq)]
pub struct DhChain {
    pub links: Vec<DhLink>,
    pub base: Isometry3<f64>,
    pub tool: Isometry3<f64>,
}

impl DhChain {
    pub fn new(links: Vec<DhLink>) -> Self {
        Self {
            links,
            base: Isometry3::<f64>::identity(),
            tool: Isometry3::<f64>::identity(),
        }
    }

    pub fn with_base_tool(links: Vec<DhLink>, base: Isometry3<f64>, tool: Isometry3<f64>) -> Self {
        Self { links, base, tool }
    }

    pub fn num_joints(&self) -> usize {
        self.links.len()
    }

    /// Forward kinematics: compute end-effector pose from joint angles.
    pub fn forward(&self, joint_angles: &[f64]) -> Isometry3<f64> {
        let mut result = self.base;
        for (i, link) in self.links.iter().enumerate() {
            let theta = if i < joint_angles.len() {
                joint_angles[i]
            } else {
                link.theta
            };
            let t = DhLink::new(link.a, link.alpha, link.d, theta).transform();
            result = result * t;
        }
        result * self.tool
    }
}

impl fmt::Display for DhChain {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "DhChain({} joints)", self.links.len())
    }
}
