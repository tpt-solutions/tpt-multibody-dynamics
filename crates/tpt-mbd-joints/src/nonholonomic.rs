//! Nonholonomic constraints for multibody systems.
//!
//! Provides:
//! - [`NonholonomicConstraint`] trait for velocity-level constraints
//! - [`RollingWithoutSlipping`] — planar rolling constraint (v = ω × r)
//! - [`GearConstraint`] — linear velocity relationship between joints

#![allow(missing_docs)]

use alloc::vec;
use core::fmt;
use num_traits::Float;
use tpt_math_linalg_fixed::Vector3;

use crate::joint::JointAxis;

/// A nonholonomic (velocity-level) constraint.
///
/// Unlike holonomic constraints which are expressed as Φ(q) = 0,
/// nonholonomic constraints are linear in velocities: `g(q, q̇) = 0`.
pub trait NonholonomicConstraint: fmt::Debug {
    /// Evaluate the velocity-level constraint `g(q, q̇)`.
    fn velocity_constraint(&self, q: &[f64], qdot: &[f64]) -> alloc::vec::Vec<f64>;

    /// Evaluate the constraint Jacobian wrt velocities `∂g/∂q̇`.
    ///
    /// Each row `[f64; 3]` represents one constraint equation's partial
    /// derivatives wrt the 3 velocity DOFs of the associated body.
    fn jacobian_qdot(&self, q: &[f64], qdot: &[f64]) -> alloc::vec::Vec<[f64; 3]>;

    /// The number of constraint equations.
    fn num_constraints(&self) -> usize;

    /// Compute the constraint violation norm ||g||₂.
    fn violation(&self, q: &[f64], qdot: &[f64]) -> f64 {
        let g = self.velocity_constraint(q, qdot);
        let sum: f64 = g.iter().map(|x| x * x).sum();
        Float::sqrt(sum)
    }
}

#[allow(dead_code)]
fn body_pos(q: &[f64], body: usize) -> [f64; 3] {
    [
        q.get(body * 6).copied().unwrap_or(0.0),
        q.get(body * 6 + 1).copied().unwrap_or(0.0),
        q.get(body * 6 + 2).copied().unwrap_or(0.0),
    ]
}

fn body_vel(qdot: &[f64], body: usize) -> [f64; 3] {
    [
        qdot.get(body * 6).copied().unwrap_or(0.0),
        qdot.get(body * 6 + 1).copied().unwrap_or(0.0),
        qdot.get(body * 6 + 2).copied().unwrap_or(0.0),
    ]
}

fn body_omega(qdot: &[f64], body: usize) -> [f64; 3] {
    [
        qdot.get(body * 6 + 3).copied().unwrap_or(0.0),
        qdot.get(body * 6 + 4).copied().unwrap_or(0.0),
        qdot.get(body * 6 + 5).copied().unwrap_or(0.0),
    ]
}

/// Rolling without slipping constraint.
///
/// Enforces that the contact point velocity is zero:
/// `v_contact = v_body + ω × r = 0`
///
/// This produces 2 in-plane velocity constraints (normal direction is
/// handled by the contact normal constraint separately).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RollingWithoutSlipping {
    pub body: usize,
    pub contact_point: [f64; 3],
    pub normal: JointAxis,
}

impl RollingWithoutSlipping {
    /// Create a new rolling without slipping constraint.
    pub fn new(body: usize, contact_point: [f64; 3], normal: JointAxis) -> Self {
        Self {
            body,
            contact_point,
            normal,
        }
    }
}

impl NonholonomicConstraint for RollingWithoutSlipping {
    fn velocity_constraint(&self, _q: &[f64], qdot: &[f64]) -> alloc::vec::Vec<f64> {
        let v = body_vel(qdot, self.body);
        let omega = body_omega(qdot, self.body);
        let r = Vector3::new(self.contact_point);
        let w_cross_r = Vector3::new([
            omega[1] * r.z() - omega[2] * r.y(),
            omega[2] * r.x() - omega[0] * r.z(),
            omega[0] * r.y() - omega[1] * r.x(),
        ]);
        let v_contact = [
            v[0] + w_cross_r.x(),
            v[1] + w_cross_r.y(),
            v[2] + w_cross_r.z(),
        ];
        match self.normal {
            JointAxis::X => vec![v_contact[1], v_contact[2]],
            JointAxis::Y => vec![v_contact[0], v_contact[2]],
            JointAxis::Z => vec![v_contact[0], v_contact[1]],
        }
    }

    fn jacobian_qdot(&self, _q: &[f64], _qdot: &[f64]) -> alloc::vec::Vec<[f64; 3]> {
        match self.normal {
            JointAxis::X => vec![[0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            JointAxis::Y => vec![[1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
            JointAxis::Z => vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        }
    }

    fn num_constraints(&self) -> usize {
        2
    }
}

/// Gear constraint: linear velocity relationship between two joints.
///
/// Enforces `q̇_i = r * q̇_j` where `r` is the gear ratio.
///
/// This is a nonholonomic constraint at the velocity level.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GearConstraint {
    pub joint_i: usize,
    pub joint_j: usize,
    pub ratio: f64,
}

impl GearConstraint {
    /// Create a new gear constraint.
    pub fn new(joint_i: usize, joint_j: usize, ratio: f64) -> Self {
        Self {
            joint_i,
            joint_j,
            ratio,
        }
    }
}

impl NonholonomicConstraint for GearConstraint {
    fn velocity_constraint(&self, q: &[f64], qdot: &[f64]) -> alloc::vec::Vec<f64> {
        let _ = q;
        let qdot_i = qdot.get(self.joint_i).copied().unwrap_or(0.0);
        let qdot_j = qdot.get(self.joint_j).copied().unwrap_or(0.0);
        vec![qdot_i - self.ratio * qdot_j]
    }

    fn jacobian_qdot(&self, _q: &[f64], _qdot: &[f64]) -> alloc::vec::Vec<[f64; 3]> {
        vec![[1.0, 0.0, 0.0]]
    }

    fn num_constraints(&self) -> usize {
        1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rolling_without_slipping() {
        let rolling = RollingWithoutSlipping::new(0, [0.0, 0.0, -0.5], JointAxis::Z);
        // Body at origin, no velocity, no angular velocity → rolling constraint satisfied
        let q = vec![0.0, 0.0, 0.5, 0.0, 0.0, 0.0];
        let qdot = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        assert!(rolling.violation(&q, &qdot) < 1e-10);
    }

    #[test]
    fn test_gear_constraint() {
        let gear = GearConstraint::new(0, 1, 2.0);
        let q = vec![0.0; 6];
        // qdot_i = 2.0 * qdot_j: set qdot_j = 1.0, qdot_i = 2.0
        let qdot = vec![2.0, 1.0, 0.0, 0.0, 0.0, 0.0];
        assert!(gear.violation(&q, &qdot) < 1e-10);
    }

    #[test]
    fn test_gear_constraint_violated() {
        let gear = GearConstraint::new(0, 1, 2.0);
        let q = vec![0.0; 6];
        let qdot = vec![2.0, 5.0, 0.0, 0.0, 0.0, 0.0];
        assert!(gear.violation(&q, &qdot) > 1e-6);
    }

    #[test]
    fn test_nonholonomic_jacobian_dimensions() {
        let rolling = RollingWithoutSlipping::new(0, [0.0, 0.0, -0.5], JointAxis::Z);
        let q = vec![0.0; 6];
        let qdot = vec![0.0; 6];
        let jac = rolling.jacobian_qdot(&q, &qdot);
        assert_eq!(jac.len(), rolling.num_constraints());
    }
}
