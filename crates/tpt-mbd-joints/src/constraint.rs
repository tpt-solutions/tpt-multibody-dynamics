//! Constraint formulation for multibody joints.
//!
//! Provides:
//! - [`JointConstraint`] trait for holonomic constraint equations Φ(q) = 0
//! - Constraint Jacobian Φ_q = ∂Φ/∂q
//! - Constraint violation metrics
//!
//! Implemented joint constraints:
//! - [`RevoluteConstraint`] — 1-DOF rotation about an axis (5 constraints)
//! - [`SphericalConstraint`] — ball-and-socket (3 position constraints)
//! - [`PrismaticConstraint`] — translation along a single axis (2 position constraints)
//! - [`UniversalConstraint`] — two orthogonal revolute axes (3 position constraints)
//! - [`CylindricalConstraint`] — revolute + prismatic along same axis (3 position constraints)
//! - [`PlanarConstraint`] — motion constrained to a plane (3 constraints)
//! - [`FixedConstraint`] — fully fixed relative pose (6 constraints)

#![allow(missing_docs)]

use alloc::vec;
use core::fmt;
use num_traits::Float;

use crate::joint::{JointAxis, JointType};

/// A holonomic constraint equation Φ(q) = 0 for a joint.
pub trait JointConstraint: fmt::Debug {
    /// Evaluate the constraint equation Φ(q).
    fn constraint(&self, q: &[f64]) -> alloc::vec::Vec<f64>;

    /// Evaluate the constraint Jacobian Φ_q = ∂Φ/∂q.
    ///
    /// Each row `[f64; 3]` represents the partial derivatives of one constraint
    /// equation with respect to the body's 3 position DOFs `[x, y, z]`.
    fn jacobian(&self, q: &[f64]) -> alloc::vec::Vec<[f64; 3]>;

    /// The number of constraint equations.
    fn num_constraints(&self) -> usize;

    /// The joint type.
    fn joint_type(&self) -> JointType;

    /// Compute the constraint violation norm ||Φ||₂.
    fn violation(&self, q: &[f64]) -> f64 {
        let phi = self.constraint(q);
        let sum: f64 = phi.iter().map(|x| x * x).sum();
        Float::sqrt(sum)
    }
}

fn body_pos(q: &[f64], body: usize) -> [f64; 3] {
    [
        q.get(body * 6).copied().unwrap_or(0.0),
        q.get(body * 6 + 1).copied().unwrap_or(0.0),
        q.get(body * 6 + 2).copied().unwrap_or(0.0),
    ]
}

fn body_rot(q: &[f64], body: usize) -> [f64; 3] {
    [
        q.get(body * 6 + 3).copied().unwrap_or(0.0),
        q.get(body * 6 + 4).copied().unwrap_or(0.0),
        q.get(body * 6 + 5).copied().unwrap_or(0.0),
    ]
}

/// Revolute joint constraint: 1-DOF rotation about an axis.
///
/// Enforces coincident joint origins with a fixed offset along the axis.
/// In the simplified model, 3 position constraints are enforced:
/// `Φ = [x_i - x_j + offset_x, offset_y, offset_z]`
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RevoluteConstraint {
    pub body_i: usize,
    pub body_j: usize,
    pub axis: JointAxis,
    pub offset: [f64; 3],
}

impl RevoluteConstraint {
    /// Create a new revolute constraint.
    pub fn new(body_i: usize, body_j: usize, axis: JointAxis) -> Self {
        Self {
            body_i,
            body_j,
            axis,
            offset: [0.0; 3],
        }
    }

    /// Set the joint axis offset.
    pub fn with_offset(mut self, offset: [f64; 3]) -> Self {
        self.offset = offset;
        self
    }
}

impl JointConstraint for RevoluteConstraint {
    fn constraint(&self, q: &[f64]) -> alloc::vec::Vec<f64> {
        let pi = body_pos(q, self.body_i);
        let pj = body_pos(q, self.body_j);
        vec![
            pi[0] - pj[0] + self.offset[0],
            pi[1] - pj[1] + self.offset[1],
            pi[2] - pj[2] + self.offset[2],
        ]
    }

    fn jacobian(&self, _q: &[f64]) -> alloc::vec::Vec<[f64; 3]> {
        let si = match self.axis {
            JointAxis::X => [1.0, 0.0, 0.0],
            JointAxis::Y => [0.0, 1.0, 0.0],
            JointAxis::Z => [0.0, 0.0, 1.0],
        };
        vec![[si[0], 0.0, 0.0], [0.0, si[1], 0.0], [0.0, 0.0, si[2]]]
    }

    fn num_constraints(&self) -> usize {
        3
    }

    fn joint_type(&self) -> JointType {
        JointType::REVOLUTE
    }
}

/// Spherical (ball-and-socket) joint constraint: 3 position constraints.
///
/// Enforces coincident origins: `Φ = p_i - p_j`
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SphericalConstraint {
    pub body_i: usize,
    pub body_j: usize,
}

impl SphericalConstraint {
    /// Create a new spherical constraint.
    pub fn new(body_i: usize, body_j: usize) -> Self {
        Self { body_i, body_j }
    }
}

impl JointConstraint for SphericalConstraint {
    fn constraint(&self, q: &[f64]) -> alloc::vec::Vec<f64> {
        let pi = body_pos(q, self.body_i);
        let pj = body_pos(q, self.body_j);
        vec![pi[0] - pj[0], pi[1] - pj[1], pi[2] - pj[2]]
    }

    fn jacobian(&self, _q: &[f64]) -> alloc::vec::Vec<[f64; 3]> {
        vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    }

    fn num_constraints(&self) -> usize {
        3
    }

    fn joint_type(&self) -> JointType {
        JointType::SPHERICAL
    }
}

/// Prismatic joint constraint: translation along a single axis.
///
/// Prevents relative translation perpendicular to the specified axis.
/// The relative translation along the axis is the free DOF.
///
/// Constraints: `Φ = [perp₁, perp₂]` (2 position constraints)
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PrismaticConstraint {
    pub body_i: usize,
    pub body_j: usize,
    pub axis: JointAxis,
    pub offset: [f64; 3],
}

impl PrismaticConstraint {
    /// Create a new prismatic constraint.
    pub fn new(body_i: usize, body_j: usize, axis: JointAxis) -> Self {
        Self {
            body_i,
            body_j,
            axis,
            offset: [0.0; 3],
        }
    }

    /// Set the joint axis offset.
    pub fn with_offset(mut self, offset: [f64; 3]) -> Self {
        self.offset = offset;
        self
    }
}

impl JointConstraint for PrismaticConstraint {
    fn constraint(&self, q: &[f64]) -> alloc::vec::Vec<f64> {
        let pi = body_pos(q, self.body_i);
        let pj = body_pos(q, self.body_j);
        match self.axis {
            JointAxis::X => vec![
                pi[1] - pj[1] + self.offset[1],
                pi[2] - pj[2] + self.offset[2],
            ],
            JointAxis::Y => vec![
                pi[0] - pj[0] + self.offset[0],
                pi[2] - pj[2] + self.offset[2],
            ],
            JointAxis::Z => vec![
                pi[0] - pj[0] + self.offset[0],
                pi[1] - pj[1] + self.offset[1],
            ],
        }
    }

    fn jacobian(&self, _q: &[f64]) -> alloc::vec::Vec<[f64; 3]> {
        match self.axis {
            JointAxis::X => vec![[0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            JointAxis::Y => vec![[1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
            JointAxis::Z => vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        }
    }

    fn num_constraints(&self) -> usize {
        2
    }

    fn joint_type(&self) -> JointType {
        JointType::PRISMATIC
    }
}

/// Universal (Cardan) joint constraint: two orthogonal revolute axes.
///
/// Enforces coincident origins. In the simplified model, 3 position
/// constraints are enforced.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UniversalConstraint {
    pub body_i: usize,
    pub body_j: usize,
    pub axis_1: JointAxis,
    pub axis_2: JointAxis,
}

impl UniversalConstraint {
    /// Create a new universal constraint.
    pub fn new(body_i: usize, body_j: usize, axis_1: JointAxis, axis_2: JointAxis) -> Self {
        Self {
            body_i,
            body_j,
            axis_1,
            axis_2,
        }
    }
}

impl JointConstraint for UniversalConstraint {
    fn constraint(&self, q: &[f64]) -> alloc::vec::Vec<f64> {
        let pi = body_pos(q, self.body_i);
        let pj = body_pos(q, self.body_j);
        vec![pi[0] - pj[0], pi[1] - pj[1], pi[2] - pj[2]]
    }

    fn jacobian(&self, _q: &[f64]) -> alloc::vec::Vec<[f64; 3]> {
        vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    }

    fn num_constraints(&self) -> usize {
        3
    }

    fn joint_type(&self) -> JointType {
        JointType::UNIVERSAL
    }
}

/// Cylindrical joint constraint: one revolute + one prismatic along the same axis.
///
/// Prevents relative translation perpendicular to the axis and all relative
/// rotation. The relative translation along the axis and rotation about it
/// are the free DOFs.
///
/// Constraints: `Φ = [perp₁, perp₂]` (2 position constraints)
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CylindricalConstraint {
    pub body_i: usize,
    pub body_j: usize,
    pub axis: JointAxis,
}

impl CylindricalConstraint {
    /// Create a new cylindrical constraint.
    pub fn new(body_i: usize, body_j: usize, axis: JointAxis) -> Self {
        Self {
            body_i,
            body_j,
            axis,
        }
    }
}

impl JointConstraint for CylindricalConstraint {
    fn constraint(&self, q: &[f64]) -> alloc::vec::Vec<f64> {
        let pi = body_pos(q, self.body_i);
        let pj = body_pos(q, self.body_j);
        match self.axis {
            JointAxis::X => vec![pi[1] - pj[1], pi[2] - pj[2]],
            JointAxis::Y => vec![pi[0] - pj[0], pi[2] - pj[2]],
            JointAxis::Z => vec![pi[0] - pj[0], pi[1] - pj[1]],
        }
    }

    fn jacobian(&self, _q: &[f64]) -> alloc::vec::Vec<[f64; 3]> {
        match self.axis {
            JointAxis::X => vec![[0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            JointAxis::Y => vec![[1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
            JointAxis::Z => vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        }
    }

    fn num_constraints(&self) -> usize {
        2
    }

    fn joint_type(&self) -> JointType {
        JointType::CYLINDRICAL
    }
}

/// Planar joint constraint: motion constrained to a plane.
///
/// Enforces that the two bodies share the same origin in the plane normal
/// direction. In-plane translation and rotation about the normal are free.
///
/// Constraints: `Φ = [n · (p_i - p_j)]` (1 position constraint in simplified model)
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanarConstraint {
    pub body_i: usize,
    pub body_j: usize,
    pub normal: JointAxis,
}

impl PlanarConstraint {
    /// Create a new planar constraint.
    pub fn new(body_i: usize, body_j: usize, normal: JointAxis) -> Self {
        Self {
            body_i,
            body_j,
            normal,
        }
    }
}

impl JointConstraint for PlanarConstraint {
    fn constraint(&self, q: &[f64]) -> alloc::vec::Vec<f64> {
        let pi = body_pos(q, self.body_i);
        let pj = body_pos(q, self.body_j);
        let diff = [pi[0] - pj[0], pi[1] - pj[1], pi[2] - pj[2]];
        match self.normal {
            JointAxis::X => vec![diff[0]],
            JointAxis::Y => vec![diff[1]],
            JointAxis::Z => vec![diff[2]],
        }
    }

    fn jacobian(&self, _q: &[f64]) -> alloc::vec::Vec<[f64; 3]> {
        match self.normal {
            JointAxis::X => vec![[1.0, 0.0, 0.0]],
            JointAxis::Y => vec![[0.0, 1.0, 0.0]],
            JointAxis::Z => vec![[0.0, 0.0, 1.0]],
        }
    }

    fn num_constraints(&self) -> usize {
        1
    }

    fn joint_type(&self) -> JointType {
        JointType::PLANAR
    }
}

/// Fixed (weld) joint constraint: fully fixes relative pose.
///
/// Enforces 6 constraints: 3 position + 3 orientation.
/// In the simplified model, 6 constraints are enforced:
/// `Φ = [p_i - p_j, ω_i - ω_j]`
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FixedConstraint {
    pub body_i: usize,
    pub body_j: usize,
}

impl FixedConstraint {
    /// Create a new fixed constraint.
    pub fn new(body_i: usize, body_j: usize) -> Self {
        Self { body_i, body_j }
    }
}

impl JointConstraint for FixedConstraint {
    fn constraint(&self, q: &[f64]) -> alloc::vec::Vec<f64> {
        let pi = body_pos(q, self.body_i);
        let pj = body_pos(q, self.body_j);
        let ri = body_rot(q, self.body_i);
        let rj = body_rot(q, self.body_j);
        vec![
            pi[0] - pj[0],
            pi[1] - pj[1],
            pi[2] - pj[2],
            ri[0] - rj[0],
            ri[1] - rj[1],
            ri[2] - rj[2],
        ]
    }

    fn jacobian(&self, _q: &[f64]) -> alloc::vec::Vec<[f64; 3]> {
        vec![
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ]
    }

    fn num_constraints(&self) -> usize {
        6
    }

    fn joint_type(&self) -> JointType {
        JointType::FIXED
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_spherical_constraint_satisfied() {
        let c = SphericalConstraint::new(0, 1);
        let q = vec![1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 0.0, 0.0, 0.0];
        assert!(c.violation(&q) < 1e-10);
    }

    #[test]
    fn test_spherical_constraint_violated() {
        let c = SphericalConstraint::new(0, 1);
        let q = vec![1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 4.0, 5.0, 6.0, 0.0, 0.0, 0.0];
        let v = c.violation(&q);
        let expected = ((3.0_f64).powi(2) * 3.0).sqrt();
        assert!((v - expected).abs() < 1e-10);
    }

    #[test]
    fn test_revolute_constraint_satisfied() {
        let c = RevoluteConstraint::new(0, 1, JointAxis::Z);
        let q = vec![1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 0.0, 0.0, 0.0];
        assert!(c.violation(&q) < 1e-10);
    }

    #[test]
    fn test_prismatic_constraint() {
        let c = PrismaticConstraint::new(0, 1, JointAxis::Z);
        let q = vec![1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 1.0, 2.0, 5.0, 0.0, 0.0, 0.0];
        // x and y are equal, z differs → 0 violation for prismatic
        assert!(c.violation(&q) < 1e-10);
    }

    #[test]
    fn test_fixed_constraint_satisfied() {
        let c = FixedConstraint::new(0, 1);
        let q = vec![1.0, 2.0, 3.0, 0.1, 0.2, 0.3, 1.0, 2.0, 3.0, 0.1, 0.2, 0.3];
        assert!(c.violation(&q) < 1e-10);
    }

    #[test]
    fn test_fixed_constraint_violated() {
        let c = FixedConstraint::new(0, 1);
        let q = vec![1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 4.0, 5.0, 6.0, 0.0, 0.0, 0.0];
        let v = c.violation(&q);
        let expected = ((3.0_f64).powi(2) * 3.0).sqrt();
        assert!((v - expected).abs() < 1e-10);
    }

    #[test]
    fn test_universal_constraint() {
        let c = UniversalConstraint::new(0, 1, JointAxis::X, JointAxis::Y);
        let q = vec![1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 0.0, 0.0, 0.0];
        assert!(c.violation(&q) < 1e-10);
    }

    #[test]
    fn test_cylindrical_constraint() {
        let c = CylindricalConstraint::new(0, 1, JointAxis::Z);
        let q = vec![1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 1.0, 2.0, 7.0, 0.0, 0.0, 0.0];
        assert!(c.violation(&q) < 1e-10);
    }

    #[test]
    fn test_planar_constraint() {
        let c = PlanarConstraint::new(0, 1, JointAxis::Z);
        let q = vec![1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 0.0, 0.0, 0.0];
        assert!(c.violation(&q) < 1e-10);
    }

    #[test]
    fn test_constraint_jacobian_dimensions() {
        let c = FixedConstraint::new(0, 1);
        let q = vec![0.0; 12];
        let jac = c.jacobian(&q);
        assert_eq!(jac.len(), c.num_constraints());
        for row in &jac {
            assert_eq!(row.len(), 3);
        }
    }
}
