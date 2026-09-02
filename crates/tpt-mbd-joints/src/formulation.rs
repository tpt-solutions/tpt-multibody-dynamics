//! Coordinate formulation support for multibody systems.
//!
//! Provides:
//! - [`MinimalCoordinateFormulation`] — reduced-coordinate (joint-space) parameterization
//! - [`MaximalCoordinateFormulation`] — full 6-DOF-per-body with Lagrange multipliers
//! - Conversion utilities between formulations

#![allow(missing_docs)]

use alloc::boxed::Box;
use alloc::vec;
use alloc::vec::Vec;

use crate::constraint::JointConstraint;
use crate::joint::JointType;
use num_traits::Float;

/// A minimal-coordinate (reduced) formulation for a multibody system.
///
/// In the minimal-coordinate formulation, the system is parameterized by
/// independent joint coordinates only. Body poses are recovered via
/// forward kinematics. This formulation is only valid for tree-topology
/// systems.
///
/// For a system with `n` bodies and `m` joints, the number of generalized
/// coordinates is `6 + Σ joint.num_dofs` (base body has 6 DOFs, each joint
/// adds its DOFs).
#[derive(Clone, Debug, PartialEq)]
pub struct MinimalCoordinateFormulation {
    /// Joint types in tree order (body 1 .. body n).
    pub joints: Vec<JointType>,
    /// Number of independent generalized coordinates.
    pub num_dofs: usize,
}

impl MinimalCoordinateFormulation {
    /// Create a new minimal-coordinate formulation.
    pub fn new(joints: Vec<JointType>) -> Self {
        let num_dofs = 6 + joints.iter().map(|j| j.num_dofs).sum::<usize>();
        Self { joints, num_dofs }
    }

    /// Get the number of joints.
    pub fn num_joints(&self) -> usize {
        self.joints.len()
    }

    /// Get the number of bodies (base + one per joint).
    pub fn num_bodies(&self) -> usize {
        1 + self.joints.len()
    }

    /// Get the DOF count for a specific joint by index.
    pub fn joint_dofs(&self, joint_idx: usize) -> usize {
        self.joints.get(joint_idx).map(|j| j.num_dofs).unwrap_or(0)
    }

    /// Map joint-space coordinates to a full 6n state vector.
    ///
    /// `q_joint` contains [base_6_dofs, joint_1_dofs, ..., joint_m_dofs].
    /// Returns a flat vector of length `6 * num_bodies`.
    pub fn joint_to_full(&self, q_joint: &[f64]) -> Vec<f64> {
        let n = self.num_bodies();
        let mut q_full = vec![0.0; n * 6];

        if q_joint.len() >= 6 {
            q_full[0..6].copy_from_slice(&q_joint[0..6]);
        }

        let mut offset = 6;
        let mut jq_offset = 6;
        for joint in &self.joints {
            let dofs = joint.num_dofs;
            if jq_offset + dofs <= q_joint.len() {
                for i in 0..dofs {
                    q_full[offset + i] = q_joint[jq_offset + i];
                }
            }
            offset += 6;
            jq_offset += dofs;
        }

        q_full
    }
}

/// A maximal-coordinate (redundant) formulation for a multibody system.
///
/// In the maximal-coordinate formulation, each body retains its full 6 DOFs
/// and constraints are enforced via Lagrange multipliers. This formulation
/// is compatible with index-3 DAE solvers and handles general constraint
/// topologies including loops.
///
/// The equations of motion take the form:
///
/// ```text
/// M·q̈ + Φ_qᵀ·λ = τ
/// Φ(q) = 0
/// ```
///
/// where `M` is the mass matrix, `Φ_q` is the constraint Jacobian, `λ` are
/// the Lagrange multipliers, and `τ` is the generalized force vector.
#[derive(Debug)]
pub struct MaximalCoordinateFormulation {
    /// Total number of DOFs (6 per body).
    pub num_dofs: usize,
    /// Number of constraint equations.
    pub num_constraints: usize,
    /// Holonomic constraints.
    pub constraints: Vec<Box<dyn JointConstraint>>,
    /// Baumgarte stabilization parameters.
    pub baumgarte: crate::stabilization::BaumgarteParams,
}

impl MaximalCoordinateFormulation {
    /// Create a new maximal-coordinate formulation.
    pub fn new(num_bodies: usize, constraints: Vec<Box<dyn JointConstraint>>) -> Self {
        let num_dofs = num_bodies * 6;
        let num_constraints = constraints.iter().map(|c| c.num_constraints()).sum();
        Self {
            num_dofs,
            num_constraints,
            constraints,
            baumgarte: crate::stabilization::BaumgarteParams::default(),
        }
    }

    /// Get the number of bodies.
    pub fn num_bodies(&self) -> usize {
        self.num_dofs / 6
    }

    /// Evaluate all constraint equations Φ(q).
    pub fn constraint_vector(&self, q: &[f64]) -> Vec<f64> {
        let mut phi = Vec::new();
        for c in &self.constraints {
            phi.extend(c.constraint(q));
        }
        phi
    }

    /// Build the full constraint Jacobian Φ_q.
    ///
    /// Returns a matrix with `num_constraints` rows and `num_dofs` columns.
    pub fn constraint_jacobian(&self, q: &[f64]) -> Vec<Vec<f64>> {
        let mut jac = Vec::new();
        for c in &self.constraints {
            let rows = c.jacobian(q);
            let body_i = 0;
            for row in rows {
                let mut full_row = vec![0.0; self.num_dofs];
                for (j, &val) in row.iter().enumerate() {
                    let idx = body_i * 6 + j;
                    if idx < full_row.len() {
                        full_row[idx] = val;
                    }
                }
                jac.push(full_row);
            }
        }
        jac
    }

    /// Compute constraint violation norm ||Φ||₂.
    pub fn violation(&self, q: &[f64]) -> f64 {
        let phi = self.constraint_vector(q);
        let sum: f64 = phi.iter().map(|x| x * x).sum();
        Float::sqrt(sum)
    }
}

/// Formulation type selector.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Formulation {
    /// Reduced-coordinate (joint-space) formulation.
    Minimal,
    /// Full-coordinate (body-space) formulation with Lagrange multipliers.
    Maximal,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constraint::SphericalConstraint;
    use alloc::boxed::Box;

    #[test]
    fn test_minimal_formulation_creation() {
        let joints = vec![JointType::REVOLUTE, JointType::REVOLUTE];
        let form = MinimalCoordinateFormulation::new(joints);
        assert_eq!(form.num_joints(), 2);
        assert_eq!(form.num_bodies(), 3);
        assert_eq!(form.num_dofs, 8);
    }

    #[test]
    fn test_minimal_formulation_joint_dofs() {
        let joints = vec![JointType::REVOLUTE, JointType::SPHERICAL];
        let form = MinimalCoordinateFormulation::new(joints);
        assert_eq!(form.joint_dofs(0), 1);
        assert_eq!(form.joint_dofs(1), 3);
    }

    #[test]
    fn test_minimal_formulation_map_to_full() {
        let joints = vec![JointType::REVOLUTE];
        let form = MinimalCoordinateFormulation::new(joints);
        let q_joint = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        let q_full = form.joint_to_full(&q_joint);
        assert_eq!(q_full.len(), 12);
        assert_eq!(q_full[0], 1.0);
        assert_eq!(q_full[6], 7.0);
    }

    #[test]
    fn test_maximal_formulation_violation() {
        let c = SphericalConstraint::new(0, 1);
        let constraints: Vec<Box<dyn JointConstraint>> = vec![Box::new(c)];
        let form = MaximalCoordinateFormulation::new(2, constraints);
        let q = vec![0.0; 12];
        assert!(form.violation(&q) < 1e-10);
    }

    #[test]
    fn test_maximal_formulation_violation_nonzero() {
        let c = SphericalConstraint::new(0, 1);
        let constraints: Vec<Box<dyn JointConstraint>> = vec![Box::new(c)];
        let form = MaximalCoordinateFormulation::new(2, constraints);
        let q = vec![1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        assert!(form.violation(&q) > 1e-6);
    }

    #[test]
    fn test_maximal_formulation_constraint_vector() {
        let c = SphericalConstraint::new(0, 1);
        let constraints: Vec<Box<dyn JointConstraint>> = vec![Box::new(c)];
        let form = MaximalCoordinateFormulation::new(2, constraints);
        let q = vec![1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 0.0, 0.0, 0.0];
        let phi = form.constraint_vector(&q);
        assert_eq!(phi.len(), 3);
    }

    #[test]
    fn test_maximal_formulation_jacobian_shape() {
        let c = SphericalConstraint::new(0, 1);
        let constraints: Vec<Box<dyn JointConstraint>> = vec![Box::new(c)];
        let form = MaximalCoordinateFormulation::new(2, constraints);
        let q = vec![0.0; 12];
        let jac = form.constraint_jacobian(&q);
        assert_eq!(jac.len(), 3);
        assert_eq!(jac[0].len(), 12);
    }

    #[test]
    fn test_formulation_enum() {
        let min = Formulation::Minimal;
        let max = Formulation::Maximal;
        assert_ne!(min, max);
    }
}
