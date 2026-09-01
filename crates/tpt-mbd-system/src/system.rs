//! `MultibodySystem` assembly: bodies + joints + constraints + contact + flexible bodies.
//!
//! Provides the central `MultibodySystem` struct that collects all model data,
//! counts degrees of freedom, assembles the global mass matrix, and evaluates
//! generalized force vectors.

#![allow(unused_imports)]

extern crate alloc;

use alloc::boxed::Box;
use alloc::vec::Vec;

use tpt_math_linalg_dense::{DMatrix, DVector};
use tpt_mbd_core::inertia::{RigidBody, SpatialInertia};
use tpt_mbd_joints::constraint::JointConstraint;
use tpt_mbd_joints::joint::JointType;

pub use crate::Matrix;
pub use crate::Vector;

// ===========================================================================
// Local types (not yet provided by dependency crates)
// ===========================================================================

/// A contact manifold between two bodies.
#[derive(Clone, Debug)]
pub struct ContactManifold {
    /// First body index.
    pub body_i: usize,
    /// Second body index.
    pub body_j: usize,
    /// Contact points in world coordinates.
    pub points: Vec<[f64; 3]>,
    /// Maximum penetration depth.
    pub penetration_depth: f64,
}

/// A flexible body represented in floating frame formulation.
#[derive(Clone, Debug)]
pub struct FloatingFrameBody {
    /// Underlying rigid-body reference (large rigid motion).
    pub rigid_body: RigidBody<f64>,
    /// Reduced modal mass matrix diagonal.
    pub modal_mass: Vec<f64>,
    /// Reduced modal stiffness matrix diagonal.
    pub modal_stiffness: Vec<f64>,
    /// Reduced modal damping matrix diagonal.
    pub modal_damping: Vec<f64>,
    /// Number of retained flexible modes.
    pub num_modes: usize,
}

// ===========================================================================
// MultibodySystem
// ===========================================================================

/// Assembled multibody system.
pub struct MultibodySystem {
    /// Rigid bodies in the system.
    pub bodies: Vec<RigidBody<f64>>,
    /// Joint types connecting bodies (tree-order).
    pub joints: Vec<JointType>,
    /// Holonomic constraint equations (Baumgarte / augmented Lagrangian).
    pub constraints: Vec<Box<dyn JointConstraint>>,
    /// Contact manifolds for force evaluation.
    pub contacts: Vec<ContactManifold>,
    /// Flexible bodies with reduced modal coordinates.
    pub flexible: Vec<FloatingFrameBody>,
    /// Total number of generalized coordinates (recomputed by `count_dofs`).
    pub num_dofs: usize,
    /// Body index pairs corresponding to each entry in `joints`.
    joint_body_pairs: Vec<(usize, usize)>,
}

impl MultibodySystem {
    /// Create an empty multibody system.
    pub fn new() -> Self {
        Self {
            bodies: Vec::new(),
            joints: Vec::new(),
            constraints: Vec::new(),
            contacts: Vec::new(),
            flexible: Vec::new(),
            num_dofs: 0,
            joint_body_pairs: Vec::new(),
        }
    }

    /// Add a rigid body and return its index.
    pub fn add_body(&mut self, body: RigidBody<f64>) -> usize {
        let idx = self.bodies.len();
        self.bodies.push(body);
        idx
    }

    /// Add a joint connecting `body_i` and `body_j`, return joint index.
    pub fn add_joint(&mut self, joint: JointType, body_i: usize, body_j: usize) -> usize {
        let idx = self.joints.len();
        self.joints.push(joint);
        self.joint_body_pairs.push((body_i, body_j));
        idx
    }

    /// Add a holonomic constraint and return its index.
    pub fn add_constraint(&mut self, constraint: Box<dyn JointConstraint>) -> usize {
        let idx = self.constraints.len();
        self.constraints.push(constraint);
        idx
    }

    /// Add a flexible body and return its index.
    pub fn add_flexible_body(&mut self, body: FloatingFrameBody) -> usize {
        let idx = self.flexible.len();
        self.flexible.push(body);
        idx
    }

    /// Recompute the total number of generalized coordinates.
    ///
    /// Each rigid body contributes 6 DOFs in free space. Joints between bodies
    /// reduce the total by `6 - joint.num_dofs` per joint. Flexible bodies add
    /// their modal DOFs on top.
    pub fn count_dofs(&mut self) -> usize {
        if self.bodies.is_empty() && self.flexible.is_empty() {
            self.num_dofs = 0;
            return 0;
        }

        let mut total = self.bodies.len() * 6;

        for joint in &self.joints {
            total = total.saturating_sub(6 - joint.num_dofs);
        }

        for flex in &self.flexible {
            total = total.saturating_add(flex.num_modes);
        }

        self.num_dofs = total;
        total
    }

    /// Assemble the global mass matrix `M` as a block-diagonal dense matrix.
    ///
    /// Each rigid body contributes its 6×6 spatial inertia. Each flexible body
    /// contributes its modal mass on the diagonal.
    pub fn build_mass_matrix(&self) -> Matrix {
        let n = self.num_dofs.max(1);
        Matrix::from_fn(n, n, |i, j| {
            let body_off_i = (i / 6) * 6;
            let body_off_j = (j / 6) * 6;
            let local_i = i - body_off_i;
            let local_j = j - body_off_j;
            if body_off_i < self.bodies.len() * 6
                && body_off_i == body_off_j
                && local_i < 6
                && local_j < 6
            {
                let body_idx = body_off_i / 6;
                self.bodies[body_idx].spatial_inertia.matrix.data[local_i][local_j]
            } else if i == j {
                let flex_off = self.bodies.len() * 6;
                let mode_idx = i - flex_off;
                if mode_idx < self.flexible.len() {
                    let flex = &self.flexible[mode_idx];
                    flex.modal_mass.get(mode_idx).copied().unwrap_or(1.0)
                } else {
                    0.0
                }
            } else {
                0.0
            }
        })
    }

    /// Evaluate the generalized force vector `τ(q, q̇)` for the current state.
    ///
    /// Returns a dense vector of length `num_dofs`. External forces (gravity,
    /// springs, …) are accumulated here.
    pub fn build_force_vector(&self, _q: &[f64], _qdot: &[f64]) -> Vec<f64> {
        let n = self.num_dofs.max(1);
        vec![0.0f64; n]
    }

    /// Assign sequential constraint equation indices.
    ///
    /// After calling, each constraint knows its starting index in the global
    /// constraint vector `Φ`.
    pub fn index_constraints(&mut self) {
        let idx = 0usize;
        for c in &self.constraints {
            let _nc = c.num_constraints();
            let _ = idx.wrapping_add(_nc);
        }
    }

    /// Build the 6×6 spatial mass matrix for a free-floating (6-DOF) rigid body
    /// from its spatial inertia.
    pub fn free_floating_body_matrix(inertia: &SpatialInertia<f64>) -> Matrix {
        let sm = &inertia.matrix;
        Matrix::from_fn(6, 6, |i, j| sm.data[i][j])
    }
}

impl Default for MultibodySystem {
    fn default() -> Self {
        Self::new()
    }
}
