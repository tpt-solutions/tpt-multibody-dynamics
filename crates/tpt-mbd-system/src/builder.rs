//! Builder pattern for `MultibodySystem`.
//!
//! Provides a fluent API for constructing multibody systems:
//!
//! ```
//! use tpt_mbd_system::builder::MultibodySystemBuilder;
//! use tpt_mbd_system::system::MultibodySystem;
//! use tpt_mbd_core::{RigidBody, SpatialInertia};
//! use tpt_mbd_core::frame::Isometry3;
//! use tpt_math_linalg_fixed::{Matrix3, Vector3};
//!
//! let si = SpatialInertia::new(
//!     1.0,
//!     Vector3::new([0.0, 0.0, 0.0]),
//!     Matrix3::new([[1.0; 3]; 3]),
//! );
//! let body = RigidBody::new(si, Isometry3::identity(), "link0", 0);
//! let system = MultibodySystemBuilder::new()
//!     .add_body(body)
//!     .build();
//! assert_eq!(system.bodies.len(), 1);
//! ```

extern crate alloc;

use alloc::boxed::Box;

use tpt_mbd_core::inertia::RigidBody;
use tpt_mbd_joints::constraint::JointConstraint;
use tpt_mbd_joints::joint::JointType;

use crate::system::{ContactManifold, FloatingFrameBody, MultibodySystem};

/// Fluent builder for [`MultibodySystem`].
///
/// Accumulates bodies, joints, constraints, contact manifolds, and flexible
/// bodies, then produces a fully initialized [`MultibodySystem`] on `build()`.
#[derive(Debug, Default)]
pub struct MultibodySystemBuilder {
    bodies: alloc::vec::Vec<RigidBody<f64>>,
    joints: alloc::vec::Vec<JointType>,
    joint_body_pairs: alloc::vec::Vec<(usize, usize)>,
    constraints: alloc::vec::Vec<Box<dyn JointConstraint>>,
    contacts: alloc::vec::Vec<ContactManifold>,
    flexible: alloc::vec::Vec<FloatingFrameBody>,
}

impl MultibodySystemBuilder {
    /// Create a new empty builder.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a rigid body to the system.
    ///
    /// Returns the builder for chaining.
    pub fn add_body(mut self, body: RigidBody<f64>) -> Self {
        self.bodies.push(body);
        self
    }

    /// Add a joint connecting two bodies.
    ///
    /// `body_i` and `body_j` are indices into the bodies vector (the order
    /// in which bodies were added via `add_body`).
    ///
    /// Returns the builder for chaining.
    pub fn add_joint(mut self, joint: JointType, body_i: usize, body_j: usize) -> Self {
        self.joints.push(joint);
        self.joint_body_pairs.push((body_i, body_j));
        self
    }

    /// Add a holonomic constraint to the system.
    ///
    /// Returns the builder for chaining.
    pub fn add_constraint(mut self, constraint: Box<dyn JointConstraint>) -> Self {
        self.constraints.push(constraint);
        self
    }

    /// Add a contact manifold to the system.
    ///
    /// Returns the builder for chaining.
    pub fn add_contact(mut self, contact: ContactManifold) -> Self {
        self.contacts.push(contact);
        self
    }

    /// Add a flexible body with reduced modal coordinates to the system.
    ///
    /// Returns the builder for chaining.
    pub fn add_flexible_body(mut self, body: FloatingFrameBody) -> Self {
        self.flexible.push(body);
        self
    }

    /// Build the [`MultibodySystem`].
    ///
    /// Computes `num_dofs` from the accumulated bodies, joints, and flexible
    /// bodies.
    pub fn build(self) -> MultibodySystem {
        let mut system = MultibodySystem {
            bodies: self.bodies,
            joints: self.joints,
            constraints: self.constraints,
            contacts: self.contacts,
            flexible: self.flexible,
            num_dofs: 0,
            joint_body_pairs: self.joint_body_pairs,
        };
        system.count_dofs();
        system
    }
}
