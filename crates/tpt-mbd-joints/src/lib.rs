//! Joint constraint formulation and stabilization for multibody systems.
//!
//! Implements holonomic joint types (revolute, prismatic, spherical, etc.),
//! constraint equations, Baumgarte stabilization, coordinate partitioning,
//! non-holonomic constraints, joint friction, limit enforcement, reaction
//! force computation, and drift detection.
//!
//! # Quick start
//!
//! ```
//! use tpt_mbd_joints::constraint::{JointConstraint, RevoluteConstraint};
//! use tpt_mbd_joints::stabilization::BaumgarteStabilization;
//!
//! // Define a revolute joint between bodies 0 and 1
//! let rev = RevoluteConstraint::new(0, 1, tpt_mbd_joints::joint::JointAxis::Z);
//!
//! // Evaluate constraint at zero configuration
//! let q = vec![0.0; 12];
//! let violation = rev.violation(&q);
//! assert!(violation < 1e-10);
//!
//! // Baumgarte stabilization
//! let stab = BaumgarteStabilization::default();
//! let acc = stab.stabilized_acceleration(0.01, 0.0);
//! ```

#![doc = include_str!("../../../README.md")]
#![allow(missing_docs)]
#![allow(clippy::all)]
#![forbid(unsafe_code)]
#![no_std]

extern crate alloc;

pub mod constraint;
pub mod drift;
pub mod formulation;
pub mod friction;
pub mod joint;
pub mod nonholonomic;
pub mod reaction;
pub mod stabilization;
