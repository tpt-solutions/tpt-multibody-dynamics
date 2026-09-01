#![doc = include_str!("../../../README.md")]
#![warn(missing_docs)]
#![forbid(unsafe_code)]

//! System assembly and time integration for complete multibody simulations.
//!
//! Implements `MultibodySystem` assembly, minimal and maximal coordinate
//! formulations, Featherstone's articulated body algorithm, explicit and
//! implicit integrators, external force application, actuator models,
//! system linearization, and parallel evaluation via Rayon.
//!
//! # Examples
//!
//! ```
//! use tpt_mbd_system::system::MultibodySystem;
//! use tpt_mbd_core::{RigidBody, SpatialInertia};
//! use tpt_mbd_core::frame::Isometry3;
//! use tpt_math_linalg_fixed::{Matrix3, Vector3};
//!
//! let mut sys = MultibodySystem::new();
//! let si = SpatialInertia::new(
//!     1.0,
//!     Vector3::new([0.0, 0.0, 0.0]),
//!     Matrix3::new([[1.0; 3]; 3]),
//! );
//! let body = RigidBody::new(si, Isometry3::identity(), "link0", 0);
//! sys.add_body(body);
//! assert_eq!(sys.bodies.len(), 1);
//! ```

pub mod actuators;
pub mod forces;
pub mod integration;
pub mod system;

/// Dense column-major matrix alias (f64).
pub type Matrix = tpt_math_linalg_dense::DMatrix<f64>;
/// Dense column-major vector alias (f64).
pub type Vector = tpt_math_linalg_dense::DVector<f64>;