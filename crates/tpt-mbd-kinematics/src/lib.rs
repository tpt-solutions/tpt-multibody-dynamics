#![doc = include_str!("../../../README.md")]
#![warn(missing_docs)]
#![forbid(unsafe_code)]

//! Forward and inverse kinematics for serial and parallel kinematic chains.
//!
//! Implements Denavit-Hartenberg (DH) parameter convention, Product of
//! Exponentials (PoE) formulation, forward/inverse kinematics solvers,
//! geometric and analytical Jacobians, singularity detection, and workspace
//! analysis.

pub mod chain;
pub mod forward;
pub mod inverse;
pub mod jacobian;
pub mod singularity;
