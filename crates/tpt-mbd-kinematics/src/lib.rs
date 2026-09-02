#![doc = include_str!("../../../README.md")]
#![allow(missing_docs)]
#![forbid(unsafe_code)]

//! Forward and inverse kinematics for serial and parallel kinematic chains.
//!
//! Implements Denavit-Hartenberg (DH) parameter convention, Product of
//! Exponentials (PoE) formulation, forward/inverse kinematics solvers,
//! geometric and analytical Jacobians, singularity detection, and workspace
//! analysis.

pub mod chain;
pub mod degrees;
pub mod forward;
pub mod inverse;
pub mod jacobian;
pub mod loop_closure;
pub mod pie;
pub mod singularity;
