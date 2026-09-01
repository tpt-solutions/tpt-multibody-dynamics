#![doc = include_str!("../README.md")]
#![warn(missing_docs)]
#![forbid(unsafe_code)]

//! Joint constraint formulation and stabilization for multibody systems.
//!
//! Implements holonomic joint types (revolute, prismatic, spherical, etc.),
//! constraint equations, Baumgarte stabilization, coordinate partitioning,
//! non-holonomic constraints, joint friction, and limit enforcement.

pub mod joint;
pub mod constraint;
pub mod stabilization;
pub mod friction;
