#![doc = include_str!("../README.md")]
#![warn(missing_docs)]
#![forbid(unsafe_code)]

//! System assembly and time integration for complete multibody simulations.
//!
//! Implements `MultibodySystem` assembly, minimal and maximal coordinate
//! formulations, Featherstone's articulated body algorithm, explicit and
//! implicit integrators, external force application, actuator models,
//! system linearization, and parallel evaluation via Rayon.

pub mod system;
pub mod integration;
pub mod forces;
pub mod actuators;
