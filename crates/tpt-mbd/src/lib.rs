#![doc = include_str!("../../../README.md")]
#![warn(missing_docs)]
#![forbid(unsafe_code)]

//! Engineering-grade multibody dynamics for robotics, vehicle dynamics,
//! biomechanics, and mechanism simulation.
//!
//! This umbrella crate re-exports capabilities from the constituent crates.
//! Compiling with no features yields a minimal crate re-exporting only
//! `tpt-mbd-core`:
//!
//! ```rust
//! # #[cfg(feature = "core")] {
//! use tpt_mbd::tpt_mbd_core;
//! # }
//! ```
//!
//! Enable individual feature flags to pull in additional solvers:
//!
//! ```toml
//! [dependencies]
//! tpt-mbd = { version = "0.1", features = ["kinematics", "joints", "system"] }
//! ```
//!
//! Available features:
//! - `core` — spatial algebra, frames, inertia (always available)
//! - `kinematics` — forward/inverse kinematics, Jacobians
//! - `joints` — joint types, constraint formulation, stabilization
//! - `contact` — collision detection, Hertzian contact, friction, impact
//! - `flexible` — Craig-Bampton CMS, modal superposition
//! - `system` — system assembly, time integration, actuators

#[cfg(feature = "core")]
pub use tpt_mbd_core;

#[cfg(feature = "kinematics")]
pub use tpt_mbd_kinematics;

#[cfg(feature = "joints")]
pub use tpt_mbd_joints;

#[cfg(feature = "contact")]
pub use tpt_mbd_contact;

#[cfg(feature = "flexible")]
pub use tpt_mbd_flexible;

#[cfg(feature = "system")]
pub use tpt_mbd_system;
