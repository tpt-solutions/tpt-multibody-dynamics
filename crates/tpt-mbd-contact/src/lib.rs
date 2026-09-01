#![doc = include_str!("../README.md")]
#![warn(missing_docs)]
#![forbid(unsafe_code)]

//! Contact mechanics and collision detection for multibody systems.
//!
//! Implements continuous and discrete collision detection (GJK/EPA),
//! Hertzian and penalty contact force models, Coulomb friction, impact
//! handling, complementarity-based contact, and Archard wear.

pub mod ccd;
pub mod detection;
pub mod contact;
pub mod friction;
pub mod impact;
pub mod wear;
