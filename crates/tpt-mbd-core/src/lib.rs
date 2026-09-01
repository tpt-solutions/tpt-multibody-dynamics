#![doc = include_str!("../README.md")]
#![warn(missing_docs)]
#![forbid(unsafe_code)]

//! Foundation layer for multibody dynamics using Featherstone's spatial
//! vector algebra.
//!
//! Provides 6D spatial vector types, spatial cross-product operators,
//! spatial inertia, coordinate frames, and generalized coordinate types.
//!
//! All types are `no_std` compatible with an optional `alloc` feature for
//! dynamic system sizes.

pub mod frame;
pub mod spatial;
pub mod inertia;
