#![doc = include_str!("../../../README.md")]
#![allow(missing_docs)]
#![allow(clippy::all)]
#![forbid(unsafe_code)]

//! Flexible multibody dynamics via component mode synthesis.
//!
//! This crate provides the building blocks for embedding flexible bodies inside
//! a multibody dynamics simulation:
//!
//! * **Craig-Bampton (CMS)** — boundary/interior DOF partitioning, fixed-
//!   interface normal modes, constraint modes, and reduced mass/stiffness
//!   assembly (`cms` module).
//! * **Rayleigh damping** — mass- and stiffness-proportional modal damping
//!   (`damping` module).
//! * **Floating frame of reference (FFR)** — large rigid-body motion plus
//!   small elastic deformation, gyroscopic coupling, and deformation-gradient
//!   evaluation (`floating_frame` module).
//! * **Absolute Nodal Coordinate Formulation (ANCF)** — consistent mass and
//!   stiffness matrices for beam/shell elements suitable for large-deformation
//!   gradient-deficient formulations (`ancf` module).
//!
//! All modules are pure Rust (no `unsafe`), use the in-house dense linear-
//! algebra backend from `tpt-math-linalg-dense`, and carry unit tests that
//! exercise each formulation against hand-computed or structural mechanics
//! reference values.

pub mod ancf;
pub mod cms;
pub mod damping;
pub mod floating_frame;
