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
//!
//! # Examples
//!
//! ```
//! use tpt_mbd_flexible::cms::{CraigBampton, ModeSelection};
//! use tpt_mbd_flexible::superposition::{ModalSuperpositionState, project_force};
//! use tpt_math_linalg_dense::DMatrix;
//!
//! // Build a simple 2-mode reduced model
//! let phi = DMatrix::from_row_slice(4, 2, &[1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]);
//! let m_red = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0]);
//! let k_red = DMatrix::from_row_slice(2, 2, &[100.0, 0.0, 0.0, 400.0]);
//! let c_red = DMatrix::from_row_slice(2, 2, &[0.5, 0.0, 0.0, 1.0]);
//! let mut state = ModalSuperpositionState::from_reduced_matrices(m_red, k_red, c_red);
//! state.q = vec![1.0, 2.0];
//! let u = state.displacement(&phi);
//! assert_eq!(u.len(), 4);
//! ```

pub mod ancf;
pub mod cms;
pub mod damping;
pub mod floating_frame;
pub mod superposition;
