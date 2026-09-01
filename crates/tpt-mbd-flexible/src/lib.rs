#![doc = include_str!("../README.md")]
#![warn(missing_docs)]
#![forbid(unsafe_code)]

//! Flexible multibody dynamics via component mode synthesis.
//!
//! Implements the Craig-Bampton method, modal superposition, floating
//! frame formulation, absolute nodal coordinate formulation (ANCF),
//! Rayleigh damping, and seamless integration with `tpt-fem-elasticity`
//! and `tpt-fem-eigen`.

pub mod cms;
pub mod floating_frame;
pub mod ancf;
pub mod damping;
