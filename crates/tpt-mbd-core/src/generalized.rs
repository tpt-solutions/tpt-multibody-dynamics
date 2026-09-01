//! Generalized coordinate types for multibody dynamics.
//!
//! - [`GeneralizedCoordinates`] — joint angles/positions `q`.
//! - [`GeneralizedVelocities`] — joint rates `q̇`.
//! - [`GeneralizedAccelerations`] — joint accelerations `q̈`.
//!
//! These are thin newtype wrappers over a fixed-size `[f64; N]` array,
//! making them `no_std`-compatible without `alloc`.
//!
//! # Examples
//!
//! ```
//! use tpt_mbd_core::GeneralizedCoordinates;
//!
//! let q = GeneralizedCoordinates::<3>::zero();
//! assert_eq!(q.len(), 3);
//! assert!(!q.is_empty());
//! ```
//!
//! ```
//! use tpt_mbd_core::GeneralizedCoordinates;
//!
//! let q = GeneralizedCoordinates::new([0.1, -0.2, 0.5]);
//! assert_eq!(q.data[1], -0.2);
//! ```
//!
//! ```
//! use tpt_mbd_core::{GeneralizedCoordinates, GeneralizedVelocities, GeneralizedAccelerations};
//!
//! let q  = GeneralizedCoordinates::<2>::new([0.0, 1.0]);
//! let qd = GeneralizedVelocities::<2>::new([0.5, -0.5]);
//! let qdd= GeneralizedAccelerations::<2>::new([1.0, 0.0]);
//! assert_eq!(q.len(), 2);
//! assert_eq!(qd.len(), 2);
//! assert_eq!(qdd.len(), 2);
//! ```

use core::fmt;

/// Generalized coordinates `q`: joint angles (revolute) or displacements
/// (prismatic).
///
/// Stored as a fixed-size `[f64; N]` array, making this type `no_std`
/// compatible without allocation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GeneralizedCoordinates<const N: usize> {
    /// Raw coordinate values.
    pub data: [f64; N],
}

impl<const N: usize> GeneralizedCoordinates<N> {
    /// Zero-initialized coordinates.
    pub fn zero() -> Self {
        Self { data: [0.0; N] }
    }

    /// Build from a raw array.
    pub fn new(data: [f64; N]) -> Self {
        Self { data }
    }

    /// Number of generalized coordinates.
    pub fn len(&self) -> usize {
        N
    }

    /// Whether there are zero coordinates.
    pub fn is_empty(&self) -> bool {
        N == 0
    }
}

/// Generalized velocities `q̇`: joint rates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GeneralizedVelocities<const N: usize> {
    /// Raw velocity values.
    pub data: [f64; N],
}

impl<const N: usize> GeneralizedVelocities<N> {
    /// Zero-initialized velocities.
    pub fn zero() -> Self {
        Self { data: [0.0; N] }
    }
    /// Build from a raw array.
    pub fn new(data: [f64; N]) -> Self {
        Self { data }
    }
    /// Number of generalized velocities.
    pub fn len(&self) -> usize {
        N
    }
    /// Whether there are zero velocities.
    pub fn is_empty(&self) -> bool {
        N == 0
    }
}

/// Generalized accelerations `q̈`: joint accelerations.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GeneralizedAccelerations<const N: usize> {
    /// Raw acceleration values.
    pub data: [f64; N],
}

impl<const N: usize> GeneralizedAccelerations<N> {
    /// Zero-initialized accelerations.
    pub fn zero() -> Self {
        Self { data: [0.0; N] }
    }
    /// Build from a raw array.
    pub fn new(data: [f64; N]) -> Self {
        Self { data }
    }
    /// Number of generalized accelerations.
    pub fn len(&self) -> usize {
        N
    }
    /// Whether there are zero accelerations.
    pub fn is_empty(&self) -> bool {
        N == 0
    }
}

impl<const N: usize> fmt::Display for GeneralizedCoordinates<N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "q = {:?}", self.data)
    }
}

impl<const N: usize> fmt::Display for GeneralizedVelocities<N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "q̇ = {:?}", self.data)
    }
}

impl<const N: usize> fmt::Display for GeneralizedAccelerations<N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "q̈ = {:?}", self.data)
    }
}
