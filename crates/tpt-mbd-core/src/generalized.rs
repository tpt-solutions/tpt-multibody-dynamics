//! Generalized coordinate types for multibody dynamics.
//!
//! - [`GeneralizedCoordinates`] — joint angles/positions `q`.
//! - [`GeneralizedVelocities`] — joint rates `q̇`.
//! - [`GeneralizedAccelerations`] — joint accelerations `q̈`.
//!
//! These are thin newtype wrappers over a fixed-size `[f64; N]` array,
//! making them `no_std`-compatible without `alloc`.

use core::fmt;

/// Generalized coordinates `q`: joint angles (revolute) or displacements
/// (prismatic).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GeneralizedCoordinates<const N: usize> {
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

    /// Number of DOFs.
    pub fn len(&self) -> usize {
        N
    }
    pub fn is_empty(&self) -> bool {
        N == 0
    }
}

/// Generalized velocities `q̇`: joint rates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GeneralizedVelocities<const N: usize> {
    pub data: [f64; N],
}

impl<const N: usize> GeneralizedVelocities<N> {
    pub fn zero() -> Self {
        Self { data: [0.0; N] }
    }
    pub fn new(data: [f64; N]) -> Self {
        Self { data }
    }
    pub fn len(&self) -> usize {
        N
    }
    pub fn is_empty(&self) -> bool {
        N == 0
    }
}

/// Generalized accelerations `q̈`: joint accelerations.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GeneralizedAccelerations<const N: usize> {
    pub data: [f64; N],
}

impl<const N: usize> GeneralizedAccelerations<N> {
    pub fn zero() -> Self {
        Self { data: [0.0; N] }
    }
    pub fn new(data: [f64; N]) -> Self {
        Self { data }
    }
    pub fn len(&self) -> usize {
        N
    }
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
