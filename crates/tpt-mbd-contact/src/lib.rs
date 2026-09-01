#![doc = include_str!("../../../README.md")]
#![allow(missing_docs)]
#![allow(clippy::all)]
#![forbid(unsafe_code)]

//! Contact mechanics and collision detection for multibody systems.
//!
//! Implements continuous and discrete collision detection (GJK/EPA),
//! Hertzian and penalty contact force models, Coulomb friction, impact
//! handling, complementarity-based contact, and Archard wear.

pub mod ccd;
pub mod contact;
pub mod detection;
pub mod friction;
pub mod impact;
pub mod wear;

use core::ops::{Add, Div, Mul, Sub};

/// 3D vector type used throughout the contact mechanics modules.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Vector3(pub [f64; 3]);

impl Vector3 {
    /// Create a new vector from x, y, z components.
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self([x, y, z])
    }

    /// Zero vector.
    pub fn zero() -> Self {
        Self([0.0; 3])
    }

    /// Dot product with another vector.
    pub fn dot(&self, other: &Self) -> f64 {
        self.0[0] * other.0[0] + self.0[1] * other.0[1] + self.0[2] * other.0[2]
    }

    /// Cross product with another vector.
    pub fn cross(&self, other: &Self) -> Self {
        Self([
            self.0[1] * other.0[2] - self.0[2] * other.0[1],
            self.0[2] * other.0[0] - self.0[0] * other.0[2],
            self.0[0] * other.0[1] - self.0[1] * other.0[0],
        ])
    }

    /// Euclidean norm (length).
    pub fn norm(&self) -> f64 {
        self.dot(self).sqrt()
    }

    /// Unit vector in the same direction, or zero if self is zero.
    pub fn normalize(&self) -> Self {
        let n = self.norm();
        if n > f64::EPSILON {
            Self([self.0[0] / n, self.0[1] / n, self.0[2] / n])
        } else {
            Self::zero()
        }
    }

    /// Component-wise multiplication.
    pub fn element_mul(&self, other: &Self) -> Self {
        Self([
            self.0[0] * other.0[0],
            self.0[1] * other.0[1],
            self.0[2] * other.0[2],
        ])
    }
}

impl Add for Vector3 {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self([
            self.0[0] + rhs.0[0],
            self.0[1] + rhs.0[1],
            self.0[2] + rhs.0[2],
        ])
    }
}

impl Sub for Vector3 {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self([
            self.0[0] - rhs.0[0],
            self.0[1] - rhs.0[1],
            self.0[2] - rhs.0[2],
        ])
    }
}

impl Mul<f64> for Vector3 {
    type Output = Self;
    fn mul(self, scalar: f64) -> Self {
        Self([self.0[0] * scalar, self.0[1] * scalar, self.0[2] * scalar])
    }
}

impl Div<f64> for Vector3 {
    type Output = Self;
    fn div(self, scalar: f64) -> Self {
        Self([self.0[0] / scalar, self.0[1] / scalar, self.0[2] / scalar])
    }
}
