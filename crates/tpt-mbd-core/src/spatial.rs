//! Domain-specific spatial vector wrappers for multibody dynamics.
//!
//! Wraps [`tpt_math_spatial`] types in newtypes so we can add `Display`
//! and domain-specific methods without violating the orphan rule:
//!
//! - [`SpatialVelocity`] — twist `(ω, v)` — wraps `MotionVector`.
//! - [`SpatialForce`] — wrench `(τ, f)` — wraps `ForceVector`.
//! - [`SpatialMomentum`] — spatial momentum `(h, l)` — wraps `ForceVector`.
//!
//! Cross-product operators follow Featherstone's `crm` / `crf` convention.

use core::fmt;

use tpt_math_linalg_fixed::Vector3;
use tpt_math_numeric::Scalar;
use tpt_math_spatial::{ForceVector, MotionVector};

// ===========================================================================
// Newtype wrappers (newtypes avoid orphan-rule issues)
// ===========================================================================

/// A spatial **velocity** (twist): angular velocity `ω` (top) and linear
/// velocity `v` (bottom), stored `[ω; v]`.
///
/// Wraps [`tpt_math_spatial::MotionVector`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SpatialVelocity<T: Copy>(pub MotionVector<T>);

/// A spatial **force** (wrench): moment `τ` (top) and force `f` (bottom),
/// stored `[τ; f]`.
///
/// Wraps [`tpt_math_spatial::ForceVector`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SpatialForce<T: Copy>(pub ForceVector<T>);

/// Spatial **momentum** `(h, l)`: angular momentum `h` (top) and linear
/// momentum `l` (bottom), stored `[h; l]`.
///
/// Wraps `ForceVector` (same storage convention as [`SpatialForce`]).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SpatialMomentum<T: Copy>(pub ForceVector<T>);

// ===========================================================================
// Constructors
// ===========================================================================

impl<T: Copy + Scalar> SpatialVelocity<T> {
    /// Build a spatial velocity from angular and linear parts.
    pub fn new(angular: Vector3<T>, linear: Vector3<T>) -> Self {
        Self(MotionVector::new(angular, linear))
    }

    /// Borrow the underlying [`MotionVector`].
    pub fn inner(&self) -> &MotionVector<T> {
        &self.0
    }

    /// Angular part (top 3 components).
    pub fn angular(&self) -> Vector3<T> {
        self.0.angular()
    }

    /// Linear part (bottom 3 components).
    pub fn linear(&self) -> Vector3<T> {
        self.0.linear()
    }
}

impl<T: Copy + Scalar> SpatialForce<T> {
    /// Build a spatial force from its angular (moment) and linear (force) parts.
    pub fn new(angular: Vector3<T>, linear: Vector3<T>) -> Self {
        Self(ForceVector::new(angular, linear))
    }

    /// Borrow the underlying [`ForceVector`].
    pub fn inner(&self) -> &ForceVector<T> {
        &self.0
    }

    /// Angular part (top 3 — moment).
    pub fn angular(&self) -> Vector3<T> {
        self.0.angular()
    }

    /// Linear part (bottom 3 — force).
    pub fn linear(&self) -> Vector3<T> {
        self.0.linear()
    }
}

impl<T: Copy + Scalar> SpatialMomentum<T> {
    /// Build spatial momentum from angular and linear momentum parts.
    pub fn new(angular: Vector3<T>, linear: Vector3<T>) -> Self {
        Self(ForceVector::new(angular, linear))
    }

    /// Angular part (top 3 — angular momentum).
    pub fn angular(&self) -> Vector3<T> {
        self.0.angular()
    }

    /// Linear part (bottom 3 — linear momentum).
    pub fn linear(&self) -> Vector3<T> {
        self.0.linear()
    }
}

// ===========================================================================
// Arithmetic forwarding
// ===========================================================================

impl<T: Scalar + Copy> Add for SpatialVelocity<T> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self(self.0 + rhs.0)
    }
}

impl<T: Scalar + Copy> Sub for SpatialVelocity<T> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self(self.0 - rhs.0)
    }
}

impl<T: Scalar + Copy> Neg for SpatialVelocity<T> {
    type Output = Self;
    fn neg(self) -> Self {
        Self(-self.0)
    }
}

impl<T: Scalar + Copy> Add for SpatialForce<T> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self(self.0 + rhs.0)
    }
}

impl<T: Scalar + Copy> Sub for SpatialForce<T> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self(self.0 - rhs.0)
    }
}

impl<T: Scalar + Copy> Neg for SpatialForce<T> {
    type Output = Self;
    fn neg(self) -> Self {
        Self(-self.0)
    }
}

impl<T: Scalar + Copy> Mul<T> for SpatialVelocity<T> {
    type Output = Self;
    fn mul(self, rhs: T) -> Self {
        Self(self.0 * rhs)
    }
}

impl<T: Scalar + Copy> Mul<T> for SpatialForce<T> {
    type Output = Self;
    fn mul(self, rhs: T) -> Self {
        Self(self.0 * rhs)
    }
}

// ===========================================================================
// Cross-product helpers (explicit names per spec)
// ===========================================================================

/// Motion cross-product `v×` (Featherstone `crm`): `motion × motion → motion`.
///
/// For two twists `v = (ω₁, v₁)` and `w = (ω₂, v₂)`:
/// ```text
/// v×w = (ω₁ × ω₂, ω₁ × v₂ + v₁ × ω₂)
/// ```
pub fn motion_cross<T: Scalar + Copy>(
    a: &SpatialVelocity<T>,
    b: &SpatialVelocity<T>,
) -> SpatialVelocity<T> {
    SpatialVelocity(a.0.cross(&b.0))
}

/// Force cross-product `v×*` (Featherstone `crf`): `force × force → motion`.
///
/// For two wrenches `f = (τ₁, f₁)` and `g = (τ₂, f₂)`:
/// ```text
/// f×*g = (τ₁ × f₂ + f₁ × τ₂, f₁ × f₂)
/// ```
pub fn force_cross<T: Scalar + Copy>(
    a: &SpatialForce<T>,
    b: &SpatialForce<T>,
) -> SpatialVelocity<T> {
    SpatialVelocity(a.0.cross(&b.0))
}

/// Mixed cross-product `motion × force → force` (Featherstone `crm` on force).
pub fn motion_cross_force<T: Scalar + Copy>(
    m: &SpatialVelocity<T>,
    f: &SpatialForce<T>,
) -> SpatialForce<T> {
    SpatialForce(m.0.cross_force(&f.0))
}

/// Mixed cross-product `force × motion → force` (Featherstone `crf`).
pub fn force_cross_motion<T: Scalar + Copy>(
    f: &SpatialForce<T>,
    m: &SpatialVelocity<T>,
) -> SpatialForce<T> {
    SpatialForce(f.0.cross_motion(&m.0))
}

// ===========================================================================
// Display
// ===========================================================================

impl<T: Copy + fmt::Debug> fmt::Display for SpatialVelocity<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "v = (ω: {:?}, v: {:?})",
            &self.0.data.data[0..3],
            &self.0.data.data[3..6]
        )
    }
}

impl<T: Copy + fmt::Debug> fmt::Display for SpatialForce<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "f = (τ: {:?}, f: {:?})",
            &self.0.data.data[0..3],
            &self.0.data.data[3..6]
        )
    }
}

impl<T: Copy + fmt::Debug> fmt::Display for SpatialMomentum<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "h = (hω: {:?}, l: {:?})",
            &self.0.data.data[0..3],
            &self.0.data.data[3..6]
        )
    }
}

use core::ops::{Add, Mul, Neg, Sub};
