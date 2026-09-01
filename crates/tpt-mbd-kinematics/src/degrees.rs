//! Degree and radian unit-safe wrappers for angular quantities.
//!
//! Provides newtype wrappers that prevent accidental mixing of degree and
//! radian values in kinematic computations.

/// Angle in degrees.
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct Degrees(pub f64);

impl Degrees {
    /// Create a new `Degrees` value.
    pub fn new(deg: f64) -> Self {
        Self(deg)
    }

    /// Convert to radians.
    pub fn to_radians(&self) -> Radians {
        Radians(self.0 * core::f64::consts::PI / 180.0)
    }

    /// Get the raw degree value.
    pub fn value(&self) -> f64 {
        self.0
    }
}

impl core::ops::Add for Degrees {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self(self.0 + rhs.0)
    }
}

impl core::ops::Sub for Degrees {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self(self.0 - rhs.0)
    }
}

impl core::ops::Mul<f64> for Degrees {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Self(self.0 * rhs)
    }
}

impl core::ops::Div<f64> for Degrees {
    type Output = Self;
    fn div(self, rhs: f64) -> Self {
        Self(self.0 / rhs)
    }
}

impl From<Degrees> for f64 {
    fn from(d: Degrees) -> f64 {
        d.0
    }
}

impl From<f64> for Degrees {
    fn from(deg: f64) -> Self {
        Self(deg)
    }
}

/// Angle in radians.
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct Radians(pub f64);

impl Radians {
    /// Create a new `Radians` value.
    pub fn new(rad: f64) -> Self {
        Self(rad)
    }

    /// Convert to degrees.
    pub fn to_degrees(&self) -> Degrees {
        Degrees(self.0 * 180.0 / core::f64::consts::PI)
    }

    /// Get the raw radian value.
    pub fn value(&self) -> f64 {
        self.0
    }
}

impl core::ops::Add for Radians {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self(self.0 + rhs.0)
    }
}

impl core::ops::Sub for Radians {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self(self.0 - rhs.0)
    }
}

impl core::ops::Mul<f64> for Radians {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Self(self.0 * rhs)
    }
}

impl core::ops::Div<f64> for Radians {
    type Output = Self;
    fn div(self, rhs: f64) -> Self {
        Self(self.0 / rhs)
    }
}

impl From<Radians> for f64 {
    fn from(r: Radians) -> f64 {
        r.0
    }
}

impl From<f64> for Radians {
    fn from(rad: f64) -> Self {
        Self(rad)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn degrees_to_radians() {
        let deg = Degrees(180.0);
        let rad = deg.to_radians();
        assert!((rad.value() - core::f64::consts::PI).abs() < 1e-12);
    }

    #[test]
    fn radians_to_degrees() {
        let rad = Radians(core::f64::consts::PI);
        let deg = rad.to_degrees();
        assert!((deg.value() - 180.0).abs() < 1e-12);
    }

    #[test]
    fn round_trip() {
        let deg = Degrees(45.0);
        let rad = deg.to_radians();
        let deg2 = rad.to_degrees();
        assert!((deg2.value() - deg.value()).abs() < 1e-12);
    }
}
