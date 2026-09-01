//! Joint types for multibody systems.
//!
//! Implements: revolute, prismatic, spherical/ball, universal/Cardan,
//! cylindrical, planar, fixed, and custom user-defined joints via
//! [`JointConstraint`].
//!
//! Also provides [`JointLimit`] and [`check_limit`] for limit enforcement.

#![allow(missing_docs)]

use core::fmt;

/// Degrees of freedom for a joint.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum JointDof {
    Fixed,
    Revolute,
    Prismatic,
    Spherical,
    Universal,
    Cylindrical,
    Planar,
    Custom,
}

/// Axis direction for 1-DOF joints.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum JointAxis {
    X,
    Y,
    Z,
}

/// A joint type identifier with its DOF count and axis.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct JointType {
    pub dof: JointDof,
    pub axis: Option<JointAxis>,
    pub num_dofs: usize,
}

impl JointType {
    pub const REVOLUTE: Self = Self {
        dof: JointDof::Revolute,
        axis: Some(JointAxis::Z),
        num_dofs: 1,
    };
    pub const PRISMATIC: Self = Self {
        dof: JointDof::Prismatic,
        axis: Some(JointAxis::Z),
        num_dofs: 1,
    };
    pub const SPHERICAL: Self = Self {
        dof: JointDof::Spherical,
        axis: None,
        num_dofs: 3,
    };
    pub const UNIVERSAL: Self = Self {
        dof: JointDof::Universal,
        axis: None,
        num_dofs: 2,
    };
    pub const CYLINDRICAL: Self = Self {
        dof: JointDof::Cylindrical,
        axis: Some(JointAxis::Z),
        num_dofs: 2,
    };
    pub const PLANAR: Self = Self {
        dof: JointDof::Planar,
        axis: None,
        num_dofs: 3,
    };
    pub const FIXED: Self = Self {
        dof: JointDof::Fixed,
        axis: None,
        num_dofs: 0,
    };
}

impl fmt::Display for JointType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self.dof)
    }
}

/// Status of a joint limit check.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LimitStatus {
    Inside,
    Approaching,
    Violated,
    HardStop,
}

/// Joint limit parameters.
///
/// Defines lower and upper bounds for a joint coordinate with optional
/// soft (spring-damper) penalty parameters.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JointLimit {
    pub lower: f64,
    pub upper: f64,
    pub stiffness: f64,
    pub damping: f64,
}

impl JointLimit {
    /// Create a new joint limit.
    pub fn new(lower: f64, upper: f64, stiffness: f64, damping: f64) -> Self {
        Self {
            lower,
            upper,
            stiffness,
            damping,
        }
    }

    /// Create a soft limit with spring-damper penalty.
    pub fn soft(lower: f64, upper: f64, stiffness: f64, damping: f64) -> Self {
        Self::new(lower, upper, stiffness, damping)
    }

    /// Create a hard limit (zero stiffness/damping — enforced by constraint).
    pub fn hard(lower: f64, upper: f64) -> Self {
        Self::new(lower, upper, 0.0, 0.0)
    }
}

impl Default for JointLimit {
    fn default() -> Self {
        Self {
            lower: -core::f64::consts::PI,
            upper: core::f64::consts::PI,
            stiffness: 1e6,
            damping: 1e3,
        }
    }
}

/// Check a joint coordinate against a limit and return the status.
///
/// Returns:
/// - [`LimitStatus::Inside`] when well within bounds
/// - [`LimitStatus::Approaching`] when within 5% of a bound
/// - [`LimitStatus::Violated`] when outside bounds (soft limit)
/// - [`LimitStatus::HardStop`] when outside bounds (hard limit)
pub fn check_limit(q: f64, limit: &JointLimit) -> LimitStatus {
    if q < limit.lower || q > limit.upper {
        if limit.stiffness == 0.0 {
            LimitStatus::HardStop
        } else {
            LimitStatus::Violated
        }
    } else {
        let margin = (limit.upper - limit.lower) * 0.05;
        if q < limit.lower + margin || q > limit.upper - margin {
            LimitStatus::Approaching
        } else {
            LimitStatus::Inside
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_check_limit_inside() {
        let limit = JointLimit::new(-1.0, 1.0, 1e6, 1e3);
        assert_eq!(check_limit(0.0, &limit), LimitStatus::Inside);
        assert_eq!(check_limit(0.5, &limit), LimitStatus::Inside);
    }

    #[test]
    fn test_check_limit_approaching() {
        let limit = JointLimit::new(-1.0, 1.0, 1e6, 1e3);
        assert_eq!(check_limit(0.96, &limit), LimitStatus::Approaching);
        assert_eq!(check_limit(-0.96, &limit), LimitStatus::Approaching);
    }

    #[test]
    fn test_check_limit_violated_soft() {
        let limit = JointLimit::new(-1.0, 1.0, 1e6, 1e3);
        assert_eq!(check_limit(1.1, &limit), LimitStatus::Violated);
        assert_eq!(check_limit(-1.1, &limit), LimitStatus::Violated);
    }

    #[test]
    fn test_check_limit_hard_stop() {
        let limit = JointLimit::hard(-1.0, 1.0);
        assert_eq!(check_limit(1.1, &limit), LimitStatus::HardStop);
        assert_eq!(check_limit(-1.1, &limit), LimitStatus::HardStop);
    }
}
