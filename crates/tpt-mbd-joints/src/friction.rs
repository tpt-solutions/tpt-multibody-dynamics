//! Joint friction models: Coulomb + viscous with Stribeck regularization.
//!
//! Provides:
//! - Coulomb friction with smooth regularization
//! - Viscous damping
//! - Stribeck friction model with static/kinetic transition
//! - Combined friction model

#![allow(missing_docs)]

use core::fmt;
use num_traits::real::Real;
use num_traits::Float;

/// Friction model parameters for a joint.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrictionParams {
    pub coulomb_coefficient: f64,
    pub viscous_coefficient: f64,
    pub regularization_velocity: f64,
}

impl Default for FrictionParams {
    fn default() -> Self {
        Self {
            coulomb_coefficient: 0.1,
            viscous_coefficient: 0.01,
            regularization_velocity: 0.001,
        }
    }
}

impl FrictionParams {
    /// Create new friction parameters.
    pub fn new(coulomb: f64, viscous: f64, reg_vel: f64) -> Self {
        Self {
            coulomb_coefficient: coulomb,
            viscous_coefficient: viscous,
            regularization_velocity: reg_vel,
        }
    }
}

/// Joint friction force/torque.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrictionForce {
    pub value: f64,
    pub regime: FrictionRegime,
}

/// Friction regime: static or kinetic (sliding).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FrictionRegime {
    Static,
    Kinetic,
}

/// Compute Coulomb friction force with smooth regularization.
///
/// Uses `F_f = μ·F_n · tanh(v / v_s)` for smooth transition between static
/// and kinetic regimes.
pub fn coulomb_friction_force(
    normal_force: f64,
    velocity: f64,
    params: &FrictionParams,
) -> FrictionForce {
    let mu = params.coulomb_coefficient;
    let vs = params.regularization_velocity;
    let tanh_v = Real::tanh(velocity);
    let tanh_vs = Real::tanh(vs);
    let friction = mu * normal_force * tanh_v / tanh_vs;
    let regime = if velocity.abs() < vs {
        FrictionRegime::Static
    } else {
        FrictionRegime::Kinetic
    };
    FrictionForce {
        value: friction,
        regime,
    }
}

/// Compute viscous friction force (proportional to velocity).
pub fn viscous_friction_force(velocity: f64, viscous_coefficient: f64) -> f64 {
    -viscous_coefficient * velocity
}

/// Combined Coulomb + viscous friction force.
pub fn combined_friction_force(normal_force: f64, velocity: f64, params: &FrictionParams) -> f64 {
    let coulomb = coulomb_friction_force(normal_force, velocity, params).value;
    let viscous = viscous_friction_force(velocity, params.viscous_coefficient);
    coulomb + viscous
}

/// Stribeck friction model with static/kinetic transition.
///
/// The Stribeck effect models the drop from static to kinetic friction
/// as velocity increases: `F = (μ_k + (μ_s - μ_k) · exp(-|v|/v_s)) · F_n · sign(v)`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StribeckFriction {
    pub static_coeff: f64,
    pub kinetic_coeff: f64,
    pub stribeck_velocity: f64,
    pub viscous_coefficient: f64,
}

impl StribeckFriction {
    /// Create a new Stribeck friction model.
    pub fn new(
        static_coeff: f64,
        kinetic_coeff: f64,
        stribeck_velocity: f64,
        viscous_coefficient: f64,
    ) -> Self {
        Self {
            static_coeff,
            kinetic_coeff,
            stribeck_velocity,
            viscous_coefficient,
        }
    }

    /// Compute friction force with Stribeck effect.
    pub fn force(&self, normal_force: f64, velocity: f64) -> f64 {
        let vs = self.stribeck_velocity.max(1e-12);
        let stribeck = self.kinetic_coeff
            + (self.static_coeff - self.kinetic_coeff) * Float::exp(-velocity.abs() / vs);
        let tanh_v = Real::tanh(velocity);
        let tanh_vs = Real::tanh(vs);
        let coulomb = stribeck * normal_force * tanh_v / tanh_vs;
        let viscous = -self.viscous_coefficient * velocity;
        coulomb + viscous
    }
}

impl fmt::Display for FrictionForce {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "FrictionForce({:.4}, {:?})", self.value, self.regime)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_coulomb_friction_static() {
        let params = FrictionParams::new(0.3, 0.01, 0.001);
        let f = coulomb_friction_force(10.0, 0.0005, &params);
        assert_eq!(f.regime, FrictionRegime::Static);
    }

    #[test]
    fn test_coulomb_friction_kinetic() {
        let params = FrictionParams::new(0.3, 0.01, 0.001);
        let f = coulomb_friction_force(10.0, 0.01, &params);
        assert_eq!(f.regime, FrictionRegime::Kinetic);
    }

    #[test]
    fn test_viscous_friction() {
        let f = viscous_friction_force(1.0, 0.5);
        assert_eq!(f, -0.5);
    }

    #[test]
    fn test_combined_friction() {
        let params = FrictionParams::new(0.2, 0.1, 0.01);
        let f = combined_friction_force(10.0, 0.0, &params);
        assert!(f.abs() < 1e-6);
    }

    #[test]
    fn test_stribeck_friction() {
        let stribeck = StribeckFriction::new(0.4, 0.3, 0.01, 0.01);
        let f_static = stribeck.force(10.0, 0.0);
        let f_kinetic = stribeck.force(10.0, 1.0);
        assert!(f_kinetic.abs() > f_static.abs());
    }
}
