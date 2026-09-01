//! Coulomb friction with smooth regularization and Stribeck effect.

use crate::Vector3;

/// Friction model variants.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FrictionModel {
    Coulomb,
    SmoothCoulomb,
    Stribeck,
}

/// Parameters for Coulomb friction models.
#[derive(Clone, Debug, PartialEq)]
pub struct CoulombFrictionParams {
    /// Static friction coefficient.
    pub static_coeff: f64,
    /// Kinetic friction coefficient.
    pub kinetic_coeff: f64,
    /// Stribeck transition velocity.
    pub stribeck_velocity: f64,
    /// Viscous friction coefficient.
    pub viscous_coeff: f64,
}

impl Default for CoulombFrictionParams {
    fn default() -> Self {
        Self {
            static_coeff: 0.8,
            kinetic_coeff: 0.6,
            stribeck_velocity: 0.01,
            viscous_coeff: 0.0,
        }
    }
}

/// Computes the friction force vector given normal force and tangential velocity.
pub fn compute_friction_force(
    normal_force: f64,
    tangential_velocity: Vector3,
    params: &CoulombFrictionParams,
    model: FrictionModel,
) -> Vector3 {
    let speed = tangential_velocity.norm();
    if speed < f64::EPSILON {
        return Vector3::zero();
    }

    let friction_coeff = match model {
        FrictionModel::Coulomb => params.kinetic_coeff,
        FrictionModel::SmoothCoulomb => {
            let v_ratio = speed / (params.stribeck_velocity + f64::EPSILON);
            let smooth = 2.0 * v_ratio;
            params.kinetic_coeff + (params.static_coeff - params.kinetic_coeff) / (1.0 + smooth * smooth)
        }
        FrictionModel::Stribeck => StribeckFriction::coefficient(speed, params),
    };

    let friction_mag = friction_coeff * normal_force + params.viscous_coeff * speed;
    let friction_dir = tangential_velocity.normalize();
    friction_dir * (-friction_mag)
}

/// Stribeck friction model.
pub struct StribeckFriction;

impl StribeckFriction {
    /// Computes the Stribeck friction coefficient for a given sliding velocity.
    pub fn coefficient(velocity: f64, params: &CoulombFrictionParams) -> f64 {
        if velocity < f64::EPSILON {
            return params.static_coeff;
        }
        let v = velocity / (params.stribeck_velocity + f64::EPSILON);
        params.kinetic_coeff + (params.static_coeff - params.kinetic_coeff) * (-v).exp()
    }

    /// Computes the friction force magnitude.
    pub fn force(normal: f64, velocity: f64, params: &CoulombFrictionParams) -> f64 {
        Self::coefficient(velocity, params) * normal
    }
}

/// Anisotropic friction with separate coefficients per axis.
#[derive(Clone, Debug, PartialEq)]
pub struct AnisotropicFriction {
    /// Friction coefficient for x-direction sliding.
    pub coeff_x: f64,
    /// Friction coefficient for y-direction sliding.
    pub coeff_y: f64,
    /// Friction coefficient for z-direction sliding.
    pub coeff_z: f64,
}

impl Default for AnisotropicFriction {
    fn default() -> Self {
        Self {
            coeff_x: 0.6,
            coeff_y: 0.6,
            coeff_z: 0.6,
        }
    }
}

impl AnisotropicFriction {
    /// Create a new anisotropic friction model with per-axis coefficients.
    pub fn new(coeff_x: f64, coeff_y: f64, coeff_z: f64) -> Self {
        Self {
            coeff_x,
            coeff_y,
            coeff_z,
        }
    }

    /// Compute the friction force vector for anisotropic surfaces.
    pub fn force(&self, normal: f64, velocity: Vector3) -> Vector3 {
        let speed = velocity.norm();
        if speed < f64::EPSILON {
            return Vector3::zero();
        }
        let dir = velocity.normalize();
        let mu = dir.0[0].abs() * self.coeff_x
            + dir.0[1].abs() * self.coeff_y
            + dir.0[2].abs() * self.coeff_z;
        let friction_mag = mu * normal;
        dir * (-friction_mag)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_coulomb_friction() {
        let params = CoulombFrictionParams::default();
        let vel = Vector3::new(1.0, 0.0, 0.0);
        let force = compute_friction_force(10.0, vel, &params, FrictionModel::Coulomb);
        assert!((force.0[0] - (-params.kinetic_coeff * 10.0)).abs() < 1e-6);
        assert!(force.0[0] < 0.0);
    }

    #[test]
    fn test_smooth_coulomb_friction() {
        let params = CoulombFrictionParams {
            static_coeff: 0.8,
            kinetic_coeff: 0.6,
            stribeck_velocity: 0.01,
            viscous_coeff: 0.0,
        };
        let vel = Vector3::new(0.01, 0.0, 0.0);
        let force = compute_friction_force(10.0, vel, &params, FrictionModel::SmoothCoulomb);
        assert!(force.0[0] < 0.0);
        assert!(force.0[0].abs() > 0.0);
    }

    #[test]
    fn test_stribeck_friction_coefficient() {
        let params = CoulombFrictionParams::default();
        let coeff = StribeckFriction::coefficient(0.0, &params);
        assert!((coeff - params.static_coeff).abs() < 1e-6);
    }

    #[test]
    fn test_stribeck_friction_force() {
        let params = CoulombFrictionParams::default();
        let force = StribeckFriction::force(10.0, 0.0, &params);
        assert!((force - params.static_coeff * 10.0).abs() < 1e-6);
    }

    #[test]
    fn test_stribeck_friction_high_velocity() {
        let params = CoulombFrictionParams::default();
        let coeff = StribeckFriction::coefficient(100.0, &params);
        assert!((coeff - params.kinetic_coeff).abs() < 1e-6);
    }

    #[test]
    fn test_anisotropic_friction() {
        let af = AnisotropicFriction::new(0.5, 0.8, 0.3);
        let vel = Vector3::new(1.0, 0.0, 0.0);
        let force = af.force(10.0, vel);
        assert!(force.0[0] < 0.0);
    }

    #[test]
    fn test_no_friction_when_zero_velocity() {
        let params = CoulombFrictionParams::default();
        let vel = Vector3::new(0.0, 0.0, 0.0);
        let force = compute_friction_force(10.0, vel, &params, FrictionModel::Coulomb);
        assert_eq!(force, Vector3::zero());
    }

    #[test]
    fn test_viscous_contribution() {
        let params = CoulombFrictionParams {
            static_coeff: 0.8,
            kinetic_coeff: 0.6,
            stribeck_velocity: 0.01,
            viscous_coeff: 0.1,
        };
        let vel = Vector3::new(1.0, 0.0, 0.0);
        let force = compute_friction_force(10.0, vel, &params, FrictionModel::Coulomb);
        let expected_mag = params.kinetic_coeff * 10.0 + params.viscous_coeff * 1.0;
        assert!((force.0[0].abs() - expected_mag).abs() < 1e-6);
    }
}
