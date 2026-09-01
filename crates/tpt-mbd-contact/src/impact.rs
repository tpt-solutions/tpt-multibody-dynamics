//! Impact handling: coefficient-of-restitution impulse and soft impact penalty springs.

use crate::Vector3;

/// Parameters for impact response models.
#[derive(Clone, Debug, PartialEq)]
pub struct ImpactParams {
    /// Coefficient of restitution (0 = perfectly inelastic, 1 = perfectly elastic).
    pub restitution_coefficient: f64,
    /// Friction coefficient during impact.
    pub friction_coeff: f64,
}

impl Default for ImpactParams {
    fn default() -> Self {
        Self {
            restitution_coefficient: 0.3,
            friction_coeff: 0.5,
        }
    }
}

/// Impulse-based impact response using the coefficient of restitution.
pub struct CoefficientOfRestitution;

impl CoefficientOfRestitution {
    /// Computes the normal impulse magnitude for a collision.
    ///
    /// j = -(1 + e) * (v_rel · n) / (1/m_i + 1/m_j)
    ///
    /// Returns 0.0 if the bodies are already separating.
    pub fn compute_impulse(
        mass_i: f64,
        mass_j: f64,
        relative_velocity: Vector3,
        normal: Vector3,
        restitution: f64,
    ) -> f64 {
        let v_rel_n = relative_velocity.dot(&normal);
        if v_rel_n >= 0.0 {
            return 0.0;
        }
        let inv_mass_sum = 1.0 / mass_i + 1.0 / mass_j;
        -(1.0 + restitution) * v_rel_n / inv_mass_sum
    }
}

/// Soft impact using a high-stiffness penalty spring.
#[derive(Clone, Debug, PartialEq)]
pub struct SoftImpact {
    /// Penalty stiffness coefficient.
    pub stiffness: f64,
    /// Damping coefficient.
    pub damping: f64,
}

impl Default for SoftImpact {
    fn default() -> Self {
        Self {
            stiffness: 1.0e7,
            damping: 1.0e4,
        }
    }
}

impl SoftImpact {
    /// Create a new soft impact model.
    pub fn new(stiffness: f64, damping: f64) -> Self {
        Self { stiffness, damping }
    }

    /// Compute the soft impact force.
    pub fn force(&self, penetration: f64, penetration_rate: f64) -> f64 {
        if penetration <= 0.0 {
            return 0.0;
        }
        self.stiffness * penetration + self.damping * penetration_rate
    }
}

/// Computes the impact impulse magnitude.
pub fn compute_impact_impulse(
    mass_i: f64,
    mass_j: f64,
    relative_velocity: Vector3,
    normal: Vector3,
    restitution: f64,
) -> f64 {
    CoefficientOfRestitution::compute_impulse(mass_i, mass_j, relative_velocity, normal, restitution)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_impact_impulse_approaching() {
        let rel_vel = Vector3::new(-5.0, 0.0, 0.0);
        let normal = Vector3::new(1.0, 0.0, 0.0);
        let impulse = compute_impact_impulse(1.0, 1.0, rel_vel, normal, 0.5);
        assert!(impulse > 0.0);
    }

    #[test]
    fn test_impact_impulse_separating() {
        let rel_vel = Vector3::new(5.0, 0.0, 0.0);
        let normal = Vector3::new(1.0, 0.0, 0.0);
        let impulse = compute_impact_impulse(1.0, 1.0, rel_vel, normal, 0.5);
        assert_eq!(impulse, 0.0);
    }

    #[test]
    fn test_impact_impulse_equal_masses() {
        let rel_vel = Vector3::new(-4.0, 0.0, 0.0);
        let normal = Vector3::new(1.0, 0.0, 0.0);
        let impulse = compute_impact_impulse(2.0, 2.0, rel_vel, normal, 0.0);
        let expected = -(1.0 + 0.0) * (-4.0) / (0.5 + 0.5);
        assert!((impulse - expected).abs() < 1e-6);
    }

    #[test]
    fn test_impact_impulse_perfectly_elastic() {
        let rel_vel = Vector3::new(-4.0, 0.0, 0.0);
        let normal = Vector3::new(1.0, 0.0, 0.0);
        let impulse = compute_impact_impulse(1.0, 1.0, rel_vel, normal, 1.0);
        let expected = -(1.0 + 1.0) * (-4.0) / (1.0 + 1.0);
        assert!((impulse - expected).abs() < 1e-6);
    }

    #[test]
    fn test_soft_impact_force() {
        let si = SoftImpact::new(1.0e7, 1.0e4);
        let force = si.force(0.001, -0.1);
        assert!((force - (1.0e7 * 0.001 + 1.0e4 * (-0.1))).abs() < 1e-6);
    }

    #[test]
    fn test_soft_impact_no_force_when_no_penetration() {
        let si = SoftImpact::default();
        assert_eq!(si.force(0.0, -0.1), 0.0);
        assert_eq!(si.force(-0.001, -0.1), 0.0);
    }
}
