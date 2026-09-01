//! Contact force models: Hertzian, Hunt-Crossley, penalty, and augmented Lagrangian.

/// Parameters for nonlinear contact force models.
#[derive(Clone, Debug, PartialEq)]
pub struct ContactParams {
    /// Hertzian stiffness coefficient.
    pub stiffness: f64,
    /// Nonlinear exponent (1.5 for spheres).
    pub exponent: f64,
    /// Damping coefficient (Hunt-Crossley).
    pub damping: f64,
}

impl Default for ContactParams {
    fn default() -> Self {
        Self {
            stiffness: 1.0e5,
            exponent: 1.5,
            damping: 0.0,
        }
    }
}

/// Hertzian contact force: F = k * δ^n
#[derive(Clone, Debug, PartialEq)]
pub struct HertzianContact {
    /// Contact parameters.
    pub params: ContactParams,
}

impl HertzianContact {
    /// Create a new Hertzian contact model.
    pub fn new(params: ContactParams) -> Self {
        Self { params }
    }

    /// Compute the Hertzian contact force for a given penetration.
    pub fn force(&self, penetration: f64, _penetration_rate: f64) -> f64 {
        if penetration <= 0.0 {
            return 0.0;
        }
        self.params.stiffness * penetration.powf(self.params.exponent)
    }
}

/// Hunt-Crossley contact force: F = k*δ^n + c*δ^n * δ̇
#[derive(Clone, Debug, PartialEq)]
pub struct HuntCrossleyContact {
    /// Contact parameters.
    pub params: ContactParams,
}

impl HuntCrossleyContact {
    /// Create a new Hunt-Crossley contact model.
    pub fn new(params: ContactParams) -> Self {
        Self { params }
    }

    /// Compute the Hunt-Crossley contact force.
    pub fn force(&self, penetration: f64, penetration_rate: f64) -> f64 {
        if penetration <= 0.0 {
            return 0.0;
        }
        let delta_n = penetration.powf(self.params.exponent);
        self.params.stiffness * delta_n + self.params.damping * delta_n * penetration_rate
    }
}

/// Linear penalty contact force: F_n = k_p * δ - c_p * δ̇
#[derive(Clone, Debug, PartialEq)]
pub struct PenaltyContact {
    /// Penalty stiffness coefficient.
    pub stiffness: f64,
    /// Damping coefficient.
    pub damping: f64,
}

impl Default for PenaltyContact {
    fn default() -> Self {
        Self {
            stiffness: 1.0e5,
            damping: 500.0,
        }
    }
}

impl PenaltyContact {
    /// Create a new penalty contact model.
    pub fn new(stiffness: f64, damping: f64) -> Self {
        Self { stiffness, damping }
    }

    /// Compute the penalty contact force.
    pub fn force(&self, penetration: f64, penetration_rate: f64) -> f64 {
        self.stiffness * penetration - self.damping * penetration_rate
    }
}

/// Augmented Lagrangian contact with iterative constraint enforcement.
#[derive(Clone, Debug, PartialEq)]
pub struct AugmentedLagrangianContact {
    /// Penalty parameter.
    pub penalty: f64,
    /// Convergence tolerance for the iterative solve.
    pub tolerance: f64,
    /// Maximum number of iterations.
    pub max_iterations: usize,
}

impl Default for AugmentedLagrangianContact {
    fn default() -> Self {
        Self {
            penalty: 1.0e5,
            tolerance: 1e-6,
            max_iterations: 20,
        }
    }
}

impl AugmentedLagrangianContact {
    /// Create a new augmented Lagrangian contact model.
    pub fn new(penalty: f64, tolerance: f64, max_iterations: usize) -> Self {
        Self {
            penalty,
            tolerance,
            max_iterations,
        }
    }

    /// Compute the constraint force via augmented Lagrangian iteration.
    pub fn compute_force(&self, penetration: f64, penetration_rate: f64, mass_eff: f64) -> f64 {
        if penetration <= 0.0 || penetration_rate >= 0.0 {
            return 0.0;
        }
        let mut lambda = 0.0;
        for _ in 0..self.max_iterations {
            let residual = self.penalty * penetration + lambda / mass_eff + penetration_rate;
            let denominator = self.penalty + 1.0 / mass_eff;
            let correction = residual / denominator;
            lambda += correction;
            if correction.abs() < self.tolerance {
                break;
            }
        }
        lambda.max(0.0)
    }
}

/// Computes a contact force using the given parameters.
///
/// Uses the Hunt-Crossley formula when damping is non-zero,
/// otherwise falls back to Hertzian.
pub fn compute_contact_force(
    params: &ContactParams,
    penetration: f64,
    penetration_rate: f64,
) -> f64 {
    if penetration <= 0.0 {
        return 0.0;
    }
    let delta_n = penetration.powf(params.exponent);
    params.stiffness * delta_n + params.damping * delta_n * penetration_rate
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hertzian_contact() {
        let hz = HertzianContact::new(ContactParams {
            stiffness: 1.0e5,
            exponent: 1.5,
            damping: 0.0,
        });
        let force = hz.force(0.01, 0.0);
        let expected = 1.0e5 * 0.01_f64.powf(1.5);
        assert!((force - expected).abs() < 1e-6);
    }

    #[test]
    fn test_hunt_crossley_contact() {
        let hc = HuntCrossleyContact::new(ContactParams {
            stiffness: 1.0e5,
            exponent: 1.5,
            damping: 100.0,
        });
        let force = hc.force(0.01, -0.1);
        let delta_n = 0.01_f64.powf(1.5);
        let expected = 1.0e5 * delta_n + 100.0 * delta_n * (-0.1);
        assert!((force - expected).abs() < 1e-6);
    }

    #[test]
    fn test_penalty_contact() {
        let pc = PenaltyContact::new(1.0e5, 500.0);
        let force = pc.force(0.01, -0.1);
        assert!((force - (1.0e5 * 0.01 - 500.0 * (-0.1))).abs() < 1e-6);
    }

    #[test]
    fn test_no_force_when_no_penetration_hertzian() {
        let hz = HertzianContact::new(ContactParams::default());
        assert_eq!(hz.force(0.0, 0.0), 0.0);
        assert_eq!(hz.force(-0.01, 0.0), 0.0);
    }

    #[test]
    fn test_no_force_when_no_penetration_hunt_crossley() {
        let hc = HuntCrossleyContact::new(ContactParams {
            stiffness: 1.0e5,
            exponent: 1.5,
            damping: 100.0,
        });
        assert_eq!(hc.force(0.0, -0.1), 0.0);
        assert_eq!(hc.force(-0.01, -0.1), 0.0);
    }

    #[test]
    fn test_augmented_lagrangian_positive_force() {
        let al = AugmentedLagrangianContact::new(1.0e5, 1e-6, 20);
        let force = al.compute_force(0.01, -0.1, 1.0);
        assert!(force > 0.0);
    }

    #[test]
    fn test_augmented_lagrangian_no_force_when_separating() {
        let al = AugmentedLagrangianContact::default();
        let force = al.compute_force(0.01, 0.1, 1.0);
        assert_eq!(force, 0.0);
    }

    #[test]
    fn test_compute_contact_force_default() {
        let params = ContactParams::default();
        let force = compute_contact_force(&params, 0.01, 0.0);
        assert!((force - 1.0e5 * 0.01_f64.powf(1.5)).abs() < 1e-6);
    }

    #[test]
    fn test_compute_contact_force_with_damping() {
        let params = ContactParams {
            stiffness: 1.0e5,
            exponent: 1.5,
            damping: 100.0,
        };
        let force = compute_contact_force(&params, 0.01, -0.1);
        let delta_n = 0.01_f64.powf(1.5);
        let expected = 1.0e5 * delta_n + 100.0 * delta_n * (-0.1);
        assert!((force - expected).abs() < 1e-6);
    }

    #[test]
    fn test_contact_params_default_exponent() {
        let params = ContactParams::default();
        assert_eq!(params.exponent, 1.5);
    }
}
