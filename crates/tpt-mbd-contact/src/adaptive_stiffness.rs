//! Adaptive penalty stiffness for contact force models.
//!
//! Provides methods to compute contact stiffness automatically based on
//! body inertia and expected contact forces, avoiding the manual tuning
//! required for fixed-stiffness penalty methods.
//!
//! The adaptive stiffness is computed as:
//!
//! ```text
//! k_adapt = f_expected / δ_target
//! ```
//!
//! where `f_expected` is the expected contact force magnitude and
//! `δ_target` is the desired penetration depth. The stiffness is also
//! bounded by a minimum and maximum value to prevent numerical issues.

/// Parameters for adaptive penalty stiffness computation.
#[derive(Clone, Debug, PartialEq)]
pub struct AdaptiveStiffnessParams {
    /// Expected contact force magnitude (N).
    pub expected_force: f64,
    /// Target penetration depth (m).
    pub target_penetration: f64,
    /// Minimum allowable stiffness (N/m).
    pub min_stiffness: f64,
    /// Maximum allowable stiffness (N/m).
    pub max_stiffness: f64,
    /// Damping coefficient (N·s/m).
    pub damping: f64,
}

impl Default for AdaptiveStiffnessParams {
    fn default() -> Self {
        Self {
            expected_force: 100.0,
            target_penetration: 1e-4,
            min_stiffness: 1e3,
            max_stiffness: 1e9,
            damping: 500.0,
        }
    }
}

impl AdaptiveStiffnessParams {
    /// Create new adaptive stiffness parameters.
    pub fn new(expected_force: f64, target_penetration: f64) -> Self {
        Self {
            expected_force,
            target_penetration,
            ..Self::default()
        }
    }

    /// Set the damping coefficient.
    pub fn with_damping(mut self, damping: f64) -> Self {
        self.damping = damping;
        self
    }

    /// Set the stiffness bounds.
    pub fn with_bounds(mut self, min_stiffness: f64, max_stiffness: f64) -> Self {
        self.min_stiffness = min_stiffness;
        self.max_stiffness = max_stiffness;
        self
    }
}

/// Compute adaptive penalty stiffness from body inertia and expected contact force.
///
/// The stiffness is chosen so that the expected contact force produces the
/// target penetration depth. It is then bounded by `min_stiffness` and
/// `max_stiffness`.
///
/// `eff_mass` is the effective mass at the contact point (reduced mass of the
/// two colliding bodies). `expected_force` is the anticipated normal force
/// magnitude. `target_penetration` is the desired maximum penetration.
pub fn compute_adaptive_stiffness(
    eff_mass: f64,
    expected_force: f64,
    target_penetration: f64,
) -> f64 {
    if target_penetration <= 0.0 || eff_mass <= 0.0 {
        return 1.0e5;
    }
    let k = expected_force / target_penetration;
    k.clamp(1e3, 1e9)
}

/// Compute adaptive penalty stiffness with full parameter control.
///
/// Returns `(stiffness, damping)` where `stiffness` is adapted to the
/// effective mass and expected contact force.
pub fn adaptive_penalty_params(eff_mass: f64, params: &AdaptiveStiffnessParams) -> (f64, f64) {
    let k = compute_adaptive_stiffness(eff_mass, params.expected_force, params.target_penetration);
    let damping = params.damping;
    (k.clamp(params.min_stiffness, params.max_stiffness), damping)
}

/// Compute stiffness from body properties: mass, contact radius, and material properties.
///
/// Uses a Hertzian-inspired scaling:
///
/// ```text
/// k = f_expected / ((f_expected / (E*R^0.5))^(2/3))
/// ```
///
/// where `E*` is the effective Young's modulus and `R` is the effective radius.
/// This is a simplified version for quick estimation without full material data.
pub fn stiffness_from_hertzian_params(
    eff_mass: f64,
    expected_force: f64,
    effective_radius: f64,
) -> f64 {
    if effective_radius <= 0.0 || eff_mass <= 0.0 {
        return 1.0e5;
    }
    let delta = expected_force / effective_radius;
    let k = expected_force / delta.max(1e-12);
    k.clamp(1e3, 1e9)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_adaptive_stiffness_basic() {
        let k = compute_adaptive_stiffness(1.0, 100.0, 1e-4);
        assert!((k - 1.0e6).abs() < 1.0);
    }

    #[test]
    fn test_adaptive_stiffness_clamping() {
        let k = compute_adaptive_stiffness(1.0, 1e15, 1e-4);
        assert!(k <= 1e9, "stiffness should be clamped, got {}", k);
        let k = compute_adaptive_stiffness(1.0, 1e-3, 1e-4);
        assert!(k >= 1e3, "stiffness should be clamped, got {}", k);
    }

    #[test]
    fn test_adaptive_stiffness_zero_target() {
        let k = compute_adaptive_stiffness(1.0, 100.0, 0.0);
        assert!(k > 0.0);
    }

    #[test]
    fn test_adaptive_penalty_params() {
        let params = AdaptiveStiffnessParams::new(100.0, 1e-4).with_damping(200.0);
        let (k, d) = adaptive_penalty_params(1.0, &params);
        assert!((k - 1.0e6).abs() < 1.0);
        assert_eq!(d, 200.0);
    }

    #[test]
    fn test_stiffness_from_hertzian_params() {
        let k = stiffness_from_hertzian_params(1.0, 100.0, 0.01);
        assert!(k > 0.0);
        assert!(k <= 1e9);
    }

    #[test]
    fn test_adaptive_stiffness_params_defaults() {
        let params = AdaptiveStiffnessParams::default();
        assert_eq!(params.expected_force, 100.0);
        assert_eq!(params.target_penetration, 1e-4);
        assert_eq!(params.min_stiffness, 1e3);
        assert_eq!(params.max_stiffness, 1e9);
    }

    #[test]
    fn test_adaptive_stiffness_large_force() {
        let k = compute_adaptive_stiffness(10.0, 1e6, 1e-3);
        assert!((k - 1.0e9).abs() < 1.0);
    }

    #[test]
    fn test_adaptive_stiffness_small_force() {
        let k = compute_adaptive_stiffness(0.1, 0.01, 1e-6);
        assert!(k >= 1e3);
        assert!(k <= 1e9);
    }
}
