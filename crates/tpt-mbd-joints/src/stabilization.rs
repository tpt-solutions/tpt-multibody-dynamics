//! Constraint stabilization methods for multibody systems.
//!
//! Implements:
//! - Baumgarte stabilization (Φ̈ + 2αΦ̇ + β²Φ = 0)
//! - Coordinate partitioning (independent/dependent DOF split)
//! - Augmented Lagrangian stabilization

#![allow(missing_docs)]

/// Baumgarte stabilization parameters.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BaumgarteParams {
    pub alpha: f64,
    pub beta: f64,
}

impl Default for BaumgarteParams {
    fn default() -> Self {
        Self {
            alpha: 20.0,
            beta: 100.0,
        }
    }
}

impl BaumgarteParams {
    pub fn new(alpha: f64, beta: f64) -> Self {
        Self { alpha, beta }
    }

    /// Auto-tune parameters based on system natural frequency and time step.
    pub fn auto_tune(natural_freq: f64, _dt: f64) -> Self {
        let w = natural_freq.max(1e-6);
        let alpha = 2.0 * w;
        let beta = w * w;
        Self { alpha, beta }
    }
}

/// Stabilized constraint acceleration using Baumgarte's method.
///
/// Instead of enforcing Φ̈ = 0, we enforce:
/// Φ̈ + 2αΦ̇ + β²Φ = 0
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BaumgarteStabilization {
    pub params: BaumgarteParams,
}

impl BaumgarteStabilization {
    pub fn new(params: BaumgarteParams) -> Self {
        Self { params }
    }

    /// Compute the stabilized constraint acceleration.
    ///
    /// Given constraint violation Φ and constraint velocity Φ̇, returns the
    /// Baumgarte-stabilized acceleration Φ̈_stab.
    pub fn stabilized_acceleration(&self, phi: f64, phi_dot: f64) -> f64 {
        let BaumgarteParams { alpha, beta } = self.params;
        -2.0 * alpha * phi_dot - beta * beta * phi
    }

    /// Compute Baumgarte stabilization force for a unit-mass system.
    pub fn stabilization_force(&self, phi: f64, phi_dot: f64, dt: f64) -> f64 {
        let acc = self.stabilized_acceleration(phi, phi_dot);
        acc * dt * dt
    }
}

impl Default for BaumgarteStabilization {
    fn default() -> Self {
        Self::new(BaumgarteParams::default())
    }
}

/// Coordinate partitioning: split DOFs into independent and dependent sets.
#[derive(Clone, Debug, PartialEq)]
pub struct CoordinatePartitioning {
    pub independent: alloc::vec::Vec<usize>,
    pub dependent: alloc::vec::Vec<usize>,
}

impl CoordinatePartitioning {
    pub fn new(total_dofs: usize, dependent: alloc::vec::Vec<usize>) -> Self {
        let mut independent = alloc::vec::Vec::new();
        for i in 0..total_dofs {
            if !dependent.contains(&i) {
                independent.push(i);
            }
        }
        Self {
            independent,
            dependent,
        }
    }

    pub fn num_independent(&self) -> usize {
        self.independent.len()
    }
    pub fn num_dependent(&self) -> usize {
        self.dependent.len()
    }
    pub fn total(&self) -> usize {
        self.independent.len() + self.dependent.len()
    }
}

/// Augmented Lagrangian stabilization: combines penalty and multiplier methods.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AugmentedLagrangian {
    pub penalty: f64,
    pub tolerance: f64,
    pub max_iterations: usize,
}

impl Default for AugmentedLagrangian {
    fn default() -> Self {
        Self {
            penalty: 1e6,
            tolerance: 1e-6,
            max_iterations: 50,
        }
    }
}

impl AugmentedLagrangian {
    pub fn new(penalty: f64, tolerance: f64, max_iterations: usize) -> Self {
        Self {
            penalty,
            tolerance,
            max_iterations,
        }
    }

    /// Stabilized constraint value: Φ + (1/k)·λ
    pub fn stabilized_constraint(&self, phi: f64, lambda: f64) -> f64 {
        phi + lambda / self.penalty
    }
}
