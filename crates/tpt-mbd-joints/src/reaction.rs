//! Reaction force computation from Lagrange multipliers.
//!
//! Provides:
//! - [`ReactionForce`] — force and torque vectors
//! - [`compute_reactions`] — extract constraint forces from Lagrange multipliers

#![allow(missing_docs)]

use num_traits::Float;

/// Reaction force and torque at a joint.
///
/// Represents the constraint forces and torques exerted by a joint on
/// the connected bodies, computed from Lagrange multipliers.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReactionForce {
    pub force: [f64; 3],
    pub torque: [f64; 3],
}

impl ReactionForce {
    /// Create a new reaction force.
    pub fn new(force: [f64; 3], torque: [f64; 3]) -> Self {
        Self { force, torque }
    }

    /// Create a zero reaction force.
    pub fn zero() -> Self {
        Self {
            force: [0.0; 3],
            torque: [0.0; 3],
        }
    }

    /// Compute the magnitude of the force vector.
    pub fn force_magnitude(&self) -> f64 {
        let sum: f64 = self.force.iter().map(|x| x * x).sum();
        Float::sqrt(sum)
    }

    /// Compute the magnitude of the torque vector.
    pub fn torque_magnitude(&self) -> f64 {
        let sum: f64 = self.torque.iter().map(|x| x * x).sum();
        Float::sqrt(sum)
    }
}

impl Default for ReactionForce {
    fn default() -> Self {
        Self::zero()
    }
}

/// Compute reaction forces from Lagrange multipliers and a constraint Jacobian.
///
/// Given Lagrange multipliers `lambda` and the constraint Jacobian `jacobian`,
/// the reaction forces are computed as:
///
/// `f = Jᵀ · λ`
///
/// where `Jᵀ` is the transpose of the constraint Jacobian and `λ` are the
/// Lagrange multipliers for each constraint equation.
///
/// # Arguments
///
/// * `lambda` — Lagrange multipliers, one per constraint equation
/// * `jacobian` — constraint Jacobian rows, each row `[f64; 3]` representing
///   derivatives wrt one body's position DOFs
///
/// # Returns
///
/// A [`ReactionForce`] with the computed force and torque vectors.
pub fn compute_reactions(lambda: &[f64], jacobian: &[[f64; 3]]) -> ReactionForce {
    let mut force = [0.0_f64; 3];
    let mut torque = [0.0_f64; 3];

    for (row, &lam) in jacobian.iter().zip(lambda.iter()) {
        force[0] += row[0] * lam;
        force[1] += row[1] * lam;
        force[2] += row[2] * lam;
    }

    for (row, &lam) in jacobian.iter().zip(lambda.iter()) {
        torque[0] += row[0] * lam;
        torque[1] += row[1] * lam;
        torque[2] += row[2] * lam;
    }

    ReactionForce::new(force, torque)
}

/// Populate reaction forces from a constraint Jacobian and multipliers.
///
/// This is a convenience wrapper around [`compute_reactions`] that also
/// handles multiple bodies by accepting separate Jacobian blocks.
pub fn populate_reactions(lambda: &[f64], jacobian_blocks: &[&[[f64; 3]]]) -> ReactionForce {
    let mut force = [0.0_f64; 3];
    let mut torque = [0.0_f64; 3];
    let mut idx = 0;

    for block in jacobian_blocks {
        for row in *block {
            if idx < lambda.len() {
                let lam = lambda[idx];
                force[0] += row[0] * lam;
                force[1] += row[1] * lam;
                force[2] += row[2] * lam;
                torque[0] += row[0] * lam;
                torque[1] += row[1] * lam;
                torque[2] += row[2] * lam;
                idx += 1;
            }
        }
    }

    ReactionForce::new(force, torque)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compute_reactions_zero_lambda() {
        let jac = &[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let lambda = &[0.0, 0.0];
        let r = compute_reactions(lambda, jac);
        assert_eq!(r.force, [0.0; 3]);
    }

    #[test]
    fn test_compute_reactions_unit_lambda() {
        let jac = &[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let lambda = &[1.0, 2.0];
        let r = compute_reactions(lambda, jac);
        assert_eq!(r.force, [1.0, 2.0, 0.0]);
    }

    #[test]
    fn test_reaction_force_magnitude() {
        let r = ReactionForce::new([3.0, 4.0, 0.0], [0.0, 0.0, 0.0]);
        assert_eq!(r.force_magnitude(), 5.0);
    }

    #[test]
    fn test_reaction_zero() {
        let r = ReactionForce::zero();
        assert_eq!(r.force_magnitude(), 0.0);
        assert_eq!(r.torque_magnitude(), 0.0);
    }
}
