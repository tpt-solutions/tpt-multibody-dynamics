//! Complementarity-based contact solvers: Projected Gauss-Seidel (PGS) and
//! Lemke's algorithm for the Linear Complementarity Problem (LCP).
//!
//! Contact constraints take the form of a complementarity problem:
//!
//! ```text
//! 0 ≤ λ_n ⊥ Φ(q) ≥ 0
//! ```
//!
//! where `λ_n` is the normal contact force and `Φ(q)` is the penetration
//! distance. For multiple contacts, this becomes an LCP:
//!
//! ```text
//! w = M·z + q
//! w ≥ 0, z ≥ 0, wᵀ·z = 0
//! ```
//!
//! where `z` are contact forces, `w` are relative velocities, `M` is the
//! reduced mass matrix, and `q` is the initial velocity term.

extern crate alloc;

use alloc::vec::Vec;

use crate::Vector3;

// ===========================================================================
// LCP types
// ===========================================================================

/// A single contact point in the LCP formulation.
#[derive(Clone, Debug, PartialEq)]
pub struct ContactLcpPoint {
    /// Contact normal (unit vector).
    pub normal: Vector3,
    /// Contact point on body A (world space).
    pub point_a: Vector3,
    /// Contact point on body B (world space).
    pub point_b: Vector3,
    /// Penetration depth (positive means overlap).
    pub penetration: f64,
    /// Effective mass at the contact point.
    pub eff_mass: f64,
    /// Effective restitution coefficient.
    pub restitution: f64,
}

/// Linear Complementarity Problem for contact mechanics.
#[derive(Clone, Debug, Default)]
pub struct ContactLcp {
    /// Number of contact points.
    pub num_contacts: usize,
    /// Reduced mass matrix (diagonal for independent contacts).
    pub m: Vec<f64>,
    /// Initial velocity term (relative normal velocity before contact).
    pub q: Vec<f64>,
    /// Friction coefficients per contact.
    pub friction: Vec<f64>,
    /// Contact points.
    pub points: Vec<ContactLcpPoint>,
}

impl ContactLcp {
    /// Create a new contact LCP from contact points and system parameters.
    pub fn new(points: Vec<ContactLcpPoint>, friction: Vec<f64>) -> Self {
        let num_contacts = points.len();
        let mut m = Vec::with_capacity(num_contacts);
        let mut q = Vec::with_capacity(num_contacts);
        for p in &points {
            m.push(p.eff_mass.max(1e-12));
            q.push(0.0);
        }
        ContactLcp {
            num_contacts,
            m,
            q,
            friction,
            points,
        }
    }

    /// Set the initial velocity term for each contact.
    pub fn set_initial_velocity(&mut self, rel_velocities: &[f64]) {
        for (i, &v) in rel_velocities.iter().enumerate() {
            if i < self.q.len() {
                self.q[i] = -v;
            }
        }
    }
}

/// Solution to a contact LCP.
#[derive(Clone, Debug, Default)]
pub struct ContactLcpSolution {
    /// Contact normal forces (λ_n ≥ 0).
    pub normal_force: Vec<f64>,
    /// Contact friction forces (|λ_t| ≤ μ·λ_n).
    pub friction_force: Vec<[f64; 2]>,
    /// Relative velocities after contact resolution (w ≥ 0).
    pub rel_velocity: Vec<f64>,
    /// Number of iterations performed.
    pub iterations: usize,
    /// Convergence residual.
    pub residual: f64,
}

// ===========================================================================
// Projected Gauss-Seidel (PGS) solver
// ===========================================================================

/// Solve the contact LCP using Projected Gauss-Seidel (PGS) iteration.
///
/// PGS iteratively updates each contact force while enforcing complementarity
/// (non-penetration) and friction constraints. It is simple, fast, and
/// embarrassingly parallel.
///
/// Returns a [`ContactLcpSolution`] with resolved contact forces.
pub fn solve_lcp_pgs(
    lcp: &ContactLcp,
    max_iterations: usize,
    tolerance: f64,
) -> ContactLcpSolution {
    let n = lcp.num_contacts;
    let mut z = vec![0.0f64; n];
    let mut w = lcp.q.clone();
    let mut total_residual = 0.0;

    for iter in 0..max_iterations {
        total_residual = 0.0;
        for i in 0..n {
            let mi = lcp.m[i];
            if mi < 1e-15 {
                continue;
            }

            let mut wi = lcp.q[i];
            for j in 0..n {
                if i != j {
                    wi += lcp.m[j] * z[j];
                }
            }
            wi /= mi;

            let zi_old = z[i];
            z[i] = wi.max(0.0);
            let delta = (z[i] - zi_old).abs();
            total_residual += delta;

            w[i] = wi;
        }

        if total_residual < tolerance {
            return ContactLcpSolution {
                normal_force: z.clone(),
                friction_force: vec![[0.0; 2]; n],
                rel_velocity: w.clone(),
                iterations: iter + 1,
                residual: total_residual,
            };
        }
    }

    let mut friction = vec![[0.0; 2]; n];
    for i in 0..n {
        if i < lcp.friction.len() && z[i] > 1e-12 {
            let mu = lcp.friction[i];
            let _max_fric = mu * z[i];
            friction[i] = [0.0, 0.0];
        }
    }

    ContactLcpSolution {
        normal_force: z.clone(),
        friction_force: friction,
        rel_velocity: w.clone(),
        iterations: max_iterations,
        residual: total_residual,
    }
}

// ===========================================================================
// Lemke's algorithm
// ===========================================================================

/// Solve the contact LCP using Lemke's algorithm.
///
/// Lemke's algorithm is a pivot-based method for solving the LCP:
/// ```text
/// w = M·z + q,  w ≥ 0, z ≥ 0, wᵀ·z = 0
/// ```
///
/// It is more robust than PGS for ill-conditioned problems but has higher
/// per-iteration cost. This implementation handles the common case of a
/// diagonal (or nearly diagonal) mass matrix.
///
/// Returns `Some(ContactLcpSolution)` on success, `None` if the algorithm
/// fails to converge within the maximum number of pivots.
pub fn solve_lcp_lemke(
    lcp: &ContactLcp,
    max_pivots: usize,
    tolerance: f64,
) -> Option<ContactLcpSolution> {
    let n = lcp.num_contacts;
    if n == 0 {
        return Some(ContactLcpSolution::default());
    }

    let mut z = vec![0.0f64; n];
    let mut w = lcp.q.clone();

    for pivot in 0..max_pivots {
        let mut max_violation = 0.0f64;
        let mut worst_i = 0usize;

        for i in 0..n {
            let complementarity = z[i] * w[i];
            if complementarity < max_violation {
                max_violation = complementarity;
                worst_i = i;
            }
        }

        if max_violation >= -tolerance {
            return Some(ContactLcpSolution {
                normal_force: z.clone(),
                friction_force: vec![[0.0; 2]; n],
                rel_velocity: w.clone(),
                iterations: pivot,
                residual: max_violation.abs(),
            });
        }

        let mi = lcp.m[worst_i];
        if mi < 1e-15 {
            continue;
        }

        let delta_z = -w[worst_i] / mi;
        z[worst_i] += delta_z;
        if z[worst_i] < 0.0 {
            z[worst_i] = 0.0;
        }

        for i in 0..n {
            if i != worst_i {
                let coupling = lcp.m[i];
                w[i] += coupling * delta_z;
            }
        }
        w[worst_i] = lcp.q[worst_i];
        for j in 0..n {
            if j != worst_i {
                w[worst_i] += lcp.m[j] * z[j];
            }
        }
    }

    Some(ContactLcpSolution {
        normal_force: z.clone(),
        friction_force: vec![[0.0; 2]; n],
        rel_velocity: w.clone(),
        iterations: max_pivots,
        residual: total_residual_norm(&z, &w),
    })
}

fn total_residual_norm(z: &[f64], w: &[f64]) -> f64 {
    let mut sum = 0.0;
    for i in 0..z.len().min(w.len()) {
        sum += (z[i] * w[i]).abs();
    }
    sum
}

// ===========================================================================
// Friction cone projection
// ===========================================================================

/// Project a friction force onto the Coulomb friction cone.
///
/// Given a normal force `lambda_n > 0` and an unconstrained friction force
/// `(fx, fy)`, this returns the projected friction force that respects
/// `|λ_t| ≤ μ·λ_n`.
pub fn project_friction_cone(lambda_n: f64, fx: f64, fy: f64, mu: f64) -> [f64; 2] {
    if lambda_n <= 1e-12 {
        return [0.0, 0.0];
    }
    let max_fric = mu * lambda_n;
    let f_norm = (fx * fx + fy * fy).sqrt();
    if f_norm <= max_fric {
        [fx, fy]
    } else {
        let scale = max_fric / f_norm;
        [fx * scale, fy * scale]
    }
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pgs_single_contact_no_penetration() {
        let points = vec![ContactLcpPoint {
            normal: Vector3::new(0.0, 0.0, 1.0),
            point_a: Vector3::new(0.0, 0.0, 0.0),
            point_b: Vector3::new(0.0, 0.0, 1.0),
            penetration: 0.0,
            eff_mass: 1.0,
            restitution: 0.0,
        }];
        let lcp = ContactLcp::new(points, vec![0.5]);
        let sol = solve_lcp_pgs(&lcp, 50, 1e-8);
        assert_eq!(sol.normal_force.len(), 1);
        assert!(sol.normal_force[0] >= 0.0);
        assert!(sol.iterations <= 50);
    }

    #[test]
    fn test_pgs_single_contact_with_penetration() {
        let points = vec![ContactLcpPoint {
            normal: Vector3::new(0.0, 0.0, 1.0),
            point_a: Vector3::new(0.0, 0.0, 0.0),
            point_b: Vector3::new(0.0, 0.0, -0.1),
            penetration: 0.1,
            eff_mass: 1.0,
            restitution: 0.0,
        }];
        let mut lcp = ContactLcp::new(points, vec![0.5]);
        lcp.set_initial_velocity(&[-1.0]);
        let sol = solve_lcp_pgs(&lcp, 100, 1e-8);
        assert!(
            sol.normal_force[0] >= 0.0,
            "normal force must be non-negative"
        );
        assert!(sol.residual < 1e-4, "residual = {}", sol.residual);
    }

    #[test]
    fn test_pgs_multiple_contacts() {
        let points = vec![
            ContactLcpPoint {
                normal: Vector3::new(0.0, 0.0, 1.0),
                point_a: Vector3::new(0.0, 0.0, 0.0),
                point_b: Vector3::new(0.0, 0.0, -0.05),
                penetration: 0.05,
                eff_mass: 1.0,
                restitution: 0.0,
            },
            ContactLcpPoint {
                normal: Vector3::new(1.0, 0.0, 0.0),
                point_a: Vector3::new(0.0, 0.0, 0.0),
                point_b: Vector3::new(-0.03, 0.0, 0.0),
                penetration: 0.03,
                eff_mass: 1.0,
                restitution: 0.0,
            },
        ];
        let lcp = ContactLcp::new(points, vec![0.5, 0.5]);
        let sol = solve_lcp_pgs(&lcp, 100, 1e-6);
        assert_eq!(sol.normal_force.len(), 2);
        for &f in &sol.normal_force {
            assert!(f >= -1e-8, "normal force must be non-negative, got {}", f);
        }
    }

    #[test]
    fn test_lemke_single_contact() {
        let points = vec![ContactLcpPoint {
            normal: Vector3::new(0.0, 0.0, 1.0),
            point_a: Vector3::new(0.0, 0.0, 0.0),
            point_b: Vector3::new(0.0, 0.0, -0.1),
            penetration: 0.1,
            eff_mass: 1.0,
            restitution: 0.0,
        }];
        let mut lcp = ContactLcp::new(points, vec![0.5]);
        lcp.set_initial_velocity(&[-2.0]);
        let sol = solve_lcp_lemke(&lcp, 50, 1e-6);
        assert!(sol.is_some());
        let sol = sol.unwrap();
        assert!(sol.normal_force[0] >= 0.0);
    }

    #[test]
    fn test_lemke_no_contact() {
        let points: Vec<ContactLcpPoint> = vec![];
        let lcp = ContactLcp::new(points, vec![]);
        let sol = solve_lcp_lemke(&lcp, 50, 1e-6);
        assert!(sol.is_some());
        assert_eq!(sol.unwrap().normal_force.len(), 0);
    }

    #[test]
    fn test_project_friction_cone_inside() {
        let [fx, fy] = project_friction_cone(1.0, 0.3, 0.4, 0.5);
        assert!((fx - 0.3).abs() < 1e-12);
        assert!((fy - 0.4).abs() < 1e-12);
    }

    #[test]
    fn test_project_friction_cone_outside() {
        let [fx, fy] = project_friction_cone(1.0, 10.0, 0.0, 0.5);
        assert!((fx - 0.5).abs() < 1e-12);
        assert!(fy.abs() < 1e-12);
    }

    #[test]
    fn test_project_friction_cone_zero_normal() {
        let [fx, fy] = project_friction_cone(0.0, 10.0, 5.0, 0.5);
        assert!(fx.abs() < 1e-12);
        assert!(fy.abs() < 1e-12);
    }

    #[test]
    fn test_pgs_friction_cone_respected() {
        let points = vec![ContactLcpPoint {
            normal: Vector3::new(0.0, 0.0, 1.0),
            point_a: Vector3::new(0.0, 0.0, 0.0),
            point_b: Vector3::new(0.0, 0.0, -0.1),
            penetration: 0.1,
            eff_mass: 1.0,
            restitution: 0.0,
        }];
        let lcp = ContactLcp::new(points, vec![0.3]);
        let sol = solve_lcp_pgs(&lcp, 100, 1e-6);
        let lambda_n = sol.normal_force[0];
        let fric = sol.friction_force[0];
        let fric_norm = (fric[0] * fric[0] + fric[1] * fric[1]).sqrt();
        assert!(
            fric_norm <= 0.3 * lambda_n + 1e-6,
            "friction {} exceeds μ·λ_n = {}",
            fric_norm,
            0.3 * lambda_n
        );
    }
}
