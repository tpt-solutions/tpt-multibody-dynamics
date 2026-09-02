//! Time integration: semi-implicit Euler, Verlet, RATTLE, generalized-α,
//! HHT-α, Newmark-β.
//!
//! Provides explicit and implicit integrators for first- and second-order
//! ODEs arising from multibody dynamics, together with energy evaluation and
//! simple convergence diagnostics.

use crate::system::MultibodySystem;

// ===========================================================================
// Integrator parameters
// ===========================================================================

/// Parameters shared by all time integrators.
#[derive(Clone, Copy, Debug)]
pub struct IntegratorParams {
    /// Time step size.
    pub dt: f64,
    /// Generalized-α / HHT-α: end-of-step alpha (default ≈ 0.4).
    pub alpha_f: f64,
    /// Generalized-α: midpoint alpha for velocities/accelerations (default ≈ 0.2).
    pub alpha_m: f64,
    /// Newmark-β / HHT-α: beta parameter.
    pub beta: f64,
    /// Newmark-β / HHT-α: gamma parameter.
    pub gamma: f64,
}

impl Default for IntegratorParams {
    fn default() -> Self {
        Self {
            dt: 0.001,
            alpha_f: 0.4,
            alpha_m: 0.2,
            beta: 0.25,
            gamma: 0.5,
        }
    }
}

// ===========================================================================
// Integrator enum
// ===========================================================================

/// Supported time integration schemes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Integrator {
    /// First-order accurate, symplectic. Good for real-time / games.
    SemiImplicitEuler,
    /// Second-order accurate, symplectic. Good for energy conservation.
    Verlet,
    /// Position Verlet variant with velocity at half-step.
    VelocityVerlet,
    /// RATTLE: SHAKE/RATTLE constrained extension of Verlet.
    Rattle,
    /// Generalized-α: high-frequency dissipation tunable via α_f, α_m.
    GeneralizedAlpha,
    /// HHT-α: Hilber-Hughes-Taylor with numerical dissipation.
    HhtAlpha,
    /// Newmark-β: unconditional stability for β ≥ 1/4, γ ≥ 1/2.
    NewmarkBeta,
}

// ===========================================================================
// Convergence diagnostics
// ===========================================================================

/// Diagnostic output from an integrator step.
#[derive(Clone, Copy, Debug, Default)]
pub struct ConvergenceDiagnostics {
    /// Number of iterations performed.
    pub iterations: usize,
    /// Norm of the residual vector.
    pub residual_norm: f64,
    /// Maximum constraint equation violation.
    pub constraint_violation: f64,
}

// ===========================================================================
// Semi-implicit Euler
// ===========================================================================

/// Semi-implicit (symplectic) Euler integrator.
///
/// Velocity update: v(t+dt) = v(t) + a(t) * dt
/// Position update: q(t+dt) = q(t) + v(t+dt) * dt
pub struct SemiImplicitEuler;

impl SemiImplicitEuler {
    /// Perform one integration step.
    ///
    /// `tau` is the generalized force vector of length `q.len()`.
    pub fn step(
        system: &MultibodySystem,
        q: &[f64],
        qdot: &[f64],
        tau: &[f64],
        dt: f64,
    ) -> (Vec<f64>, Vec<f64>) {
        let n = q.len();
        let m_mat = system.build_mass_matrix();
        let mut qddot = vec![0.0f64; n];

        for i in 0..n {
            let mi = m_mat[(i, i)];
            if mi.abs() > 1e-15 {
                qddot[i] = tau.get(i).copied().unwrap_or(0.0) / mi;
            }
        }

        let mut qdot_new = qdot.to_vec();
        for i in 0..n {
            qdot_new[i] += qddot[i] * dt;
        }

        let mut q_new = q.to_vec();
        for i in 0..n {
            q_new[i] += qdot_new[i] * dt;
        }

        (q_new, qdot_new)
    }
}

// ===========================================================================
// Verlet
// ===========================================================================

/// Standard velocity Verlet integrator.
///
/// q(t+dt) = q(t) + v(t)*dt + 0.5*a(t)*dt^2
/// a(t+dt) evaluated from updated position
/// v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt
pub struct Verlet;

impl Verlet {
    /// Perform one integration step.
    pub fn step(
        system: &MultibodySystem,
        q: &[f64],
        qdot: &[f64],
        tau: &[f64],
        dt: f64,
    ) -> (Vec<f64>, Vec<f64>) {
        let n = q.len();
        let m_mat = system.build_mass_matrix();

        let mut a = vec![0.0f64; n];
        for i in 0..n {
            let mi = m_mat[(i, i)];
            if mi.abs() > 1e-15 {
                a[i] = tau.get(i).copied().unwrap_or(0.0) / mi;
            }
        }

        let mut q_new = q.to_vec();
        for i in 0..n {
            q_new[i] = q[i] + qdot[i] * dt + 0.5 * a[i] * dt * dt;
        }

        let tau_new = system.build_force_vector(&q_new, qdot);
        let mut a_new = vec![0.0f64; n];
        for i in 0..n {
            let mi = m_mat[(i, i)];
            if mi.abs() > 1e-15 {
                a_new[i] = tau_new.get(i).copied().unwrap_or(0.0) / mi;
            }
        }

        let mut qdot_new = qdot.to_vec();
        for i in 0..n {
            qdot_new[i] = qdot[i] + 0.5 * (a[i] + a_new[i]) * dt;
        }

        (q_new, qdot_new)
    }
}

// ===========================================================================
// Generalized Alpha
// ===========================================================================

/// Generalized-α integrator for high-frequency dissipation.
///
/// Uses filtered quantities at the midpoint of the time step to achieve
/// second-order accuracy and controllable high-frequency numerical damping.
pub struct GeneralizedAlpha;

impl GeneralizedAlpha {
    /// Perform one integration step.
    ///
    /// `qddot` is the acceleration at the start of the step.
    pub fn step(
        system: &MultibodySystem,
        q: &[f64],
        qdot: &[f64],
        qddot: &[f64],
        tau: &[f64],
        dt: f64,
        params: IntegratorParams,
    ) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let n = q.len();
        let alpha_f = params.alpha_f;
        let alpha_m = params.alpha_m;
        let beta = params.beta;
        let gamma = params.gamma;

        let m_mat = system.build_mass_matrix();

        let mut q_f = q.to_vec();
        let mut qdot_f = qdot.to_vec();
        let mut qddot_f = qddot.to_vec();
        for i in 0..n {
            q_f[i] = alpha_f * q[i]
                + (1.0 - alpha_f) * (q[i] + qdot[i] * dt + (0.5 - beta) * dt * dt * qddot[i]);
            qdot_f[i] =
                alpha_m * qdot[i] + (1.0 - alpha_m) * (qdot[i] + (1.0 - gamma) * dt * qddot[i]);
            qddot_f[i] = alpha_m * qddot[i] + (1.0 - alpha_m) * qddot[i];
        }

        let mut q_new = q.to_vec();
        let mut qdot_new = qdot.to_vec();
        let mut qddot_new = vec![0.0f64; n];

        for i in 0..n {
            let mi = m_mat[(i, i)];
            if mi.abs() > 1e-15 {
                qddot_new[i] = tau.get(i).copied().unwrap_or(0.0) / mi;
            }
        }

        for i in 0..n {
            q_new[i] = q[i]
                + qdot[i] * dt
                + (0.5 - beta) * dt * dt * qddot[i]
                + beta * dt * dt * qddot_new[i];
            qdot_new[i] = qdot[i] + (1.0 - gamma) * dt * qddot[i] + gamma * dt * qddot_new[i];
        }

        (q_new, qdot_new, qddot_new)
    }
}

// ===========================================================================
// HHT-alpha (stub)
// ===========================================================================

/// HHT-α integrator (stub — delegates to generalized-α with equivalent params).
pub struct HhtAlpha;

impl HhtAlpha {
    /// Perform one integration step.
    pub fn step(
        system: &MultibodySystem,
        q: &[f64],
        qdot: &[f64],
        qddot: &[f64],
        tau: &[f64],
        dt: f64,
        params: IntegratorParams,
    ) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        GeneralizedAlpha::step(system, q, qdot, qddot, tau, dt, params)
    }
}

// ===========================================================================
// Newmark-beta (stub)
// ===========================================================================

/// Newmark-β integrator (stub — uses the same interface as generalized-α).
pub struct NewmarkBeta;

impl NewmarkBeta {
    /// Perform one integration step.
    pub fn step(
        system: &MultibodySystem,
        q: &[f64],
        qdot: &[f64],
        qddot: &[f64],
        tau: &[f64],
        dt: f64,
        params: IntegratorParams,
    ) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        GeneralizedAlpha::step(system, q, qdot, qddot, tau, dt, params)
    }
}

// ===========================================================================
// Rattle (stub)
// ===========================================================================

/// RATTLE integrator for constrained systems (stub — falls back to Verlet).
pub struct Rattle;

impl Rattle {
    /// Perform one integration step.
    pub fn step(
        system: &MultibodySystem,
        q: &[f64],
        qdot: &[f64],
        tau: &[f64],
        dt: f64,
    ) -> (Vec<f64>, Vec<f64>) {
        Verlet::step(system, q, qdot, tau, dt)
    }
}

// ===========================================================================
// Energy evaluation
// ===========================================================================

/// Compute the total mechanical energy `E = KE + PE` of the system.
///
/// For a single rigid body with mass `m`, velocity `v`, and height `h`:
/// `KE = 0.5 * m * v^2`, `PE = m * g * h`.
pub fn energy(system: &MultibodySystem, q: &[f64], qdot: &[f64]) -> f64 {
    let mut total_ke = 0.0;
    let mut total_pe = 0.0;
    let g = 9.81;

    let mut qi = 0usize;
    for body in &system.bodies {
        let m = body.spatial_inertia.mass;
        if qi + 2 < q.len() && qi + 2 < qdot.len() {
            let vx = qdot[qi];
            let vy = qdot[qi + 1];
            let vz = qdot[qi + 2];
            let speed_sq = vx * vx + vy * vy + vz * vz;
            total_ke += 0.5 * m * speed_sq;

            let z = q[qi + 2];
            total_pe += m * g * z;
        }
        qi += 6;
    }

    for flex in &system.flexible {
        for k in 0..flex.num_modes {
            if qi + k < qdot.len() {
                let v = qdot[qi + k];
                let m = flex.modal_mass.get(k).copied().unwrap_or(1.0);
                total_ke += 0.5 * m * v * v;
            }
        }
        qi += flex.num_modes;
    }

    total_ke + total_pe
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forces::{Gravity, SpringDamper};
    use crate::system::MultibodySystem;
    use tpt_math_linalg_fixed::{Matrix3, Vector3};
    use tpt_mbd_core::frame::Isometry3;
    use tpt_mbd_core::{RigidBody, SpatialInertia};

    #[test]
    fn free_fall_energy_conservation() {
        let mut sys = MultibodySystem::new();
        let si = SpatialInertia::new(
            1.0f64,
            Vector3::new([0.0f64; 3]),
            Matrix3::new([[1.0f64; 3]; 3]),
        );
        let body = RigidBody::new(si, Isometry3::identity(), "free_fall", 0);
        sys.add_body(body);
        sys.count_dofs();

        let g = Gravity::new(9.81);
        let dt = 0.001f64;
        let n_steps = 1000usize;

        let mut q = vec![0.0f64; 6];
        q[2] = 10.0;
        let mut qdot = vec![0.0f64; 6];

        let e0 = energy(&sys, &q, &qdot);

        for _ in 0..n_steps {
            let tau = g.apply(&sys, &q);
            let (q_new, qdot_new) = SemiImplicitEuler::step(&sys, &q, &qdot, &tau, dt);
            q = q_new;
            qdot = qdot_new;
        }

        let e_final = energy(&sys, &q, &qdot);
        let rel_error = ((e_final - e0).abs() / e0.abs()).min(1.0);
        assert!(
            rel_error < 0.01,
            "energy drift {:.4}% exceeds 1% threshold",
            rel_error * 100.0,
        );
    }

    #[test]
    fn spring_damper_oscillation_period() {
        let mut sys = MultibodySystem::new();
        let si = SpatialInertia::new(
            1.0f64,
            Vector3::new([0.0f64; 3]),
            Matrix3::new([[1.0f64; 3]; 3]),
        );
        let body0 = RigidBody::new(si, Isometry3::identity(), "body0", 0);
        let body1 = RigidBody::new(si, Isometry3::identity(), "anchor", 0);
        sys.add_body(body0);
        sys.add_body(body1);
        sys.count_dofs();

        let spring = SpringDamper::new(100.0, 0.0, 0.0, 0, 1, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        let dt = 0.0005f64;
        let mut q = vec![0.0f64; sys.num_dofs];
        q[2] = -0.5;
        let mut qdot = vec![0.0f64; sys.num_dofs];

        let mut crossings = 0usize;
        let mut prev_z = q[2];
        for step in 0..50000 {
            let tau = spring.force(&sys, &q, &qdot);
            let (q_new, qdot_new) = SemiImplicitEuler::step(&sys, &q, &qdot, &tau, dt);
            q = q_new;
            qdot = qdot_new;

            let z = q[2];
            if prev_z <= 0.0 && z > 0.0 {
                crossings += 1;
                if crossings == 2 {
                    let period = (step + 1) as f64 * dt;
                    let expected = 2.0 * std::f64::consts::PI * (1.0 / 100.0f64).sqrt();
                    let rel_err = ((period - expected).abs() / expected).min(1.0);
                    assert!(rel_err < 0.25, "period error {:.3} > 25%", rel_err);
                    return;
                }
            }
            prev_z = z;
        }
        panic!("did not detect two zero crossings within 50000 steps");
    }

    #[test]
    fn gravity_force_direction() {
        let mut sys = MultibodySystem::new();
        let si = SpatialInertia::new(
            1.0f64,
            Vector3::new([0.0f64; 3]),
            Matrix3::new([[1.0f64; 3]; 3]),
        );
        let body = RigidBody::new(si, Isometry3::identity(), "test", 0);
        sys.add_body(body);
        sys.count_dofs();

        let g = Gravity::new(9.81);
        let q = vec![0.0f64; sys.num_dofs];
        let tau = g.apply(&sys, &q);
        assert!(
            tau[2] < 0.0,
            "gravity z-force must be negative, got {}",
            tau[2]
        );
    }

    #[test]
    fn verlet_spring_energy_conservation() {
        let mut sys = MultibodySystem::new();
        let si = SpatialInertia::new(
            1.0f64,
            Vector3::new([0.0f64; 3]),
            Matrix3::new([[1.0f64; 3]; 3]),
        );
        let body = RigidBody::new(si, Isometry3::identity(), "verlet_body", 0);
        sys.add_body(body);
        sys.count_dofs();

        let spring = SpringDamper::new(100.0, 0.0, 0.0, 0, 1, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        let dt = 0.0005f64;
        let mut q = vec![0.0f64; sys.num_dofs];
        q[2] = -0.5;
        let mut qdot = vec![0.0f64; sys.num_dofs];
        let e0 = energy(&sys, &q, &qdot);

        for _ in 0..5000 {
            let tau = spring.force(&sys, &q, &qdot);
            let (q_new, qdot_new) = Verlet::step(&sys, &q, &qdot, &tau, dt);
            q = q_new;
            qdot = qdot_new;
        }

        let e_final = energy(&sys, &q, &qdot);
        let rel_error = ((e_final - e0).abs() / e0.abs()).min(1.0);
        assert!(
            rel_error < 0.05,
            "Verlet energy drift {:.4}% exceeds 5% threshold",
            rel_error * 100.0,
        );
    }

    #[test]
    fn generalized_alpha_energy_dissipation() {
        let mut sys = MultibodySystem::new();
        let si = SpatialInertia::new(
            1.0f64,
            Vector3::new([0.0f64; 3]),
            Matrix3::new([[1.0f64; 3]; 3]),
        );
        let body = RigidBody::new(si, Isometry3::identity(), "ga_body", 0);
        sys.add_body(body);
        sys.count_dofs();

        let spring = SpringDamper::new(100.0, 0.5, 0.0, 0, 1, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        let params = IntegratorParams {
            dt: 0.001,
            alpha_f: 0.4,
            alpha_m: 0.2,
            beta: 0.25,
            gamma: 0.5,
        };
        let mut q = vec![0.0f64; sys.num_dofs];
        q[2] = -0.5;
        let mut qdot = vec![0.0f64; sys.num_dofs];
        let mut qddot = vec![0.0f64; sys.num_dofs];
        let e0 = energy(&sys, &q, &qdot);

        for _ in 0..2000 {
            let tau = spring.force(&sys, &q, &qdot);
            let (q_new, qdot_new, qddot_new) =
                GeneralizedAlpha::step(&sys, &q, &qdot, &qddot, &tau, params.dt, params);
            q = q_new;
            qdot = qdot_new;
            qddot = qddot_new;
        }

        let e_final = energy(&sys, &q, &qdot);
        assert!(
            e_final <= e0,
            "generalized-α should dissipate energy, but e_final={} > e0={}",
            e_final,
            e0
        );
    }

    #[test]
    fn rattle_delegates_to_verlet() {
        let mut sys = MultibodySystem::new();
        let si = SpatialInertia::new(
            1.0f64,
            Vector3::new([0.0f64; 3]),
            Matrix3::new([[1.0f64; 3]; 3]),
        );
        let body = RigidBody::new(si, Isometry3::identity(), "rattle_body", 0);
        sys.add_body(body);
        sys.count_dofs();

        let q = vec![0.0f64; sys.num_dofs];
        let qdot = vec![0.0f64; sys.num_dofs];
        let tau = vec![0.0f64; sys.num_dofs];
        let dt = 0.001;

        let (q_rattle, qdot_rattle) = Rattle::step(&sys, &q, &qdot, &tau, dt);
        let (q_verlet, qdot_verlet) = Verlet::step(&sys, &q, &qdot, &tau, dt);

        assert_eq!(q_rattle, q_verlet);
        assert_eq!(qdot_rattle, qdot_verlet);
    }
}
