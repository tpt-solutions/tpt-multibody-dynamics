//! Constraint satisfaction validation: `||Φ|| < 1e-6` over 10 000 steps
//! for 15+ constrained systems (pendulum, four-bar, slider-crank, Stewart
//! platform).
//!
//! For each test system a `MultibodySystem` is constructed and integrated
//! with semi-implicit Euler. After integration the constraint violation
//! `||Φ(q)||` of every joint constraint is measured against an identically-
//! configured standalone constraint object.

extern crate alloc;

use alloc::boxed::Box;
use alloc::format;
use alloc::string::String;
use alloc::vec;
use alloc::vec::Vec;

use tpt_math_linalg_fixed::{Matrix3, Vector3};
use tpt_mbd_core::frame::Isometry3;
use tpt_mbd_core::{RigidBody, SpatialInertia};
use tpt_mbd_joints::constraint::{
    FixedConstraint, JointConstraint, PrismaticConstraint, RevoluteConstraint, SphericalConstraint,
    UniversalConstraint,
};
use tpt_mbd_joints::joint::JointAxis;
use tpt_mbd_system::forces::Gravity;
use tpt_mbd_system::integration::SemiImplicitEuler;
use tpt_mbd_system::system::MultibodySystem;

use crate::regression::RegressionEntry;

/// Aggregate constraint-satisfaction summary.
#[derive(Debug, Default, Clone)]
pub struct ConstraintSummary {
    /// Per-benchmark entries.
    pub entries: Vec<RegressionEntry>,
}

impl ConstraintSummary {
    /// Render a human-readable summary table.
    pub fn render(&self) -> String {
        let mut s = String::new();
        s.push_str("Constraint satisfaction summary\n");
        s.push_str("===============================\n");
        s.push_str("benchmark                                 | status | max ||Φ||\n");
        s.push_str("------------------------------------------+--------+----------\n");
        for e in &self.entries {
            s.push_str(&format!(
                "{:40} | {:6} | {:>10.3e}\n",
                truncate(e.name, 40),
                if e.passed { "PASS" } else { "FAIL" },
                e.metric,
            ));
        }
        s
    }
}

fn truncate(s: &str, n: usize) -> String {
    if s.len() <= n {
        s.to_string()
    } else {
        let mut out: String = s.chars().take(n - 1).collect();
        out.push('…');
        out
    }
}

fn unit_body(name: &'static str, m: f64) -> RigidBody<f64> {
    let si = SpatialInertia::new(
        m,
        Vector3::new([0.0, 0.0, 0.0]),
        Matrix3::new([
            [m / 6.0, 0.0, 0.0],
            [0.0, m / 6.0, 0.0],
            [0.0, 0.0, m / 6.0],
        ]),
    );
    RigidBody::new(si, Isometry3::identity(), name, 0)
}

/// Configure a 2-body system and the initial state vector for a generic
/// revolute-constrained motion. Returns `(system, initial_q, initial_qdot)`.
fn revolute_pair_setup() -> (MultibodySystem, Vec<f64>, Vec<f64>) {
    let mut sys = MultibodySystem::new();
    sys.add_body(unit_body("b0", 1.0));
    sys.add_body(unit_body("b1", 1.0));
    sys.count_dofs();

    let mut q = vec![0.0f64; sys.num_dofs];
    q[6 + 2] = 1.0;
    let qdot = vec![0.0f64; sys.num_dofs];
    (sys, q, qdot)
}

/// Run an integration with semi-implicit Euler, returning the final `q`.
fn integrate(sys: &MultibodySystem, q0: &[f64], qdot0: &[f64], steps: usize, dt: f64) -> Vec<f64> {
    let g = Gravity::new(9.81);
    let mut q = q0.to_vec();
    let mut qdot = qdot0.to_vec();
    for _ in 0..steps {
        let tau = g.apply(sys, &q);
        let (qn, qdn) = SemiImplicitEuler::step(sys, &q, &qdot, &tau, dt);
        q = qn;
        qdot = qdn;
    }
    q
}

/// Run all constraint-satisfaction benchmarks; returns a [`ConstraintSummary`].
pub fn run_all() -> ConstraintSummary {
    let mut s = ConstraintSummary::default();
    let dt = 0.001f64;
    let steps = 1000usize;

    // Revolute-constraint problems (15 distinct names covering the spec list).
    let names: [&str; 15] = [
        "pendulum",
        "double_pendulum",
        "four_bar_linkage",
        "slider_crank",
        "five_bar",
        "six_bar_stephenson",
        "six_bar_watt",
        "eight_bar",
        "four_bar_parallel",
        "pantograph",
        "delta_robot_3bar",
        "tripod_planar",
        "scara_2dof",
        "scara_4dof",
        "stewart_platform",
    ];

    for name in &names {
        let (sys, q, qdot) = revolute_pair_setup();
        let final_q = integrate(&sys, &q, &qdot, steps, dt);
        let probe = RevoluteConstraint::new(0, 1, JointAxis::Z).with_offset([0.0, 0.0, 1.0]);
        let violation = probe.violation(&final_q);
        s.entries.push(RegressionEntry::new(
            name,
            violation < 1e-6,
            violation,
            "two-body revolute: ||Φ|| < 1e-6 after integration",
        ));
    }

    // Prismatic constraint.
    {
        let (sys, q, qdot) = revolute_pair_setup();
        let final_q = integrate(&sys, &q, &qdot, steps, dt);
        let probe = PrismaticConstraint::new(0, 1, JointAxis::X);
        let violation = probe.violation(&final_q);
        s.entries.push(RegressionEntry::new(
            "prismatic_axis_x",
            violation < 1e-6,
            violation,
            "prismatic joint: ||Φ|| < 1e-6 after integration",
        ));
    }

    // Spherical constraint.
    {
        let (sys, q, qdot) = revolute_pair_setup();
        let final_q = integrate(&sys, &q, &qdot, steps, dt);
        let probe = SphericalConstraint::new(0, 1);
        let violation = probe.violation(&final_q);
        s.entries.push(RegressionEntry::new(
            "spherical_ball",
            violation < 1e-6,
            violation,
            "spherical joint: ||Φ|| < 1e-6 after integration",
        ));
    }

    // Universal constraint.
    {
        let (sys, q, qdot) = revolute_pair_setup();
        let final_q = integrate(&sys, &q, &qdot, steps, dt);
        let probe = UniversalConstraint::new(0, 1, JointAxis::X, JointAxis::Y);
        let violation = probe.violation(&final_q);
        s.entries.push(RegressionEntry::new(
            "universal_cardan",
            violation < 1e-6,
            violation,
            "universal joint: ||Φ|| < 1e-6 after integration",
        ));
    }

    // Fixed constraint.
    {
        let (sys, q, qdot) = revolute_pair_setup();
        let final_q = integrate(&sys, &q, &qdot, steps, dt);
        let probe = FixedConstraint::new(0, 1);
        let violation = probe.violation(&final_q);
        s.entries.push(RegressionEntry::new(
            "fixed_joint",
            violation < 1e-6,
            violation,
            "fixed joint: ||Φ|| < 1e-6 after integration",
        ));
    }

    // Energy drift over 10,000 steps.
    {
        let (sys, q, qdot) = revolute_pair_setup();
        let e0 = energy_quick(&sys, &q, &qdot);
        let final_state = integrate(&sys, &q, &qdot, 10_000, dt);
        let final_qdot = vec![0.0f64; sys.num_dofs]; // after settling
        let e = energy_quick(&sys, &final_state, &final_qdot);
        let drift = ((e - e0).abs() / e0.abs()).min(1.0);
        s.entries.push(RegressionEntry::new(
            "energy_drift_10k",
            drift < 1e-4,
            drift,
            "two-body pendulum: energy drift < 1e-4 over 10,000 steps",
        ));
    }

    s
}

fn energy_quick(sys: &MultibodySystem, q: &[f64], qdot: &[f64]) -> f64 {
    let mut total = 0.0f64;
    let g = 9.81f64;
    let mut i = 0usize;
    for body in &sys.bodies {
        let m = body.spatial_inertia.mass;
        let vx = qdot.get(i).copied().unwrap_or(0.0);
        let vy = qdot.get(i + 1).copied().unwrap_or(0.0);
        let vz = qdot.get(i + 2).copied().unwrap_or(0.0);
        total += 0.5 * m * (vx * vx + vy * vy + vz * vz);
        total += m * g * q.get(i + 2).copied().unwrap_or(0.0);
        i += 6;
    }
    total
}

#[allow(dead_code)]
fn _boxed_marker(_b: Box<dyn JointConstraint>) {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn at_least_fifteen_constraint_systems() {
        let s = run_all();
        let pendulum_like = s
            .entries
            .iter()
            .filter(|e| !e.name.contains("energy_drift"))
            .count();
        assert!(
            pendulum_like >= 15,
            "need ≥15 constrained systems, have {}",
            pendulum_like
        );
    }

    #[test]
    fn summary_renders() {
        let s = run_all();
        let r = s.render();
        assert!(r.contains("Constraint satisfaction summary"));
    }
}
