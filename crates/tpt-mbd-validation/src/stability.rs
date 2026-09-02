//! Numerical-stability validation: stiff systems stable over 100 000+ time
//! steps.
//!
//! A stiff system is approximated with a high-stiffness spring (k = 1e6)
//! acting on a unit mass with semi-implicit Euler at a step size Δt = 1e-5.
//! Semi-implicit Euler is unconditionally stable for linear systems; this
//! test asserts that 100 000 steps do not produce non-finite values.

extern crate alloc;

use alloc::format;
use alloc::string::String;
use alloc::vec;
use alloc::vec::Vec;

use tpt_math_linalg_fixed::{Matrix3, Vector3};
use tpt_mbd_core::frame::Isometry3;
use tpt_mbd_core::{RigidBody, SpatialInertia};
use tpt_mbd_system::forces::SpringDamper;
use tpt_mbd_system::integration::SemiImplicitEuler;
use tpt_mbd_system::system::MultibodySystem;

use crate::regression::RegressionEntry;

/// Aggregate stability summary.
#[derive(Debug, Default, Clone)]
pub struct StabilitySummary {
    /// Per-test entries.
    pub entries: Vec<RegressionEntry>,
}

impl StabilitySummary {
    /// Render a human-readable summary table.
    pub fn render(&self) -> String {
        let mut s = String::new();
        s.push_str("Stability summary\n");
        s.push_str("=================\n");
        s.push_str("benchmark                                 | status | max |q|\n");
        s.push_str("------------------------------------------+--------+-------\n");
        for e in &self.entries {
            s.push_str(&format!(
                "{:40} | {:6} | {:>7.3e}\n",
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

/// Stiff oscillator: k = 1e6, m = 1.0 → ω ≈ 1 000 rad/s, period ≈ 6.3 ms.
/// With dt = 1e-5 s, the natural period spans ~630 steps — enough to
/// validate that the integrator survives 100 000 steps without blow-up.
fn stiff_oscillator(n_steps: usize, k: f64, dt: f64) -> f64 {
    let mut sys = MultibodySystem::new();
    sys.add_body(unit_body("m", 1.0));
    sys.add_body(unit_body("anchor", 1.0e9));
    sys.count_dofs();

    let spring = SpringDamper::new(k, 0.1, 0.0, 0, 1, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]);

    let mut q = vec![0.0f64; sys.num_dofs];
    q[2] = 0.001;
    let mut qdot = vec![0.0f64; sys.num_dofs];

    let mut max_abs_q = 0.0f64;
    for _ in 0..n_steps {
        let tau = spring.force(&sys, &q, &qdot);
        let (qn, qdn) = SemiImplicitEuler::step(&sys, &q, &qdot, &tau, dt);
        q = qn;
        qdot = qdn;
        let abs = q.iter().map(|x| x.abs()).fold(0.0f64, f64::max);
        if abs > max_abs_q {
            max_abs_q = abs;
        }
        if !q.iter().all(|x| x.is_finite()) {
            return f64::INFINITY;
        }
    }
    max_abs_q
}

/// Run all stability benchmarks.
pub fn run_all() -> StabilitySummary {
    let mut s = StabilitySummary::default();

    let cases: [(&str, f64, f64, usize); 5] = [
        ("stiff_k1e4_100k_dt1e5", 1.0e4, 1.0e-5, 100_000),
        ("stiff_k1e5_100k_dt1e5", 1.0e5, 1.0e-5, 100_000),
        ("stiff_k1e6_100k_dt1e5", 1.0e6, 1.0e-5, 100_000),
        ("stiff_k1e6_50k_dt1e6", 1.0e6, 1.0e-6, 50_000),
        ("stiff_k1e7_20k_dt1e6", 1.0e7, 1.0e-6, 20_000),
    ];

    for (name, k, dt, n) in &cases {
        let max_abs_q = stiff_oscillator(*n, *k, *dt);
        // Acceptance: state remains finite and bounded.
        let passed = max_abs_q.is_finite() && max_abs_q < 1.0e3;
        s.entries.push(RegressionEntry::new(
            name,
            passed,
            max_abs_q,
            "stiff oscillator survives N steps without blow-up",
        ));
    }

    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_stability_cases_pass() {
        let s = run_all();
        for e in &s.entries {
            assert!(e.passed, "{} failed (max |q|={})", e.name, e.metric);
        }
    }

    #[test]
    fn summary_renders() {
        let s = run_all();
        let r = s.render();
        assert!(r.contains("Stability summary"));
    }
}
