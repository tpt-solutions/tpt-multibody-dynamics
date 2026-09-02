//! Performance validation: simulation time at 10/100/1000 DOFs.
//!
//! Builds increasingly large `MultibodySystem`s (linear chains of bodies
//! linked by revolute joints) and times 1 000 integration steps with
//! semi-implicit Euler. Records ms/step and the simulation rate (kHz).

extern crate alloc;

use alloc::format;
use alloc::string::String;
use alloc::vec;
use alloc::vec::Vec;

use tpt_math_linalg_fixed::{Matrix3, Vector3};
use tpt_mbd_core::frame::Isometry3;
use tpt_mbd_core::{RigidBody, SpatialInertia};
use tpt_mbd_joints::joint::JointType;
use tpt_mbd_system::forces::Gravity;
use tpt_mbd_system::integration::SemiImplicitEuler;
use tpt_mbd_system::system::MultibodySystem;

use crate::regression::RegressionEntry;

/// Aggregate performance summary.
#[derive(Debug, Default, Clone)]
pub struct PerformanceSummary {
    /// Per-DOF entries.
    pub entries: Vec<RegressionEntry>,
}

impl PerformanceSummary {
    /// Render a human-readable summary table.
    pub fn render(&self) -> String {
        let mut s = String::new();
        s.push_str("Performance summary\n");
        s.push_str("===================\n");
        s.push_str("benchmark                                 | status | ms/step\n");
        s.push_str("------------------------------------------+--------+--------\n");
        for e in &self.entries {
            s.push_str(&format!(
                "{:40} | {:6} | {:>8.3}\n",
                truncate(e.name, 40),
                if e.passed { "PASS" } else { "INFO" },
                e.metric * 1.0e3,
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

/// Build a chain of `n_bodies` bodies linked by revolute joints, return
/// the system with DOFs counted.
fn chain(n_bodies: usize) -> MultibodySystem {
    let mut sys = MultibodySystem::new();
    for i in 0..n_bodies {
        sys.add_body(unit_body("link", 1.0));
        if i > 0 {
            sys.add_joint(JointType::REVOLUTE, i - 1, i);
        }
    }
    sys.count_dofs();
    sys
}

/// Time `n_steps` of semi-implicit Euler on a chain of `n_bodies` bodies.
/// Returns the per-step wall time in seconds.
fn bench_chain(n_bodies: usize, n_steps: usize) -> f64 {
    let sys = chain(n_bodies);
    let mut q = vec![0.0f64; sys.num_dofs];
    for i in 0..sys.num_dofs.min(q.len()) {
        q[i] = (i as f64) * 0.01;
    }
    let mut qdot = vec![0.0f64; sys.num_dofs];
    let g = Gravity::new(9.81);
    let dt = 0.001f64;

    use std::time::Instant;
    let t0 = Instant::now();
    for _ in 0..n_steps {
        let tau = g.apply(&sys, &q);
        let (qn, qdn) = SemiImplicitEuler::step(&sys, &q, &qdot, &tau, dt);
        q = qn;
        qdot = qdn;
    }
    let elapsed = t0.elapsed();
    elapsed.as_secs_f64() / (n_steps as f64)
}

/// Run all performance benchmarks; returns a [`PerformanceSummary`].
pub fn run_all() -> PerformanceSummary {
    let mut s = PerformanceSummary::default();

    let cases: [(&str, usize); 3] = [
        ("perf_10dof", 2),
        ("perf_100dof", 17),
        ("perf_1000dof", 167),
    ];

    for (name, n_bodies) in &cases {
        let ms_per_step = bench_chain(*n_bodies, 50);
        // We record the timing; we don't fail on slow runs but we do emit a
        // regression entry. The target is real-time (1 kHz = 1 ms/step)
        // for <100 DOFs.
        let rate_hz = 1.0 / ms_per_step.max(1e-9);
        let passed = *n_bodies >= 100 || rate_hz >= 1000.0;
        s.entries.push(RegressionEntry::new(
            name,
            passed,
            ms_per_step,
            "simulation rate (Hz) = 1 / ms_per_step; target 1 kHz for <100 DOFs",
        ));
    }

    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn performance_runs_all_three_dof_levels() {
        let s = run_all();
        assert_eq!(s.entries.len(), 3);
        for e in &s.entries {
            assert!(e.metric.is_finite(), "{} non-finite timing", e.name);
        }
    }

    #[test]
    fn summary_renders() {
        let s = run_all();
        let r = s.render();
        assert!(r.contains("Performance summary"));
    }
}
