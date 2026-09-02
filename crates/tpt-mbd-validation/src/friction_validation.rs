//! Friction validation: force vs. analytical solutions for block-on-plane,
//! rolling wheel, brake pad. Spec: force < 10% error, correct stick-slip
//! transition.
//!
//! Uses the Coulomb friction kernel from `tpt-mbd-contact` and compares to
//! closed-form expectations: `F = μ * N` for sliding, and `F → 0` as the
//! tangential velocity vanishes (smooth regularisation).

extern crate alloc;

use alloc::format;
use alloc::string::String;
use alloc::vec::Vec;

use tpt_mbd_contact::friction::{
    compute_friction_force, AnisotropicFriction, CoulombFrictionParams, FrictionModel,
    StribeckFriction,
};
use tpt_mbd_contact::Vector3;

use crate::regression::RegressionEntry;

/// Aggregate friction-validation summary.
#[derive(Debug, Default, Clone)]
pub struct FrictionSummary {
    /// Per-benchmark entries.
    pub entries: Vec<RegressionEntry>,
}

impl FrictionSummary {
    /// Render a human-readable summary table.
    pub fn render(&self) -> String {
        let mut s = String::new();
        s.push_str("Friction validation summary\n");
        s.push_str("===========================\n");
        s.push_str("benchmark                                 | status | rel. error\n");
        s.push_str("------------------------------------------+--------+-----------\n");
        for e in &self.entries {
            s.push_str(&format!(
                "{:40} | {:6} | {:>11.4}\n",
                truncate(e.name, 40),
                if e.passed { "PASS" } else { "FAIL" },
                e.metric * 100.0,
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

/// Run all friction-validation benchmarks.
pub fn run_all() -> FrictionSummary {
    let mut s = FrictionSummary::default();

    // ---- Block-on-plane: sliding velocity sweep, Coulomb model.
    let params = CoulombFrictionParams {
        static_coeff: 0.5,
        kinetic_coeff: 0.4,
        stribeck_velocity: 0.01,
        viscous_coeff: 0.0,
    };

    let sweep: [(&str, f64); 10] = [
        ("block_plane_v0p01", 0.01),
        ("block_plane_v0p05", 0.05),
        ("block_plane_v0p1", 0.1),
        ("block_plane_v0p5", 0.5),
        ("block_plane_v1", 1.0),
        ("block_plane_v2", 2.0),
        ("block_plane_v5", 5.0),
        ("block_plane_v10", 10.0),
        ("block_plane_v50", 50.0),
        ("block_plane_v100", 100.0),
    ];

    let normal = 100.0f64;
    for (name, v) in &sweep {
        let vel = Vector3::new(*v, 0.0, 0.0);
        let f = compute_friction_force(normal, vel, &params, FrictionModel::Coulomb);
        // For pure Coulomb, F = μ_k * N, direction opposite to velocity.
        let expected_mag = params.kinetic_coeff * normal;
        let measured_mag = f.norm();
        let rel = ((measured_mag - expected_mag).abs() / expected_mag).min(1.0);
        s.entries.push(RegressionEntry::new(
            name,
            rel < 0.10,
            rel,
            "Coulomb friction force vs. analytical μ_k * N, <10% error",
        ));
    }

    // ---- Stick-slip transition: zero velocity → zero force.
    {
        let vel = Vector3::new(0.0, 0.0, 0.0);
        let f = compute_friction_force(normal, vel, &params, FrictionModel::Coulomb);
        s.entries.push(RegressionEntry::new(
            "stick_zero_velocity",
            f == Vector3::new(0.0, 0.0, 0.0),
            if f == Vector3::new(0.0, 0.0, 0.0) {
                0.0
            } else {
                1.0
            },
            "stick regime: zero velocity → zero friction",
        ));
    }

    // ---- Rolling wheel: tangential velocity at very low values, Stribeck
    //      coefficient must be ≥ kinetic coefficient and ≤ static coefficient.
    {
        let coeff_v0 = StribeckFriction::coefficient(0.0, &params);
        let coeff_v_inf = StribeckFriction::coefficient(1.0e6, &params);
        let within_bounds = coeff_v0 >= params.kinetic_coeff - 1e-9
            && coeff_v0 <= params.static_coeff + 1e-9
            && (coeff_v_inf - params.kinetic_coeff).abs() < 1e-3;
        s.entries.push(RegressionEntry::new(
            "stribeck_bounds",
            within_bounds,
            (coeff_v0 - coeff_v_inf).abs(),
            "Stribeck coefficient bounded by [μ_k, μ_s]",
        ));
    }

    // ---- Brake pad: viscous friction contribution.
    let brake_params = CoulombFrictionParams {
        static_coeff: 0.6,
        kinetic_coeff: 0.5,
        stribeck_velocity: 0.01,
        viscous_coeff: 0.05,
    };
    {
        let vel = Vector3::new(2.0, 0.0, 0.0);
        let f = compute_friction_force(200.0, vel, &brake_params, FrictionModel::Coulomb);
        let expected_mag = brake_params.kinetic_coeff * 200.0 + brake_params.viscous_coeff * 2.0;
        let rel = ((f.norm() - expected_mag).abs() / expected_mag).min(1.0);
        s.entries.push(RegressionEntry::new(
            "brake_pad_viscous",
            rel < 0.10,
            rel,
            "Coulomb+viscous: F = μ_k N + c v, <10% error",
        ));
    }

    // ---- Anisotropic friction: per-axis coefficient.
    {
        let af = AnisotropicFriction::new(0.5, 0.8, 0.3);
        let vel_x = Vector3::new(1.0, 0.0, 0.0);
        let f_x = af.force(10.0, vel_x);
        let expected_x = 0.5 * 10.0;
        let rel_x = ((f_x.norm() - expected_x).abs() / expected_x).min(1.0);
        s.entries.push(RegressionEntry::new(
            "anisotropic_x",
            rel_x < 0.10,
            rel_x,
            "anisotropic friction x-axis: F = μ_x * N",
        ));

        let vel_y = Vector3::new(0.0, 1.0, 0.0);
        let f_y = af.force(10.0, vel_y);
        let expected_y = 0.8 * 10.0;
        let rel_y = ((f_y.norm() - expected_y).abs() / expected_y).min(1.0);
        s.entries.push(RegressionEntry::new(
            "anisotropic_y",
            rel_y < 0.10,
            rel_y,
            "anisotropic friction y-axis: F = μ_y * N",
        ));

        let vel_z = Vector3::new(0.0, 0.0, 1.0);
        let f_z = af.force(10.0, vel_z);
        let expected_z = 0.3 * 10.0;
        let rel_z = ((f_z.norm() - expected_z).abs() / expected_z).min(1.0);
        s.entries.push(RegressionEntry::new(
            "anisotropic_z",
            rel_z < 0.10,
            rel_z,
            "anisotropic friction z-axis: F = μ_z * N",
        ));
    }

    // ---- Smooth Coulomb transition: velocity just above stribeck_velocity
    //      should yield a coefficient between μ_k and μ_s.
    {
        let sm = compute_friction_force(
            10.0,
            Vector3::new(0.05, 0.0, 0.0),
            &params,
            FrictionModel::SmoothCoulomb,
        );
        let measured_coeff = sm.norm() / 10.0;
        let within = measured_coeff >= params.kinetic_coeff - 1e-6
            && measured_coeff <= params.static_coeff + 1e-6;
        s.entries.push(RegressionEntry::new(
            "smooth_coulomb_transition",
            within,
            (measured_coeff - params.kinetic_coeff).abs(),
            "Smooth Coulomb: μ_k ≤ μ ≤ μ_s in transition",
        ));
    }

    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_friction_cases_pass() {
        let s = run_all();
        for e in &s.entries {
            assert!(e.passed, "{} failed (metric {})", e.name, e.metric);
        }
    }

    #[test]
    fn summary_renders() {
        let s = run_all();
        let r = s.render();
        assert!(r.contains("Friction validation summary"));
    }
}
