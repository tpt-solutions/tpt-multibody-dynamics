//! External force application: gravity, spring-damper, prescribed motion, bushings.
//!
//! Provides force elements that can be assembled into the generalized force
//! vector `τ(q, q̇)` for a `MultibodySystem`.

use crate::system::MultibodySystem;
pub use crate::Vector;

use std::sync::Arc;

// ===========================================================================
// Gravity
// ===========================================================================

/// Uniform gravitational field.
#[derive(Clone, Copy, Debug, Default)]
pub struct Gravity {
    /// Gravitational acceleration (m/s², default 9.81).
    pub g: f64,
}

impl Gravity {
    /// Create a gravity model with the given acceleration magnitude.
    pub fn new(g: f64) -> Self {
        Self { g }
    }

    /// Apply gravity to every rigid body in the system.
    ///
    /// Returns the generalized force vector `τ` of length `num_dofs`.
    /// Gravity acts downward (−z): `f = [0, 0, −m·g]` per body.
    pub fn apply(&self, system: &MultibodySystem, _q: &[f64]) -> Vec<f64> {
        let n = system.num_dofs.max(1);
        let mut tau = vec![0.0f64; n];

        let mut offset = 0usize;
        for body in &system.bodies {
            let m = body.spatial_inertia.mass;
            let fz = -m * self.g;
            if offset + 2 < n {
                tau[offset + 2] = fz;
            }
            offset += 6;
        }

        for flex in &system.flexible {
            for k in 0..flex.num_modes {
                let idx = offset + k;
                if idx < n {
                    tau[idx] = 0.0;
                }
            }
            offset += flex.num_modes;
        }

        tau
    }
}

// ===========================================================================
// Spring-Damper
// ===========================================================================

/// Translational spring-damper element between two attachment points.
#[derive(Clone, Debug)]
pub struct SpringDamper {
    /// Spring stiffness (N/m).
    pub k: f64,
    /// Damping coefficient (N·s/m).
    pub c: f64,
    /// Rest length (m).
    pub rest_length: f64,
    /// Index of the first body.
    pub body_i: usize,
    /// Index of the second body.
    pub body_j: usize,
    /// Attachment point on body `i` expressed in body-local frame.
    pub attachment_point_i: [f64; 3],
    /// Attachment point on body `j` expressed in body-local frame.
    pub attachment_point_j: [f64; 3],
}

impl SpringDamper {
    /// Build a new spring-damper element.
    pub fn new(
        k: f64,
        c: f64,
        rest_length: f64,
        body_i: usize,
        body_j: usize,
        attachment_point_i: [f64; 3],
        attachment_point_j: [f64; 3],
    ) -> Self {
        Self {
            k,
            c,
            rest_length,
            body_i,
            body_j,
            attachment_point_i,
            attachment_point_j,
        }
    }

    /// Evaluate spring-damper force for the current configuration.
    ///
    /// Returns the generalized force vector contribution (length `num_dofs`).
    pub fn force(&self, system: &MultibodySystem, q: &[f64], qdot: &[f64]) -> Vec<f64> {
        let n = system.num_dofs.max(1);
        let mut tau = vec![0.0f64; n];

        let bi = self.body_i.min(system.bodies.len().saturating_sub(1));
        let bj = self.body_j.min(system.bodies.len().saturating_sub(1));

        let pi = self.attachment_point_i;
        let pj = self.attachment_point_j;

        let qi_off = bi * 6;
        let qj_off = bj * 6;

        let dx = q.get(qi_off).copied().unwrap_or(0.0) - q.get(qj_off).copied().unwrap_or(0.0)
            + pi[0]
            - pj[0];
        let dy = q.get(qi_off + 1).copied().unwrap_or(0.0)
            - q.get(qj_off + 1).copied().unwrap_or(0.0)
            + pi[1]
            - pj[1];
        let dz = q.get(qi_off + 2).copied().unwrap_or(0.0)
            - q.get(qj_off + 2).copied().unwrap_or(0.0)
            + pi[2]
            - pj[2];

        let dist = (dx * dx + dy * dy + dz * dz).sqrt().max(1e-12);
        let direction = 1.0 / dist;

        let vi_x =
            qdot.get(qi_off).copied().unwrap_or(0.0) - qdot.get(qj_off).copied().unwrap_or(0.0);
        let vi_y = qdot.get(qi_off + 1).copied().unwrap_or(0.0)
            - qdot.get(qj_off + 1).copied().unwrap_or(0.0);
        let vi_z = qdot.get(qi_off + 2).copied().unwrap_or(0.0)
            - qdot.get(qj_off + 2).copied().unwrap_or(0.0);

        let rel_vel = (vi_x * dx + vi_y * dy + vi_z * dz) * direction;
        let stretch = self.rest_length - dist;

        let force_mag = self.k * stretch + self.c * rel_vel;
        let fx = force_mag * dx * direction;
        let fy = force_mag * dy * direction;
        let fz = force_mag * dz * direction;

        if qi_off + 2 < n {
            tau[qi_off] += fx;
            tau[qi_off + 1] += fy;
            tau[qi_off + 2] += fz;
        }
        if qj_off + 2 < n {
            tau[qj_off] -= fx;
            tau[qj_off + 1] -= fy;
            tau[qj_off + 2] -= fz;
        }

        tau
    }
}

// ===========================================================================
// Prescribed motion
// ===========================================================================

/// Kinematic driver that prescribes the trajectory of a single body.
#[derive(Clone)]
pub struct PrescribedMotion {
    /// Index of the body being driven.
    pub body_index: usize,
    /// Prescribed position as a function of time: `(x, y, z)`.
    pub prescribed_q: Arc<dyn Fn(f64) -> [f64; 3]>,
    /// Prescribed velocity as a function of time: `(vx, vy, vz)`.
    pub prescribed_qdot: Arc<dyn Fn(f64) -> [f64; 3]>,
}

impl PrescribedMotion {
    /// Create a new prescribed motion driver.
    pub fn new<F, G>(body_index: usize, prescribed_q: F, prescribed_qdot: G) -> Self
    where
        F: Fn(f64) -> [f64; 3] + 'static,
        G: Fn(f64) -> [f64; 3] + 'static,
    {
        Self {
            body_index,
            prescribed_q: Arc::new(prescribed_q),
            prescribed_qdot: Arc::new(prescribed_qdot),
        }
    }

    /// Evaluate the driver at time `t`, returning `(q, qdot)`.
    pub fn driver(&self, t: f64) -> ([f64; 3], [f64; 3]) {
        ((self.prescribed_q)(t), (self.prescribed_qdot)(t))
    }
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::system::MultibodySystem;
    use tpt_math_linalg_fixed::{Matrix3, Vector3};
    use tpt_mbd_core::frame::Isometry3;
    use tpt_mbd_core::{RigidBody, SpatialInertia};

    #[test]
    fn gravity_force_direction_check() {
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
        assert_eq!(tau[0], 0.0);
        assert_eq!(tau[1], 0.0);
    }

    #[test]
    fn gravity_scales_with_mass() {
        let mut sys = MultibodySystem::new();
        let si = SpatialInertia::new(
            2.5f64,
            Vector3::new([0.0f64; 3]),
            Matrix3::new([[1.0f64; 3]; 3]),
        );
        let body = RigidBody::new(si, Isometry3::identity(), "heavy", 0);
        sys.add_body(body);
        sys.count_dofs();

        let g = Gravity::new(9.81);
        let tau = g.apply(&sys, &[0.0f64; 6]);
        assert!((tau[2] + 2.5 * 9.81).abs() < 1e-9);
    }
}
