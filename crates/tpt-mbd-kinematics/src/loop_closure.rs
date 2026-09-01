//! Loop-closure and parallel mechanism kinematics.
//!
//! Provides:
//! - [`ParallelChain`] — a collection of serial chains sharing a common
//!   end-effector pose
//! - [`loop_closure_error`] — pose discrepancy for a single chain against the
//!   target
//! - [`parallel_forward_kinematics`] — compute end-effector pose from all
//!   chains
//! - [`solve_loop_closure`] — simple coordinate descent to minimize loop-
//!   closure error across all chains

extern crate alloc;

use alloc::vec::Vec;

use tpt_math_geometry::Isometry3;
use tpt_math_linalg_fixed::Vector3;

use crate::chain::DhLink;
use crate::forward::forward_kinematics;

/// A single serial chain in a parallel mechanism.
#[derive(Clone, Debug)]
pub struct Chain {
    /// DH parameters for this chain.
    pub links: Vec<DhLink>,
    /// Current joint angles.
    pub joint_angles: Vec<f64>,
}

impl Chain {
    /// Create a new chain from DH links and initial joint angles.
    pub fn new(links: Vec<DhLink>, joint_angles: Vec<f64>) -> Self {
        Self {
            links,
            joint_angles,
        }
    }

    /// Evaluate forward kinematics for this chain.
    pub fn forward(&self) -> Isometry3<f64> {
        forward_kinematics(&self.links, &self.joint_angles)
    }
}

/// A parallel mechanism consisting of multiple serial chains sharing a common
/// end-effector.
///
/// All chains are rooted at the same base frame and must reach the same
/// end-effector pose.  Loop-closure constraints enforce this consistency.
#[derive(Clone, Debug)]
pub struct ParallelChain {
    /// Serial chains that form the parallel mechanism.
    pub chains: Vec<Chain>,
    /// Target end-effector pose (shared by all chains).
    pub target: Isometry3<f64>,
    /// Tolerance for loop-closure error (meters + radians).
    pub tolerance: f64,
}

impl ParallelChain {
    /// Create a new parallel mechanism.
    ///
    /// * `chains` — serial chains that connect the base to the common
    ///   end-effector
    /// * `target` — desired end-effector pose
    pub fn new(chains: Vec<Chain>, target: Isometry3<f64>) -> Self {
        Self {
            chains,
            target,
            tolerance: 1e-6,
        }
    }

    /// Compute loop-closure error for a specific chain.
    ///
    /// Returns a 6-vector `[translation_error; rotation_error]` where
    /// translation error is in meters and rotation error is in radians.
    pub fn loop_closure_error(&self, chain_idx: usize) -> [f64; 6] {
        let chain = &self.chains[chain_idx];
        let actual = chain.forward();
        let dt = actual.translation.vector - self.target.translation.vector;
        let tx = dt.x();
        let ty = dt.y();
        let tz = dt.z();

        let rot_err = rotation_error(&actual.rotation, &self.target.rotation);

        [tx, ty, tz, rot_err[0], rot_err[1], rot_err[2]]
    }

    /// Compute the maximum loop-closure error across all chains.
    ///
    /// Returns `(max_translation_error, max_rotation_error)`.
    pub fn max_error(&self) -> (f64, f64) {
        let mut max_t = 0.0f64;
        let mut max_r = 0.0f64;
        for i in 0..self.chains.len() {
            let err = self.loop_closure_error(i);
            let t = (err[0] * err[0] + err[1] * err[1] + err[2] * err[2]).sqrt();
            let r = (err[3] * err[3] + err[4] * err[4] + err[5] * err[5]).sqrt();
            max_t = max_t.max(t);
            max_r = max_r.max(r);
        }
        (max_t, max_r)
    }

    /// Check whether all chains satisfy the loop-closure tolerance.
    pub fn is_converged(&self) -> bool {
        let (t, r) = self.max_error();
        t < self.tolerance && r < self.tolerance
    }

    /// Solve loop-closure via simple coordinate-descent.
    ///
    /// Iteratively adjusts joint angles for each chain to minimize the
    /// loop-closure error.  Returns the number of iterations performed.
    pub fn solve_loop_closure(&mut self, max_iters: usize) -> usize {
        for iter in 0..max_iters {
            if self.is_converged() {
                return iter;
            }
            let errors: alloc::vec::Vec<[f64; 6]> = (0..self.chains.len())
                .map(|i| self.loop_closure_error(i))
                .collect();
            for (chain, err) in self.chains.iter_mut().zip(errors.iter()) {
                let scale = 0.1;
                for (a, e) in chain.joint_angles.iter_mut().zip(err.iter().cycle()) {
                    *a -= scale * e;
                }
            }
        }
        max_iters
    }
}

/// Compute the axis-angle rotation error between two orientations.
fn rotation_error(
    actual: &tpt_math_geometry::Rotation3<f64>,
    target: &tpt_math_geometry::Rotation3<f64>,
) -> [f64; 3] {
    let rel = *target * actual.inverse();
    let trace = rel.matrix().data[0][0] + rel.matrix().data[1][1] + rel.matrix().data[2][2];
    let clamped = (trace - 1.0).clamp(-1.0, 1.0);
    let angle = clamped.acos();
    let axis = Vector3::new([
        rel.matrix().data[2][1] - rel.matrix().data[1][2],
        rel.matrix().data[0][2] - rel.matrix().data[2][0],
        rel.matrix().data[1][0] - rel.matrix().data[0][1],
    ]);
    let norm = (axis.x() * axis.x() + axis.y() * axis.y() + axis.z() * axis.z()).sqrt();
    if norm < 1e-12 {
        [0.0, 0.0, 0.0]
    } else {
        let s = angle / norm;
        [axis.x() * s, axis.y() * s, axis.z() * s]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_math_geometry::Isometry3;

    #[test]
    fn parallel_chain_single_chain_zero_error() {
        let links = vec![DhLink::new(0.0, 0.0, 0.0, 0.0)];
        let chain = Chain::new(links, vec![0.0]);
        let target = Isometry3::<f64>::identity();
        let parallel = ParallelChain::new(vec![chain], target);
        let err = parallel.loop_closure_error(0);
        for &e in &err {
            assert!(e.abs() < 1e-10, "error too large: {}", e);
        }
    }

    #[test]
    fn parallel_chain_converges_to_target() {
        let links = vec![DhLink::new(0.0, 0.0, 0.0, 0.0)];
        let chain = Chain::new(links, vec![0.1]);
        let target = Isometry3::<f64>::identity();
        let mut parallel = ParallelChain::new(vec![chain], target);
        let iters = parallel.solve_loop_closure(100);
        assert!(
            parallel.is_converged(),
            "did not converge after {} iters",
            iters
        );
    }

    #[test]
    fn parallel_chain_two_chains_consistency() {
        let links = vec![DhLink::new(0.0, 0.0, 0.0, 0.0)];
        let chain0 = Chain::new(links.clone(), vec![0.05]);
        let chain1 = Chain::new(links, vec![-0.05]);
        let target = Isometry3::<f64>::identity();
        let mut parallel = ParallelChain::new(vec![chain0, chain1], target);
        let iters = parallel.solve_loop_closure(100);
        assert!(
            parallel.is_converged(),
            "did not converge after {} iters",
            iters
        );
    }
}
