//! Inverse kinematics solvers for serial and parallel manipulators.
//!
//! Implements:
//! - Newton-Raphson with damped least squares (Levenberg-Marquardt)
//! - Jacobian transpose method
//! - Closed-form 6-DOF IK for spherical-wrist manipulators

use core::fmt;

use tpt_math_geometry::Isometry3;
use tpt_math_linalg_fixed::Vector3;

use crate::forward::{geometric_jacobian, DhLink};

/// IK solver configuration.
#[derive(Clone, Debug, PartialEq)]
pub struct IkOptions {
    pub max_iterations: usize,
    pub tolerance_position: f64,
    pub tolerance_orientation: f64,
    pub lambda: f64,
    pub damping: f64,
}

impl Default for IkOptions {
    fn default() -> Self {
        Self {
            max_iterations: 100,
            tolerance_position: 1e-6,
            tolerance_orientation: 1e-6,
            lambda: 0.01,
            damping: 1e-6,
        }
    }
}

/// Result of an IK solve.
#[derive(Clone, Debug, PartialEq)]
pub struct IkResult {
    pub solution: Vec<f64>,
    pub converged: bool,
    pub iterations: usize,
    pub error_position: f64,
    pub error_orientation: f64,
}

impl IkResult {
    pub fn success(
        solution: Vec<f64>,
        iterations: usize,
        error_position: f64,
        error_orientation: f64,
    ) -> Self {
        Self {
            solution,
            converged: true,
            iterations,
            error_position,
            error_orientation,
        }
    }

    pub fn failure(
        solution: Vec<f64>,
        iterations: usize,
        error_position: f64,
        error_orientation: f64,
    ) -> Self {
        Self {
            solution,
            converged: false,
            iterations,
            error_position,
            error_orientation,
        }
    }
}

/// Solve IK using Newton-Raphson with damped least squares (Levenberg-Marquardt).
pub fn solve_newton_lm(
    links: &[DhLink],
    target: &Isometry3<f64>,
    initial_q: &[f64],
    opts: &IkOptions,
) -> IkResult {
    let n = links.len();
    let mut q = initial_q.to_vec();
    let mut q = vec![0.0; n];

    for i in 0..n {
        q[i] = initial_q.get(i).copied().unwrap_or(0.0);
    }

    let identity = Isometry3::<f64>::identity();

    for iter in 0..opts.max_iterations {
        let current_pose = crate::forward::forward_kinematics(links, &q);
        let jac = geometric_jacobian(links, &q);

        let (pos_err, rot_err) = pose_error(&current_pose, target);
        if pos_err < opts.tolerance_position && rot_err < opts.tolerance_orientation {
            return IkResult::success(q, iter + 1, pos_err, rot_err);
        }

        let mut error = [0.0f64; 6];
        error[0] = pos_err;
        error[1] = rot_err;
        error[2] = rot_err;
        error[3] = 0.0;
        error[4] = 0.0;
        error[5] = 0.0;

        let mut jt_j = [[0.0; 6]; 6];
        for col in 0..n {
            let ang = jac.angular_column(col);
            let lin = jac.linear_column(col);
            for row in 0..6 {
                let jval = if row < 3 { ang[row] } else { lin[row - 3] };
                for k in 0..6 {
                    jt_j[k][row] += jval * jval;
                }
            }
        }

        for i in 0..6 {
            jt_j[i][i] += opts.damping;
        }

        let mut delta_q = [0.0f64; 6];
        for i in 0..6.min(n) {
            let mut s = 0.0;
            for k in 0..6 {
                s += jt_j[i][k] * error[k];
            }
            delta_q[i] = s * opts.lambda;
        }

        for i in 0..n {
            q[i] += delta_q[i];
        }
    }

    let final_pose = crate::forward::forward_kinematics(links, &q);
    let (pos_err, rot_err) = pose_error(&final_pose, target);
    IkResult::failure(q, opts.max_iterations, pos_err, rot_err)
}

/// Solve IK using the Jacobian transpose method (gradient descent on pose error).
pub fn solve_jacobian_transpose(
    links: &[DhLink],
    target: &Isometry3<f64>,
    initial_q: &[f64],
    opts: &IkOptions,
) -> IkResult {
    let n = links.len();
    let mut q = vec![0.0; n];
    for i in 0..n {
        q[i] = initial_q.get(i).copied().unwrap_or(0.0);
    }

    let alpha = 0.5;

    for iter in 0..opts.max_iterations {
        let current_pose = crate::forward::forward_kinematics(links, &q);
        let (pos_err, rot_err) = pose_error(&current_pose, target);

        if pos_err < opts.tolerance_position && rot_err < opts.tolerance_orientation {
            return IkResult::success(q, iter + 1, pos_err, rot_err);
        }

        let jac = geometric_jacobian(links, &q);
        for i in 0..n {
            let ang = jac.angular_column(i);
            let lin = jac.linear_column(i);
            let ang_norm = ang[0] * ang[0] + ang[1] * ang[1] + ang[2] * ang[2];
            let lin_norm = lin[0] * lin[0] + lin[1] * lin[1] + lin[2] * lin[2];
            let j_norm_sq = ang_norm + lin_norm;
            if j_norm_sq > 1e-12 {
                let scale = alpha * (pos_err + rot_err) / j_norm_sq.sqrt();
                q[i] += scale * (ang[0] + lin[0]);
            }
        }
    }

    let final_pose = crate::forward::forward_kinematics(links, &q);
    let (pos_err, rot_err) = pose_error(&final_pose, target);
    IkResult::failure(q, opts.max_iterations, pos_err, rot_err)
}

/// Closed-form IK for 6-DOF manipulators with a spherical wrist.
///
/// Uses PUMA-style kinematic decoupling: solve wrist center position first,
/// then spherical wrist orientation.
pub fn solve_spherical_wrist(
    links: &[DhLink],
    target: &Isometry3<f64>,
    elbow_config: ElbowConfig,
) -> Option<IkResult> {
    if links.len() < 6 {
        return None;
    }

    let n = links.len();
    let mut q = vec![0.0; n];

    let t_target = *target;
    let r = t_target.rotation.matrix();
    let t = t_target.translation.vector;

    let d6 = links[5].d;
    let a3 = links[2].a;

    let wrist_x = t.data[0] - d6 * r.data[0][2];
    let wrist_y = t.data[1] - d6 * r.data[1][2];
    let wrist_z = t.data[2] - d6 * r.data[2][2];

    let q1 = wrist_y.atan2(wrist_x);

    let r_sq = wrist_x * wrist_x + wrist_y * wrist_y;
    let q2_arg = (wrist_z - links[0].d) / r_sq.sqrt();
    let q2 = q2_arg.atan2(1.0);

    let q3_arg = (r_sq + (wrist_z - links[0].d).powi(2) - a3.powi(2) - links[1].a.powi(2))
        / (2.0 * a3 * links[1].a);
    let q3 = q3_arg.clamp(-1.0, 1.0).acos();

    match elbow_config {
        ElbowConfig::Down => {
            q[2] = q3;
        }
        ElbowConfig::Up => {
            q[2] = -q3;
        }
    }

    q[0] = q1;
    q[1] = q2;

    let pose = crate::forward::forward_kinematics(links, &q);
    let (pos_err, rot_err) = pose_error(&pose, target);

    if pos_err < 1e-6 && rot_err < 1e-6 {
        Some(IkResult::success(q, 1, pos_err, rot_err))
    } else {
        Some(IkResult::failure(q, 1, pos_err, rot_err))
    }
}

/// Elbow configuration for spherical-wrist IK.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ElbowConfig {
    Up,
    Down,
}

/// Compute pose error between current and target.
fn pose_error(current: &Isometry3<f64>, target: &Isometry3<f64>) -> (f64, f64) {
    let dp = current.translation.vector - target.translation.vector;
    let pos_err =
        (dp.data[0] * dp.data[0] + dp.data[1] * dp.data[1] + dp.data[2] * dp.data[2]).sqrt();

    let rot_err = current.rotation.matrix().data[0][0]
        + current.rotation.matrix().data[1][1]
        + current.rotation.matrix().data[2][2];

    let trace = (rot_err - 1.0) / 2.0;
    let trace = trace.clamp(-1.0, 1.0);
    let rot_err = trace.acos();

    (pos_err, rot_err)
}

impl fmt::Display for IkResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "IkResult(converged={}, iter={}, pos_err={:.6}, rot_err={:.6})",
            self.converged, self.iterations, self.error_position, self.error_orientation
        )
    }
}
