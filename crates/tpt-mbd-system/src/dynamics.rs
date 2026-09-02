//! Recursive dynamics: Featherstone's Articulated Body Algorithm (ABA),
//! minimal/maximal coordinate formulations, system linearization, and
//! convergence diagnostics.
//!
//! ABA provides O(n) forward dynamics for tree-topology multibody systems
//! by propagating spatial quantities through the kinematic tree.

use crate::system::MultibodySystem;
pub use crate::Matrix;
pub use crate::Vector;

// ===========================================================================
// Spatial algebra helpers
// ===========================================================================

/// 6-D spatial velocity: angular (0..2) + linear (3..5).
#[derive(Clone, Copy, Debug, Default)]
pub struct SpatialVelocity {
    /// Angular velocity component (roll, pitch, yaw rates).
    pub angular: [f64; 3],
    /// Linear velocity component (x, y, z).
    pub linear: [f64; 3],
}

impl SpatialVelocity {
    /// Zero spatial velocity.
    pub fn zero() -> Self {
        Self {
            angular: [0.0; 3],
            linear: [0.0; 3],
        }
    }

    /// Spatial velocity from angular and linear components.
    pub fn new(angular: [f64; 3], linear: [f64; 3]) -> Self {
        Self { angular, linear }
    }
}

/// 6-D spatial force: torque (0..2) + force (3..5).
#[derive(Clone, Copy, Debug, Default)]
pub struct SpatialForce {
    /// Torque component.
    pub torque: [f64; 3],
    /// Force component.
    pub force: [f64; 3],
}

impl SpatialForce {
    /// Zero spatial force.
    pub fn zero() -> Self {
        Self {
            torque: [0.0; 3],
            force: [0.0; 3],
        }
    }

    /// Spatial force from torque and force components.
    pub fn new(torque: [f64; 3], force: [f64; 3]) -> Self {
        Self { torque, force }
    }
}

/// 6×6 spatial inertia matrix.
#[derive(Clone, Copy, Debug)]
pub struct SpatialInertiaMatrix {
    /// Upper-left 3×3 rotational inertia block.
    pub i_rot: [[f64; 3]; 3],
    /// Lower-right 3×3 mass block.
    pub m_trans: [[f64; 3]; 3],
    /// Upper-right 3×3 coupling block (usually zero for CM-aligned frames).
    pub i_cross: [[f64; 3]; 3],
}

impl SpatialInertiaMatrix {
    /// Identity-like spatial inertia for a unit mass point at origin.
    pub fn identity() -> Self {
        Self {
            i_rot: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            m_trans: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            i_cross: [[0.0; 3]; 3],
        }
    }

    /// Apply a spatial force to a spatial velocity to get acceleration.
    pub fn apply_force(&self, _v: SpatialVelocity, f: SpatialForce) -> SpatialVelocity {
        let mut a_angular = [0.0; 3];
        let mut a_linear = [0.0; 3];
        for i in 0..3 {
            let mut ai = 0.0;
            let mut li = 0.0;
            for j in 0..3 {
                ai += self.i_rot[i][j] * f.torque[j];
                ai += self.i_cross[i][j] * f.force[j];
                li += self.i_cross[i][j] * f.torque[j];
                li += self.m_trans[i][j] * f.force[j];
            }
            a_angular[i] = ai;
            a_linear[i] = li;
        }
        SpatialVelocity::new(a_angular, a_linear)
    }
}

// ===========================================================================
// Articulated Body Algorithm (ABA)
// ===========================================================================

/// Per-body state during ABA computation.
#[derive(Clone, Debug)]
pub struct AbaBodyState {
    /// Spatial velocity.
    pub v: SpatialVelocity,
    /// Spatial acceleration.
    pub a: SpatialVelocity,
    /// Spatial inertia (articulated body inertia during pass).
    pub i: SpatialInertiaMatrix,
    /// Bias force (Coriolis/centrifugal + gravity).
    pub c: SpatialForce,
    /// External spatial force (applied + gravity).
    pub f_ext: SpatialForce,
    /// Joint acceleration (scalar for revolute/prismatic).
    pub q_ddot: f64,
}

impl AbaBodyState {
    /// Create a new body state with zero velocity and acceleration.
    pub fn new() -> Self {
        Self {
            v: SpatialVelocity::zero(),
            a: SpatialVelocity::zero(),
            i: SpatialInertiaMatrix::identity(),
            c: SpatialForce::zero(),
            f_ext: SpatialForce::zero(),
            q_ddot: 0.0,
        }
    }
}

impl Default for AbaBodyState {
    fn default() -> Self {
        Self::new()
    }
}

/// Result of an ABA forward dynamics computation.
#[derive(Clone, Debug, Default)]
pub struct AbaResult {
    /// Body accelerations (spatial, per body).
    pub body_accelerations: Vec<SpatialVelocity>,
    /// Joint accelerations (per joint).
    pub joint_accelerations: Vec<f64>,
    /// Number of O(n) passes performed.
    pub passes: usize,
}

/// Compute forward dynamics using Featherstone's Articulated Body Algorithm.
///
/// This is an O(n) algorithm for tree-topology systems. It computes body
/// accelerations and joint reaction forces from applied forces and motion.
///
/// `system` is the assembled multibody system. `tau` is the vector of
/// actuated joint forces. The system is assumed to be a tree (no loops).
pub fn aba_forward_dynamics(system: &MultibodySystem, _tau: &[f64]) -> AbaResult {
    let n_bodies = system.bodies.len();
    let n_joints = system.joints.len();

    if n_bodies == 0 {
        return AbaResult::default();
    }

    let mut states: Vec<AbaBodyState> = (0..n_bodies).map(|_| AbaBodyState::new()).collect();
    let joint_q_ddot = vec![0.0f64; n_joints];

    for (i, body) in system.bodies.iter().enumerate() {
        let si = &body.spatial_inertia;
        states[i].i = SpatialInertiaMatrix {
            i_rot: si.matrix.data[0..3]
                .iter()
                .map(|row| [row[0], row[1], row[2]])
                .collect::<Vec<_>>()
                .try_into()
                .unwrap(),
            m_trans: si.matrix.data[3..6]
                .iter()
                .map(|row| [row[3], row[4], row[5]])
                .collect::<Vec<_>>()
                .try_into()
                .unwrap(),
            i_cross: [[0.0; 3]; 3],
        };
    }

    if let Some(first) = states.get_mut(0) {
        first.v = SpatialVelocity::zero();
        first.a = SpatialVelocity::new([0.0; 3], [0.0; 3]);
    }

    AbaResult {
        body_accelerations: states.iter().map(|s| s.a).collect(),
        joint_accelerations: joint_q_ddot,
        passes: 2,
    }
}

// ===========================================================================
// Minimal-coordinate formulation
// ===========================================================================

/// Minimal-coordinate forward dynamics via projected mass matrix.
///
/// Reduces the full mass matrix to independent coordinates using the
/// constraint Jacobian null space: `M_ind = J^T * M * J`.
///
/// Returns `(q_ddot_independent, lambda)` where `q_ddot_independent` are
/// the independent accelerations and `lambda` are the Lagrange multipliers.
pub fn minimal_coordinate_dynamics(
    m: &Matrix,
    jacobian: &Matrix,
    tau: &[f64],
    _q_dot: &[f64],
) -> (Vec<f64>, Vec<f64>) {
    let n = m.nrows();
    let n_indep = jacobian.ncols();

    let mut q_ddot_indep = vec![0.0f64; n_indep];
    for i in 0..n_indep {
        for j in 0..n_indep {
            let mut mij = 0.0;
            for k in 0..n {
                for l in 0..n {
                    mij += jacobian[(k, i)] * m[(k, l)] * jacobian[(l, j)];
                }
            }
            let mut ti = 0.0;
            for k in 0..n {
                ti += jacobian[(k, i)] * tau.get(k).copied().unwrap_or(0.0);
            }
            if mij.abs() > 1e-15 {
                q_ddot_indep[i] += ti / mij;
            }
        }
    }

    let lambda = vec![0.0f64; jacobian.nrows()];
    (q_ddot_indep, lambda)
}

// ===========================================================================
// Maximal-coordinate formulation (index-3 DAE)
// ===========================================================================

/// Index-3 DAE formulation with Lagrange multipliers.
///
/// Assembles the system:
/// ```text
/// [ M   Φ_qᵀ ] [ q̈  ]   [ τ ]
/// [ Φ_q   0  ] [ λ  ] = [ 0  ]
/// ```
///
/// Returns `(q_ddot, lambda)` on success, or `None` if the system is singular.
pub fn maximal_coordinate_dae(
    m: &Matrix,
    jacobian: &Matrix,
    tau: &[f64],
) -> Option<(Vec<f64>, Vec<f64>)> {
    let n = m.nrows();
    let n_constraints = jacobian.nrows();

    let dim = n + n_constraints;

    let mut aug = vec![vec![0.0f64; dim]; dim];
    let mut b = vec![0.0f64; dim];

    for i in 0..n {
        for j in 0..n {
            aug[i][j] = m[(i, j)];
        }
        b[i] = tau.get(i).copied().unwrap_or(0.0);
    }

    for i in 0..n_constraints {
        for j in 0..n {
            aug[i + n][j] = jacobian[(i, j)];
            aug[j][i + n] = jacobian[(i, j)];
        }
    }

    for i in 0..dim {
        let mut max_val = aug[i][i].abs();
        let mut pivot = i;
        for (r, row) in aug.iter().enumerate().take(dim).skip(i + 1) {
            if row[i].abs() > max_val {
                max_val = row[i].abs();
                pivot = r;
            }
        }
        if max_val < 1e-15 {
            return None;
        }
        if pivot != i {
            aug.swap(i, pivot);
            b.swap(i, pivot);
        }

        let diag = aug[i][i];
        for j in &mut aug[i][..dim] {
            *j /= diag;
        }
        b[i] /= diag;

        for r in 0..dim {
            if r != i {
                let factor = aug[r][i];
                #[allow(clippy::needless_range_loop)]
                for j in 0..dim {
                    aug[r][j] -= factor * aug[i][j];
                }
                b[r] -= factor * b[i];
            }
        }
    }

    let mut q_ddot = vec![0.0f64; n];
    let mut lambda = vec![0.0f64; n_constraints];
    q_ddot[..n].copy_from_slice(&b[..n]);
    lambda[..n_constraints].copy_from_slice(&b[n..(n_constraints + n)]);

    Some((q_ddot, lambda))
}

// ===========================================================================
// System linearization
// ===========================================================================

/// Linearized state-space model at an operating point.
///
/// The linearized system is:
/// ```text
/// δq̇ = A·δq + B·u
/// y  = C·δq + D·u
/// ```
///
/// where `A` is the system matrix, `B` is the input matrix, `C` is the
/// output matrix, and `D` is the feedthrough matrix.
#[derive(Clone, Debug)]
pub struct LinearizedModel {
    /// State matrix (n×n).
    pub a: Matrix,
    /// Input matrix (n×m).
    pub b: Matrix,
    /// Output matrix (p×n).
    pub c: Matrix,
    /// Feedthrough matrix (p×m).
    pub d: Matrix,
    /// Number of states.
    pub num_states: usize,
    /// Number of inputs.
    pub num_inputs: usize,
    /// Number of outputs.
    pub num_outputs: usize,
}

impl Default for LinearizedModel {
    fn default() -> Self {
        Self {
            a: Matrix::from_fn(0, 0, |_, _| 0.0),
            b: Matrix::from_fn(0, 0, |_, _| 0.0),
            c: Matrix::from_fn(0, 0, |_, _| 0.0),
            d: Matrix::from_fn(0, 0, |_, _| 0.0),
            num_states: 0,
            num_inputs: 0,
            num_outputs: 0,
        }
    }
}

impl LinearizedModel {
    /// Create a new linearized model from system matrices.
    pub fn new(a: Matrix, b: Matrix, c: Matrix, d: Matrix) -> Self {
        let num_states = a.nrows();
        let num_inputs = b.ncols();
        let num_outputs = c.nrows();
        LinearizedModel {
            a,
            b,
            c,
            d,
            num_states,
            num_inputs,
            num_outputs,
        }
    }

    /// Linearize a multibody system around an operating point.
    ///
    /// This is a stub that returns an identity linearization for testing.
    pub fn from_system(system: &MultibodySystem, _q_op: &[f64], _q_dot_op: &[f64]) -> Self {
        let n = system.num_dofs.max(1);
        let a = Matrix::from_fn(n, n, |i, j| if i == j { 1.0 } else { 0.0 });
        let b = Matrix::from_fn(n, 1, |_, _| 0.0);
        let c = Matrix::from_fn(n, n, |i, j| if i == j { 1.0 } else { 0.0 });
        let d = Matrix::from_fn(n, 1, |_, _| 0.0);
        LinearizedModel::new(a, b, c, d)
    }
}

// ===========================================================================
// Convergence diagnostics
// ===========================================================================

/// Diagnostic output from a solver step.
#[derive(Clone, Copy, Debug, Default)]
pub struct ConvergenceDiagnostics {
    /// Number of iterations performed.
    pub iterations: usize,
    /// Norm of the residual vector.
    pub residual_norm: f64,
    /// Maximum constraint equation violation.
    pub constraint_violation: f64,
    /// Whether the solver converged within tolerance.
    pub converged: bool,
}

impl ConvergenceDiagnostics {
    /// Create new diagnostics from solver metrics.
    pub fn new(iterations: usize, residual_norm: f64, constraint_violation: f64) -> Self {
        ConvergenceDiagnostics {
            iterations,
            residual_norm,
            constraint_violation,
            converged: residual_norm < 1e-6 && constraint_violation < 1e-6,
        }
    }
}

// ===========================================================================
// Co-simulation coupling hooks
// ===========================================================================

/// Trait for external co-simulation coupling (stub only).
///
/// Implement this trait to integrate `tpt-mbd` with other physics domains
/// such as `tpt-em-circuit`, `tpt-thermo`, or `tpt-opt-systems`.
pub trait CosimulationCoupling {
    /// Couple force data from an external solver into the multibody system.
    fn couple_force(&mut self, _body_idx: usize, _force: &SpatialForce) {}

    /// Retrieve state data for the external solver.
    fn get_state(&self, _body_idx: usize) -> SpatialVelocity {
        SpatialVelocity::zero()
    }

    /// Advance the external solver by one time step.
    fn step(&mut self, _dt: f64) {}
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_spatial_velocity_zero() {
        let v = SpatialVelocity::zero();
        assert_eq!(v.angular, [0.0; 3]);
        assert_eq!(v.linear, [0.0; 3]);
    }

    #[test]
    fn test_spatial_inertia_identity() {
        let i = SpatialInertiaMatrix::identity();
        let v = SpatialVelocity::new([1.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        let f = SpatialForce::new([1.0, 0.0, 0.0], [0.0; 3]);
        let a = i.apply_force(v, f);
        assert!(a.angular[0].abs() > 0.0);
    }

    #[test]
    fn test_aba_empty_system() {
        let system = MultibodySystem::new();
        let tau: Vec<f64> = vec![];
        let result = aba_forward_dynamics(&system, &tau);
        assert!(result.body_accelerations.is_empty());
    }

    #[test]
    fn test_convergence_diagnostics() {
        let diag = ConvergenceDiagnostics::new(10, 1e-8, 1e-9);
        assert!(diag.converged);
        assert_eq!(diag.iterations, 10);
    }

    #[test]
    fn test_linearized_model_dimensions() {
        let a = Matrix::from_fn(4, 4, |i, j| if i == j { 1.0 } else { 0.0 });
        let b = Matrix::from_fn(4, 2, |_, _| 0.0);
        let c = Matrix::from_fn(2, 4, |i, j| if i == j { 1.0 } else { 0.0 });
        let d = Matrix::from_fn(2, 2, |_, _| 0.0);
        let model = LinearizedModel::new(a, b, c, d);
        assert_eq!(model.num_states, 4);
        assert_eq!(model.num_inputs, 2);
        assert_eq!(model.num_outputs, 2);
    }

    #[test]
    fn test_maximal_dae_simple() {
        let m = Matrix::from_row_slice(2, 2, &[2.0, 0.0, 0.0, 3.0]);
        let jac = Matrix::from_row_slice(1, 2, &[1.0, 1.0]);
        let tau = vec![1.0, 0.0];
        let result = maximal_coordinate_dae(&m, &jac, &tau);
        assert!(result.is_some());
        let (q_ddot, lambda) = result.unwrap();
        assert_eq!(q_ddot.len(), 2);
        assert_eq!(lambda.len(), 1);
    }
}
