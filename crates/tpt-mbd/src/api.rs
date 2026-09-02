use crate::error::{MbdError, Result};
#[cfg(feature = "kinematics")]
use tpt_math_geometry::Isometry3;
#[cfg(feature = "kinematics")]
use tpt_mbd_kinematics::chain::DhChain;
#[cfg(feature = "kinematics")]
use tpt_mbd_kinematics::inverse::{solve_newton_lm, IkOptions, IkResult};
#[cfg(feature = "system")]
use tpt_mbd_system::system::MultibodySystem;

/// Result of a full multibody simulation.
#[derive(Debug, Clone, Default)]
pub struct SimulationResult {
    /// Time stamps for each recorded state.
    pub times: Vec<f64>,
    /// Generalized coordinates at each time step.
    pub states: Vec<Vec<f64>>,
}

/// Compute forward kinematics for a serial chain.
///
/// When the `kinematics` feature is disabled this returns a stub.
///
/// # Examples
///
/// ```
/// use tpt_mbd::api::forward_kinematics;
/// use tpt_mbd_kinematics::chain::DhLink;
///
/// let links = vec![DhLink::new(0.0, 0.0, 0.0, 0.0)];
/// let pose = forward_kinematics(&links, &[0.0]);
/// assert!(pose.rotation.matrix().data[0][0].is_finite());
/// ```
#[cfg(feature = "kinematics")]
pub fn forward_kinematics(chain: &DhChain, joint_angles: &[f64]) -> Isometry3<f64> {
    chain.forward(joint_angles)
}

/// Compute forward kinematics for a serial chain (stub, no `kinematics` feature).
#[cfg(not(feature = "kinematics"))]
pub fn forward_kinematics(_chain: &(), _joint_angles: &[f64]) {}

/// Solve inverse kinematics for a serial chain.
///
/// Returns `Ok(IkResult)` on convergence, or `Err(MbdError)` with a
/// `KinematicsErrorKind::IkNotConverged` / `InvalidChain` detail on failure.
///
/// Requires the `kinematics` feature.
#[cfg(feature = "kinematics")]
pub fn inverse_kinematics(
    chain: &DhChain,
    target: &Isometry3<f64>,
    initial: &[f64],
) -> Result<IkResult> {
    let links = &chain.links;
    if links.is_empty() {
        return Err(MbdError::KinematicsError {
            message: "chain has no links".to_string(),
            kind: crate::error::KinematicsErrorKind::InvalidChain,
            source: None,
        });
    }
    let opts = IkOptions::default();
    let result = solve_newton_lm(links, target, initial, &opts);
    if result.converged {
        Ok(result)
    } else {
        Err(MbdError::KinematicsError {
            message: format!(
                "IK did not converge after {} iterations (pos_err={:.3e}, rot_err={:.3e})",
                result.iterations, result.error_position, result.error_orientation
            ),
            kind: crate::error::KinematicsErrorKind::IkNotConverged,
            source: None,
        })
    }
}

/// Solve inverse kinematics (stub, no `kinematics` feature).
#[cfg(not(feature = "kinematics"))]
pub fn inverse_kinematics(_chain: &(), _target: &(), _initial: &[f64]) -> Result<()> {
    Err(MbdError::KinematicsError {
        message: "kinematics feature not enabled".to_string(),
        kind: crate::error::KinematicsErrorKind::InvalidChain,
        source: None,
    })
}

/// Compute inverse dynamics (constraint forces for given accelerations).
///
/// Requires the `system` feature.  Returns the generalized forces `τ = M·q̈`
/// computed from the system mass matrix and the provided accelerations.
#[cfg(feature = "system")]
pub fn inverse_dynamics(
    system: &MultibodySystem,
    _q: &[f64],
    _qdot: &[f64],
    qddot: &[f64],
    _tau: &[f64],
) -> Result<Vec<f64>> {
    Ok(system.inverse_dynamics(qddot))
}

/// Compute inverse dynamics (stub, no `system` feature).
#[cfg(not(feature = "system"))]
pub fn inverse_dynamics(
    _system: &(),
    _q: &[f64],
    _qdot: &[f64],
    _qddot: &[f64],
    _tau: &[f64],
) -> Result<Vec<f64>> {
    Err(MbdError::SystemError {
        message: "system feature not enabled".to_string(),
        kind: crate::error::SystemErrorKind::Unconstrained,
    })
}

/// Run a full multibody simulation from `t = 0` to `t_final`.
///
/// Requires the `system` feature.  Returns an empty `SimulationResult` until
/// `tpt-mbd-system` exposes an integrator.
///
/// # Examples
///
/// ```
/// use tpt_mbd::api::simulate;
/// use tpt_mbd_system::builder::MultibodySystemBuilder;
///
/// let system = MultibodySystemBuilder::new("test")
///     .add_body("ground")
///     .build()
///     .expect("builder should succeed");
/// let result = simulate(&system, 1.0, 0.01);
/// assert!(result.times.is_empty());
/// ```
#[cfg(feature = "system")]
pub fn simulate(system: &MultibodySystem, _t_final: f64, _dt: f64) -> SimulationResult {
    let _ = system;
    SimulationResult::default()
}

/// Run a full multibody simulation (stub, no `system` feature).
#[cfg(not(feature = "system"))]
pub fn simulate(_system: &(), _t_final: f64, _dt: f64) -> SimulationResult {
    SimulationResult::default()
}
