//! Geometric nonlinearity for flexible multibody dynamics.
//!
//! Provides stress stiffening (centrifugal stiffening) and spin softening
//! effects for rotating flexible bodies. These geometric nonlinearities
//! arise from the coupling between large rigid-body rotations and small
//! elastic deformations in the floating frame formulation.
//!
//! # Stress stiffening
//!
//! When a flexible body rotates, the centrifugal forces stiffen the structure:
//!
//! ```text
//! K_geo = ω² · (∫ ρ·(N'ᵀ·N') dV)
//! ```
//!
//! where `ω` is the angular velocity and `N'` are the shape function gradients
//! in the rotating frame.
//!
//! # Spin softening
//!
//! For rotating systems, the stiffness in the direction of rotation is reduced:
//!
//! ```text
//! K_eff = K_elastic - K_geo
//! ```

use tpt_math_linalg_dense::DMatrix;

/// Parameters for geometric nonlinearity effects.
#[derive(Clone, Debug, PartialEq)]
pub struct GeometricNonlinearityParams {
    /// Angular velocity magnitude (rad/s).
    pub angular_velocity: f64,
    /// Rotation axis (unit vector).
    pub rotation_axis: [f64; 3],
    /// Material density (kg/m³).
    pub density: f64,
    /// Stress stiffening coefficient (default: 1.0, full effect).
    pub stress_stiffening_factor: f64,
    /// Spin softening coefficient (default: 1.0, full effect).
    pub spin_softening_factor: f64,
}

impl Default for GeometricNonlinearityParams {
    fn default() -> Self {
        Self {
            angular_velocity: 0.0,
            rotation_axis: [0.0, 0.0, 1.0],
            density: 7850.0,
            stress_stiffening_factor: 1.0,
            spin_softening_factor: 1.0,
        }
    }
}

impl GeometricNonlinearityParams {
    /// Create new geometric nonlinearity parameters.
    pub fn new(angular_velocity: f64, rotation_axis: [f64; 3]) -> Self {
        Self {
            angular_velocity,
            rotation_axis,
            ..Self::default()
        }
    }

    /// Set the material density.
    pub fn with_density(mut self, density: f64) -> Self {
        self.density = density;
        self
    }

    /// Set the stress stiffening factor (0 = no stiffening, 1 = full effect).
    pub fn with_stress_stiffening(mut self, factor: f64) -> Self {
        self.stress_stiffening_factor = factor;
        self
    }

    /// Set the spin softening factor (0 = no softening, 1 = full effect).
    pub fn with_spin_softening(mut self, factor: f64) -> Self {
        self.spin_softening_factor = factor;
        self
    }
}

/// Compute the geometric stiffness matrix due to stress stiffening.
///
/// For a rotating body with angular velocity `ω`, the centrifugal stiffening
/// adds a positive contribution to the stiffness matrix:
///
/// ```text
/// K_stiffening = ω² · M_modal · (Φᵀ·G·Φ)
/// ```
///
/// where `G` is the gyroscopic coupling matrix and `Φ` is the modal basis.
///
/// Returns the geometric stiffness matrix of size `n_modes × n_modes`.
pub fn stress_stiffening_matrix(
    omega: f64,
    modal_mass: &DMatrix<f64>,
    mode_shapes: &DMatrix<f64>,
    rotation_axis: &[f64; 3],
) -> DMatrix<f64> {
    let n_modes = modal_mass.nrows();
    if n_modes == 0 || omega.abs() < 1e-12 {
        return DMatrix::zeros(n_modes, n_modes);
    }

    let omega_sq = omega * omega;
    let mut k_geo = DMatrix::zeros(n_modes, n_modes);

    for i in 0..n_modes {
        for j in 0..n_modes {
            let mut sum = 0.0;
            for k in 0..mode_shapes.nrows() {
                let phi_ki = mode_shapes[(k, i)];
                let phi_kj = mode_shapes[(k, j)];
                let proj = phi_ki * phi_kj;
                let axis_dot = rotation_axis[0] * rotation_axis[0]
                    + rotation_axis[1] * rotation_axis[1]
                    + rotation_axis[2] * rotation_axis[2];
                sum += proj * axis_dot;
            }
            let mi = modal_mass[(i, i)].max(1e-15);
            let mj = modal_mass[(j, j)].max(1e-15);
            k_geo[(i, j)] = omega_sq * sum * (mi + mj) * 0.5;
        }
    }

    k_geo
}

/// Compute the spin softening matrix for rotating bodies.
///
/// When a body spins, the stiffness in the direction of rotation is reduced.
/// This effect is modeled as a negative contribution to the stiffness matrix:
///
/// ```text
/// K_softening = -ω² · diag(M_modal)
/// ```
///
/// Returns the spin softening matrix of size `n_modes × n_modes`.
pub fn spin_softening_matrix(omega: f64, modal_mass: &DMatrix<f64>) -> DMatrix<f64> {
    let n = modal_mass.nrows();
    if n == 0 || omega.abs() < 1e-12 {
        return DMatrix::zeros(n, n);
    }

    let omega_sq = omega * omega;
    let mut k_soft = DMatrix::zeros(n, n);
    for i in 0..n {
        let mi = modal_mass[(i, i)].max(1e-15);
        k_soft[(i, i)] = -omega_sq * mi;
    }
    k_soft
}

/// Compute the effective geometric stiffness matrix.
///
/// Combines stress stiffening and spin softening:
///
/// ```text
/// K_geo = K_stiffening + K_softening
/// ```
///
/// Returns the effective geometric stiffness matrix.
pub fn effective_geometric_stiffness(
    params: &GeometricNonlinearityParams,
    modal_mass: &DMatrix<f64>,
    mode_shapes: &DMatrix<f64>,
) -> DMatrix<f64> {
    let k_stiff = stress_stiffening_matrix(
        params.angular_velocity,
        modal_mass,
        mode_shapes,
        &params.rotation_axis,
    );
    let k_soft = spin_softening_matrix(params.angular_velocity, modal_mass);

    let n = modal_mass.nrows();
    let mut k_geo = DMatrix::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            k_geo[(i, j)] = params.stress_stiffening_factor * k_stiff[(i, j)]
                + params.spin_softening_factor * k_soft[(i, j)];
        }
    }
    k_geo
}

/// Compute the stress stiffening factor for a rotating body.
///
/// Returns a dimensionless factor `α` such that the stiffened stiffness is:
///
/// ```text
/// K_stiffened = K_elastic + α · K_geo
/// ```
///
/// The factor depends on the ratio of centrifugal force to elastic restoring force.
pub fn stress_stiffening_factor(
    angular_velocity: f64,
    radius: f64,
    elastic_modulus: f64,
) -> f64 {
    if elastic_modulus <= 0.0 {
        return 0.0;
    }
    let centrifugal = angular_velocity * angular_velocity * radius;
    (centrifugal / elastic_modulus).abs()
}

/// Compute the spin softening factor for a rotating body.
///
/// Returns a dimensionless factor `β` such that the softened stiffness is:
///
/// ```text
/// K_softened = K_elastic - β · K_geo
/// ```
pub fn spin_softening_factor(angular_velocity: f64, natural_frequency: f64) -> f64 {
    if natural_frequency <= 0.0 {
        return 0.0;
    }
    (angular_velocity / natural_frequency).abs()
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_math_linalg_dense::DMatrix;

    #[test]
    fn test_stress_stiffening_zero_omega() {
        let m = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0]);
        let phi = DMatrix::from_row_slice(4, 2, &[1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]);
        let k = stress_stiffening_matrix(0.0, &m, &phi, &[0.0, 0.0, 1.0]);
        assert_eq!(k.nrows(), 2);
        assert_eq!(k.ncols(), 2);
        for i in 0..2 {
            for j in 0..2 {
                assert!(k[(i, j)].abs() < 1e-12);
            }
        }
    }

    #[test]
    fn test_stress_stiffening_positive() {
        let m = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0]);
        let phi = DMatrix::from_row_slice(4, 2, &[1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]);
        let k = stress_stiffening_matrix(10.0, &m, &phi, &[0.0, 0.0, 1.0]);
        assert!(k[(0, 0)] > 0.0 || k[(1, 1)] > 0.0);
    }

    #[test]
    fn test_spin_softening_negative_diagonal() {
        let m = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0]);
        let k = spin_softening_matrix(10.0, &m);
        assert!(k[(0, 0)] < 0.0);
        assert!(k[(1, 1)] < 0.0);
    }

    #[test]
    fn test_spin_softening_zero_omega() {
        let m = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0]);
        let k = spin_softening_matrix(0.0, &m);
        for i in 0..2 {
            for j in 0..2 {
                assert!(k[(i, j)].abs() < 1e-12);
            }
        }
    }

    #[test]
    fn test_effective_geometric_stiffness_combined() {
        let m = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0]);
        let phi = DMatrix::from_row_slice(4, 2, &[1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]);
        let params = GeometricNonlinearityParams::new(10.0, [0.0, 0.0, 1.0]);
        let k_geo = effective_geometric_stiffness(&params, &m, &phi);
        assert_eq!(k_geo.nrows(), 2);
        assert_eq!(k_geo.ncols(), 2);
    }

    #[test]
    fn test_stress_stiffening_factor() {
        let factor = stress_stiffening_factor(10.0, 0.1, 2.0e11);
        assert!(factor >= 0.0);
        assert!(factor < 1.0);
    }

    #[test]
    fn test_spin_softening_factor() {
        let factor = spin_softening_factor(100.0, 50.0);
        assert!((factor - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_geometric_params_defaults() {
        let params = GeometricNonlinearityParams::default();
        assert_eq!(params.angular_velocity, 0.0);
        assert_eq!(params.rotation_axis, [0.0, 0.0, 1.0]);
    }

    #[test]
    fn test_geometric_params_with_density() {
        let params = GeometricNonlinearityParams::default().with_density(7800.0);
        assert_eq!(params.density, 7800.0);
    }

    #[test]
    fn test_adaptive_stiffness_basic() {
        let k = compute_adaptive_stiffness(1.0, 100.0, 1e-4);
        assert!((k - 1.0e6).abs() < 1.0);
    }
}
