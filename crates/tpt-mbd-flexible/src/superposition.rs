//! Modal superposition for flexible multibody dynamics.
//!
//! Provides the time-integration framework for flexible bodies represented in
//! modal coordinates. The displacement of a flexible body is expressed as:
//!
//! ```text
//! u = Φ · q_modal
//! ```
//!
//! where `Φ` is the modal basis (from Craig-Bampton CMS) and `q_modal` are
//! the time-varying generalized coordinates. The reduced equations of motion
//! are:
//!
//! ```text
//! M_red · q̈ + C_red · q̇ + K_red · q = f_modal
//! ```
//!
//! This module handles the assembly of reduced matrices, modal force
//! projection, and time integration of the modal coordinates.

use tpt_math_linalg_dense::DMatrix;

use crate::cms::{select_modes, CraigBampton, ModalBasis, ModeSelection};
use crate::damping::{compute_damping_matrix, RayleighDamping};

// ===========================================================================
// Modal superposition state
// ===========================================================================

/// State of a flexible body in modal coordinates.
#[derive(Clone, Debug)]
pub struct ModalSuperpositionState {
    /// Current modal coordinates (displacement amplitudes).
    pub q: Vec<f64>,
    /// Current modal velocities.
    pub q_dot: Vec<f64>,
    /// Current modal accelerations.
    pub q_ddot: Vec<f64>,
    /// Reduced mass matrix (from CMS projection).
    pub m_red: DMatrix<f64>,
    /// Reduced stiffness matrix (from CMS projection).
    pub k_red: DMatrix<f64>,
    /// Reduced damping matrix (from Rayleigh damping).
    pub c_red: DMatrix<f64>,
    /// Number of retained modes.
    pub num_modes: usize,
}

impl Default for ModalSuperpositionState {
    fn default() -> Self {
        Self {
            q: Vec::new(),
            q_dot: Vec::new(),
            q_ddot: Vec::new(),
            m_red: DMatrix::zeros(0, 0),
            k_red: DMatrix::zeros(0, 0),
            c_red: DMatrix::zeros(0, 0),
            num_modes: 0,
        }
    }
}

impl ModalSuperpositionState {
    /// Create a new modal superposition state from a Craig-Bampton reduced model.
    pub fn from_craig_bampton(
        cb: &CraigBampton,
        m_full: DMatrix<f64>,
        k_full: DMatrix<f64>,
        rayleigh: &RayleighDamping,
    ) -> Self {
        let (m_red, k_red) = crate::cms::reduce_system(cb, m_full, k_full);
        let c_red = compute_damping_matrix(rayleigh, m_red.clone(), k_red.clone());
        let num_modes = m_red.nrows();
        ModalSuperpositionState {
            q: vec![0.0; num_modes],
            q_dot: vec![0.0; num_modes],
            q_ddot: vec![0.0; num_modes],
            m_red,
            k_red,
            c_red,
            num_modes,
        }
    }

    /// Create a modal superposition state from pre-computed reduced matrices.
    pub fn from_reduced_matrices(
        m_red: DMatrix<f64>,
        k_red: DMatrix<f64>,
        c_red: DMatrix<f64>,
    ) -> Self {
        let num_modes = m_red.nrows();
        ModalSuperpositionState {
            q: vec![0.0; num_modes],
            q_dot: vec![0.0; num_modes],
            q_ddot: vec![0.0; num_modes],
            m_red,
            k_red,
            c_red,
            num_modes,
        }
    }

    /// Compute the full displacement vector from modal coordinates.
    ///
    /// `phi` is the modal basis matrix (`n_dofs × n_modes`). The returned
    /// vector has length `n_dofs`.
    pub fn displacement(&self, phi: &DMatrix<f64>) -> Vec<f64> {
        let n = phi.nrows();
        let mut u = vec![0.0; n];
        for i in 0..n {
            for j in 0..self.num_modes {
                u[i] += phi[(i, j)] * self.q[j];
            }
        }
        u
    }

    /// Advance the modal state by one time step using semi-implicit Euler.
    ///
    /// `f_modal` is the generalized force vector in modal coordinates.
    pub fn step_semi_implicit_euler(&mut self, f_modal: &[f64], dt: f64) {
        let n = self.num_modes;
        if n == 0 {
            return;
        }

        for i in 0..n {
            let mi = self.m_red[(i, i)].max(1e-15);
            let ci = self.c_red[(i, i)];
            let ki = self.k_red[(i, i)];
            let fi = f_modal.get(i).copied().unwrap_or(0.0);
            let qdd = (fi - ci * self.q_dot[i] - ki * self.q[i]) / mi;
            self.q_ddot[i] = qdd;
            self.q_dot[i] += qdd * dt;
            self.q[i] += self.q_dot[i] * dt;
        }
    }

    /// Compute the kinetic energy of the flexible body in modal coordinates.
    pub fn kinetic_energy(&self) -> f64 {
        let mut ke = 0.0;
        for i in 0..self.num_modes {
            let mi = self.m_red[(i, i)];
            ke += 0.5 * mi * self.q_dot[i] * self.q_dot[i];
        }
        ke
    }

    /// Compute the potential (strain) energy of the flexible body.
    pub fn potential_energy(&self) -> f64 {
        let mut pe = 0.0;
        for i in 0..self.num_modes {
            let ki = self.k_red[(i, i)];
            pe += 0.5 * ki * self.q[i] * self.q[i];
        }
        pe
    }
}

// ===========================================================================
// Modal force projection
// ===========================================================================

/// Project a distributed force vector onto the modal basis.
///
/// Given a force vector `f` expressed in full FE coordinates, compute the
/// equivalent generalized forces in modal coordinates:
///
/// ```text
/// f_modal = Φᵀ · f
/// ```
///
/// where `Φ` is the modal basis matrix.
pub fn project_force(phi: &DMatrix<f64>, f: &[f64]) -> Vec<f64> {
    let n_modes = phi.ncols();
    let mut f_modal = vec![0.0; n_modes];
    for j in 0..n_modes {
        let mut sum = 0.0;
        for i in 0..phi.nrows().min(f.len()) {
            sum += phi[(i, j)] * f[i];
        }
        f_modal[j] = sum;
    }
    f_modal
}

// ===========================================================================
// Mode selection for flexible bodies
// ===========================================================================

/// Select modes from a full modal basis according to the given criteria.
///
/// This is a convenience wrapper around [`select_modes`] from the `cms` module.
pub fn select_flexible_modes(
    basis: &ModalBasis,
    frequency_cutoff: Option<f64>,
    participation_threshold: Option<f64>,
) -> ModalBasis {
    let selection = ModeSelection {
        frequency_cutoff,
        participation_factor_threshold: participation_threshold,
    };
    select_modes(&basis.eigenvalues, basis.eigenvectors.clone(), selection)
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_math_linalg_dense::DMatrix;

    #[test]
    fn test_modal_superposition_displacement() {
        let phi = DMatrix::from_row_slice(4, 2, &[1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]);
        let m_red = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0]);
        let k_red = DMatrix::from_row_slice(2, 2, &[100.0, 0.0, 0.0, 400.0]);
        let c_red = DMatrix::from_row_slice(2, 2, &[0.5, 0.0, 0.0, 1.0]);
        let mut state = ModalSuperpositionState::from_reduced_matrices(m_red, k_red, c_red);
        state.q = vec![1.0, 2.0];
        let u = state.displacement(&phi);
        assert_eq!(u.len(), 4);
        assert!((u[0] - 1.0).abs() < 1e-12);
        assert!((u[1] - 2.0).abs() < 1e-12);
        assert!(u[2].abs() < 1e-12);
        assert!(u[3].abs() < 1e-12);
    }

    #[test]
    fn test_modal_superposition_energy() {
        let m_red = DMatrix::from_row_slice(2, 2, &[2.0, 0.0, 0.0, 3.0]);
        let k_red = DMatrix::from_row_slice(2, 2, &[200.0, 0.0, 0.0, 800.0]);
        let c_red = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 2.0]);
        let mut state = ModalSuperpositionState::from_reduced_matrices(m_red, k_red, c_red);
        state.q = vec![0.1, 0.2];
        state.q_dot = vec![1.0, -0.5];
        let ke = state.kinetic_energy();
        let expected_ke = 0.5 * 2.0 * 1.0 + 0.5 * 3.0 * 0.25;
        assert!((ke - expected_ke).abs() < 1e-12);
    }

    #[test]
    fn test_modal_step_semi_implicit_euler() {
        let m_red = DMatrix::from_row_slice(1, 1, &[1.0]);
        let k_red = DMatrix::from_row_slice(1, 1, &[100.0]);
        let c_red = DMatrix::from_row_slice(1, 1, &[0.0]);
        let mut state = ModalSuperpositionState::from_reduced_matrices(m_red, k_red, c_red);
        state.q = vec![0.0];
        state.q_dot = vec![0.0];
        state.step_semi_implicit_euler(&[1.0], 0.001);
        assert!(state.q_dot[0] > 0.0);
        assert!(state.q[0] > 0.0);
    }

    #[test]
    fn test_project_force() {
        let phi = DMatrix::from_row_slice(4, 2, &[1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]);
        let f = vec![10.0, 20.0, 0.0, 0.0];
        let f_modal = project_force(&phi, &f);
        assert_eq!(f_modal.len(), 2);
        assert!((f_modal[0] - 10.0).abs() < 1e-12);
        assert!((f_modal[1] - 20.0).abs() < 1e-12);
    }

    #[test]
    fn test_select_modes_preserves_order() {
        let eig = vec![1.0, 4.0, 9.0, 16.0];
        let evec = DMatrix::from_row_slice(
            4,
            4,
            &[
                1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
            ],
        );
        let sel = select_flexible_modes(
            &ModalBasis {
                eigenvalues: eig.clone(),
                eigenvectors: evec.clone(),
            },
            Some(2.5),
            None,
        );
        assert_eq!(sel.eigenvalues.len(), 2);
        assert_eq!(sel.eigenvalues, vec![1.0, 4.0]);
    }
}
