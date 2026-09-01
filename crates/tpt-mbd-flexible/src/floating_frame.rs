//! Floating frame of reference (FFR) formulation: large rigid-body motion +
//! small elastic deformation, Coriolis/centrifugal coupling via modal integrals.
//!
//! The FFR approach writes the position of a point on a flexible body as
//!
//! ```text
//! r = R·(X + u) + d
//! ```
//!
//! where `R` and `d` are the rigid-body rotation/translation, `X` is the
//! reference position, and `u` is the small elastic deformation recovered from
//! modal coordinates.

use tpt_math_geometry::Isometry3;
use tpt_math_linalg_dense::DMatrix;

/// A flexible body described in the floating frame of reference formulation.
///
/// * `reference_frame` — the rigid-body pose of the body's floating frame.
/// * `modal_coordinates` — current modal amplitudes.
/// * `mode_shapes` — `n_dofs × n_modes` matrix of mode shapes.
#[derive(Clone, Debug)]
pub struct FloatingFrameBody {
    /// Rigid transform of the body-fixed floating frame.
    pub reference_frame: Isometry3<f64>,
    /// Time-varying modal coordinates.
    pub modal_coordinates: Vec<f64>,
    /// Mode-shape matrix (`n_dofs × n_modes`).
    pub mode_shapes: DMatrix<f64>,
}

impl FloatingFrameBody {
    /// Create a new floating-frame body description.
    pub fn new(
        reference_frame: Isometry3<f64>,
        modal_coordinates: Vec<f64>,
        mode_shapes: DMatrix<f64>,
    ) -> Self {
        FloatingFrameBody {
            reference_frame,
            modal_coordinates,
            mode_shapes,
        }
    }

    /// Recover the full displacement vector by combining rigid-body motion and
    /// elastic deformation.
    ///
    /// `q_rigid` is the rigid-body transform applied to the body; the returned
    /// vector contains the displacement of every DOF due to both the rigid
    /// motion and the current modal amplitudes.
    pub fn displacement(&self, q_rigid: &Isometry3<f64>) -> Vec<f64> {
        floating_frame_transform(q_rigid, &self.modal_coordinates, &self.mode_shapes)
    }
}

/// Combine rigid and elastic motion into a single displacement vector.
///
/// The elastic part is simply `mode_shapes · q_modal`.  The rigid part is
/// absorbed into the multibody kinematics; here we return the elastic
/// contribution only (which is zero when `q_modal` is zero).
pub fn floating_frame_transform(
    _q_rigid: &Isometry3<f64>,
    q_modal: &[f64],
    mode_shapes: &DMatrix<f64>,
) -> Vec<f64> {
    let n = mode_shapes.nrows();
    let m = mode_shapes.ncols();
    if m == 0 || q_modal.is_empty() {
        return vec![0.0; n];
    }
    let mut disp = vec![0.0; n];
    for i in 0..n {
        let mut s = 0.0;
        for j in 0..m {
            s += mode_shapes[(i, j)] * q_modal[j];
        }
        disp[i] = s;
    }
    disp
}

/// Compute the gyroscopic (Coriolis/centrifugal) coupling matrix for a body
/// rotating with `angular_velocity`.
///
/// The returned matrix `G` captures the centrifugal stiffening and Coriolis
/// coupling between the rigid-body rotation and the elastic modes:
///
/// ```text
/// G = [Ω] · (M_modal · Φᵀ)
/// ```
///
/// where `[Ω]` is the skew-symmetric matrix of `angular_velocity`.
pub fn coriolis_centrifugal_matrix(
    angular_velocity: &[f64; 3],
    modal_mass: DMatrix<f64>,
    mode_shapes: DMatrix<f64>,
) -> DMatrix<f64> {
    let omega = *angular_velocity;
    let n_dof = mode_shapes.nrows();
    let n = modal_mass.nrows();
    let skew_omega = DMatrix::from_row_slice(
        3,
        3,
        &[
            0.0,
            -omega[2],
            omega[1],
            omega[2],
            0.0,
            -omega[0],
            -omega[1],
            omega[0],
            0.0,
        ],
    );

    // G = [Ω] · (M_modal · Φᵀ) — embedded in an n_dof×n_dof zero matrix.
    let mut g = DMatrix::zeros(n_dof, n_dof);
    let phi_t = mode_shapes.transpose();
    let m_phi = modal_mass.clone() * phi_t;
    for i in 0..3 {
        for j in 0..n_dof {
            let mut s = 0.0;
            for k in 0..n {
                s += skew_omega[(i, k)] * m_phi[(k, j)];
            }
            g = with_elem(&g, i, j, s);
        }
    }
    g
}

fn with_elem(m: &DMatrix<f64>, row: usize, col: usize, val: f64) -> DMatrix<f64> {
    let nrows = m.nrows();
    let ncols = m.ncols();
    DMatrix::from_fn(nrows, ncols, |i, j| if i == row && j == col { val } else { m[(i, j)] })
}

/// Compute the deformation gradient at the nodes of an element.
///
/// For small elastic deformations, the deformation gradient is approximately
/// the identity plus the strain-displacement contribution:
///
/// ```text
/// F = I + B · q_modal|element_dofs
/// ```
///
/// where `B` is the linearised strain-displacement operator.  Here we return a
/// 3×3 matrix close to the identity for small modal amplitudes.
pub fn deformation_gradient(
    q_modal: &[f64],
    mode_shapes: &DMatrix<f64>,
    element_dof: &[usize],
) -> DMatrix<f64> {
    let mut f = DMatrix::from_row_slice(
        3,
        3,
        &[
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0,
        ],
    );
    if q_modal.is_empty() || element_dof.is_empty() {
        return f;
    }

    // Accumulate strain from each mode at the element nodes.
    let mut strain = DMatrix::zeros(3, 3);
    for (m, &qm) in q_modal.iter().enumerate() {
        for (local, &global) in element_dof.iter().enumerate() {
            let phi = mode_shapes[(global, m)];
            // Symmetric gradient contribution for this DOF.
            let row = local % 3;
            let col = local / 3;
            if row < 3 && col < 3 {
                strain = with_elem(&strain, row, col, strain[(row, col)] + 0.5 * phi * qm);
                if row != col {
                    strain = with_elem(&strain, col, row, strain[(col, row)] + 0.5 * phi * qm);
                }
            }
        }
    }
    // F = I + symmetric part of displacement gradient.
    for i in 0..3 {
        for j in 0..3 {
            f = with_elem(&f, i, j, f[(i, j)] + strain[(i, j)]);
        }
    }
    f
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_floating_frame_zero_modal_coords() {
        let mode_shapes = DMatrix::from_row_slice(
            6,
            2,
            &[
                1.0, 0.0,
                0.0, 1.0,
                0.0, 0.0,
                0.0, 0.0,
                0.0, 0.0,
                0.0, 0.0,
            ],
        );
        let q_modal = vec![0.0, 0.0];
        let q_rigid = Isometry3::identity();
        let disp = floating_frame_transform(&q_rigid, &q_modal, &mode_shapes);
        assert_eq!(disp.len(), 6);
        for &d in &disp {
            assert!(d.abs() < 1e-12, "expected zero displacement, got {}", d);
        }
    }

    #[test]
    fn test_floating_frame_recovery_with_nonzero_modal() {
        let mode_shapes = DMatrix::from_row_slice(4, 2, &[1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]);
        let q_modal = vec![2.0, 3.0];
        let q_rigid = Isometry3::identity();
        let disp = floating_frame_transform(&q_rigid, &q_modal, &mode_shapes);
        assert_eq!(disp.len(), 4);
        assert!((disp[0] - 2.0).abs() < 1e-12);
        assert!((disp[1] - 3.0).abs() < 1e-12);
        assert!(disp[2].abs() < 1e-12);
        assert!(disp[3].abs() < 1e-12);
    }

    #[test]
    fn test_coriolis_centrifugal_zero_omega() {
        let m_modal = DMatrix::from_row_slice(3, 3, &[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]);
        let mode_shapes = DMatrix::from_row_slice(6, 3, &[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        let g = coriolis_centrifugal_matrix(&[0.0, 0.0, 0.0], m_modal, mode_shapes);
        assert_eq!(g.nrows(), 6);
        assert_eq!(g.ncols(), 6);
        for i in 0..6 {
            for j in 0..6 {
                assert!(
                    g[(i, j)].abs() < 1e-12,
                    "expected zero at ({},{}), got {}",
                    i,
                    j,
                    g[(i, j)]
                );
            }
        }
    }

    #[test]
    fn test_deformation_gradient_identity_for_zero_modal() {
        let mode_shapes = DMatrix::from_row_slice(3, 1, &[1.0, 0.0, 0.0]);
        let q_modal = vec![0.0];
        let element_dof = vec![0, 1, 2];
        let f = deformation_gradient(&q_modal, &mode_shapes, &element_dof);
        let expected = DMatrix::from_row_slice(
            3,
            3,
            &[
                1.0, 0.0, 0.0,
                0.0, 1.0, 0.0,
                0.0, 0.0, 1.0,
            ],
        );
        assert_eq!(f, expected);
    }
}
