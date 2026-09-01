//! Rayleigh modal damping: αM + βK projected onto modal coordinates.
//!
//! Rayleigh damping assumes the damping matrix is a linear combination of the
//! mass and stiffness matrices. When projected onto the modal basis, each mode
//! acquires a damping ratio ζ_i = (α + β·ω_i²) / (2·ω_i).

use tpt_math_linalg_dense::DMatrix;

/// Rayleigh damping parameters.
///
/// The physical damping matrix is approximated as `C = α·M + β·K`.
#[derive(Clone, Debug, Copy)]
pub struct RayleighDamping {
    /// Mass-proportional damping coefficient (`α`).
    pub alpha: f64,
    /// Stiffness-proportional damping coefficient (`β`).
    pub beta: f64,
}

impl RayleighDamping {
    /// Build a new Rayleigh-damping specification.
    pub fn new(alpha: f64, beta: f64) -> Self {
        RayleighDamping { alpha, beta }
    }
}

/// Damping ratio associated with a single mode.
#[derive(Clone, Debug, Copy)]
pub struct ModalDampingRatio {
    /// Zero-based index of the mode.
    pub mode_index: usize,
    /// Damping ratio `ζ` (dimensionless, typically 0.01–0.1 for lightly damped structures).
    pub damping_ratio: f64,
}

/// Compute the modal damping matrix `C_red = α·M_red + β·K_red`.
///
/// Both `modal_mass` and `modal_stiffness` must be square and of the same size.
pub fn compute_damping_matrix(
    rayleigh: &RayleighDamping,
    modal_mass: DMatrix<f64>,
    modal_stiffness: DMatrix<f64>,
) -> DMatrix<f64> {
    let n = modal_mass.nrows();
    let p = modal_mass.ncols();
    let alpha = rayleigh.alpha;
    let beta = rayleigh.beta;
    DMatrix::from_fn(n, p, |i, j| alpha * modal_mass[(i, j)] + beta * modal_stiffness[(i, j)])
}

/// Compute the damping ratio for each mode given the Rayleigh coefficients.
///
/// For mode `i` with natural frequency `ω_i` (rad/s), the damping ratio is
/// `ζ_i = (α + β·ω_i²) / (2·ω_i)`.  Returns one [`ModalDampingRatio`] per
/// eigenvalue in the order provided.
pub fn compute_damping_ratios(eigenvalues: &[f64], alpha: f64, beta: f64) -> Vec<ModalDampingRatio> {
    eigenvalues
        .iter()
        .enumerate()
        .map(|(i, &omega_sq)| {
            let omega = omega_sq.sqrt();
            let zeta = (alpha + beta * omega_sq) / (2.0 * omega);
            ModalDampingRatio {
                mode_index: i,
                damping_ratio: zeta,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_damping_ratio_formula() {
        // ω = 2π rad/s  ->  ω² = 4π²
        // α = 0.1, β = 0.01
        // ζ = (0.1 + 0.01 * 4π²) / (2 * 2π)
        let omega = 2.0 * std::f64::consts::PI;
        let eig = vec![omega * omega];
        let ratios = compute_damping_ratios(&eig, 0.1, 0.01);
        assert_eq!(ratios.len(), 1);
        let expected = (0.1 + 0.01 * omega * omega) / (2.0 * omega);
        assert!(
            (ratios[0].damping_ratio - expected).abs() < 1e-12,
            "got {} expected {}",
            ratios[0].damping_ratio,
            expected
        );
        assert_eq!(ratios[0].mode_index, 0);
    }

    #[test]
    fn test_damping_matrix_linear_combination() {
        let m = DMatrix::from_row_slice(2, 2, &[2.0, 0.0, 0.0, 3.0]);
        let k = DMatrix::from_row_slice(2, 2, &[4.0, -2.0, -2.0, 4.0]);
        let rayleigh = RayleighDamping::new(0.1, 0.01);
        let c = compute_damping_matrix(&rayleigh, m.clone(), k.clone());
        let expected = DMatrix::from_row_slice(2, 2, &[0.1 * 2.0 + 0.01 * 4.0, 0.1 * 0.0 + 0.01 * (-2.0), 0.1 * 0.0 + 0.01 * (-2.0), 0.1 * 3.0 + 0.01 * 4.0]);
        for i in 0..2 {
            for j in 0..2 {
                assert!(
                    (c[(i, j)] - expected[(i, j)]).abs() < 1e-12,
                    "C[{},{}] = {} expected {}",
                    i,
                    j,
                    c[(i, j)],
                    expected[(i, j)]
                );
            }
        }
    }
}
