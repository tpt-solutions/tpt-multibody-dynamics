//! Craig-Bampton component mode synthesis (CMS): boundary/interior DOF partitioning,
//! fixed-interface normal modes, constraint modes, reduced mass/stiffness assembly,
//! and mode selection.
//!
//! The Craig-Bampton method reduces a full FE model to a small set of generalized
//! coordinates by projecting the system onto a basis of fixed-interface modes
//! (free-vibration modes with boundary DOFs fixed) and constraint modes (static
//! shapes produced by unit displacements at the boundary).
// Dense numeric kernels are clearest with explicit indexing; the indexed-loop
// lint does not fit this code.
#![allow(clippy::needless_range_loop)]

use tpt_math_linalg_dense::{DMatrix, DVector};

/// Fixed-interface mode shapes plus eigenvalues.
#[derive(Clone, Debug)]
pub struct ModalBasis {
    /// Natural frequencies squared (eigenvalues) of the retained modes.
    pub eigenvalues: Vec<f64>,
    /// Mode-shape matrix (`n_dofs × n_modes`); each column is one mode.
    pub eigenvectors: DMatrix<f64>,
}

/// Criteria for selecting a subset of modes from a full basis.
#[derive(Clone, Debug)]
pub struct ModeSelection {
    /// Keep only modes with frequency below this cutoff (rad/s). `None` disables
    /// the cutoff.
    pub frequency_cutoff: Option<f64>,
    /// Keep only modes whose participation factor exceeds this threshold.
    /// `None` disables the threshold.
    pub participation_factor_threshold: Option<f64>,
}

/// Craig-Bampton reduced model.
///
/// Holds the boundary DOF list and the assembled modal basis together with the
/// projected (reduced) mass and stiffness matrices.
#[derive(Clone, Debug)]
pub struct CraigBampton {
    /// Boundary (interface) DOF indices in the full FE model.
    pub boundary_dofs: Vec<usize>,
    /// Fixed-interface normal modes (`n_dofs × n_fixed_modes`).
    pub fixed_modes: DMatrix<f64>,
    /// Constraint modes (`n_dofs × n_boundary`); each column corresponds to one
    /// unit displacement imposed at the matching boundary DOF.
    pub constraint_modes: DMatrix<f64>,
    /// Reduced mass matrix (`(n_fixed_modes + n_boundary) × (n_fixed_modes + n_boundary)`).
    pub modal_mass: DMatrix<f64>,
    /// Reduced stiffness matrix (`(n_fixed_modes + n_boundary) × (n_fixed_modes + n_boundary)`).
    pub modal_stiffness: DMatrix<f64>,
}

impl CraigBampton {
    /// Build a Craig-Bampton reduced model from a full FE stiffness matrix.
    ///
    /// * `mesh` — provides the total DOF count.
    /// * `fem_stiffness` — full FE stiffness matrix (square, symmetric).
    /// * `fixed_dofs` — boundary/interface DOF indices (remain active in the
    ///   reduced model).
    /// * `num_modes` — number of fixed-interface modes to retain.
    ///
    /// The interior stiffness matrix `K_ii` is solved for its smallest
    /// eigenpairs (with unit mass `M = I`) to obtain the fixed-interface modes.
    /// Constraint modes are computed statically from `K_ii ψ_i = -K_ib`.
    pub fn from_fem(
        _mesh: &tpt_fem_mesh::Mesh,
        fem_stiffness: DMatrix<f64>,
        fixed_dofs: &[usize],
        num_modes: usize,
    ) -> Self {
        let n = fem_stiffness.nrows();
        let n_b = fixed_dofs.len();
        let n_i = n - n_b;

        let interior_dofs: Vec<usize> = (0..n).filter(|d| !fixed_dofs.contains(d)).collect();

        // Partition K into blocks.
        let k_ii = block(&fem_stiffness, &interior_dofs, &interior_dofs);
        let k_ib = block(&fem_stiffness, &interior_dofs, fixed_dofs);

        // Fixed-interface modes: smallest eigenpairs of K_ii (with M_ii = I).
        let basis = if n_i > 0 && num_modes > 0 {
            solve_eigenvalue_solve(k_ii.clone(), DMatrix::from_fn(n_i, n_i, |_, _| 1.0), num_modes)
        } else {
            ModalBasis {
                eigenvalues: Vec::new(),
                eigenvectors: DMatrix::zeros(n_i, 0),
            }
        };

        let fixed_modes = expand_to_full(&basis.eigenvectors, &interior_dofs, n);

        // Constraint modes: static shapes from unit boundary displacement.
        // K_ii ψ_ij = -K_ib[:, j]
        let constraint_modes = if n_i > 0 && n_b > 0 {
            let mut psi_i = DMatrix::zeros(n_i, n_b);
            for j in 0..n_b {
                let rhs: DVector<f64> = DVector::from_vec((0..n_i).map(|i| -k_ib[(i, j)]).collect());
                if let Ok(sol) = k_ii.clone().solve(&rhs) {
                    for i in 0..n_i {
                        psi_i = with_elem(&psi_i, i, j, sol[i]);
                    }
                }
            }
            expand_to_full(&psi_i, &interior_dofs, n)
        } else {
            DMatrix::zeros(n, n_b)
        };

        // Assemble full modal basis: [fixed_modes | constraint_modes_with_boundary_identity]
        let n_modes = basis.eigenvectors.ncols();
        let phi = if n_modes > 0 || n_b > 0 {
            let total_cols = n_modes + n_b;
            let mut phi = DMatrix::zeros(n, total_cols);
            // Fixed-mode columns: [Φ_ii; 0]
            for col in 0..n_modes {
                for i in 0..n_i {
                    phi = with_elem(&phi, interior_dofs[i], col, fixed_modes[(interior_dofs[i], col)]);
                }
            }
            // Constraint-mode columns: [Ψ_ib; I]
            for j in 0..n_b {
                for i in 0..n_i {
                    phi = with_elem(&phi, interior_dofs[i], n_modes + j, constraint_modes[(interior_dofs[i], j)]);
                }
                phi = with_elem(&phi, fixed_dofs[j], n_modes + j, 1.0);
            }
            phi
        } else {
            DMatrix::zeros(n, 0)
        };

        // Project full matrices onto modal basis.
        let (m_red, k_red) = if phi.ncols() > 0 {
            let mt = phi.transpose();
            let m_full = DMatrix::from_fn(n, n, |i, j| if i == j { 1.0 } else { 0.0 });
            let m_red = mt.clone() * m_full * phi.clone();
            let k_red = mt * fem_stiffness * phi;
            (m_red, k_red)
        } else {
            (DMatrix::zeros(0, 0), DMatrix::zeros(0, 0))
        };

        CraigBampton {
            boundary_dofs: fixed_dofs.to_vec(),
            fixed_modes,
            constraint_modes,
            modal_mass: m_red,
            modal_stiffness: k_red,
        }
    }
}

fn with_elem(m: &DMatrix<f64>, row: usize, col: usize, val: f64) -> DMatrix<f64> {
    let nrows = m.nrows();
    let ncols = m.ncols();
    DMatrix::from_fn(nrows, ncols, |i, j| if i == row && j == col { val } else { m[(i, j)] })
}

/// Project the full FE mass and stiffness matrices onto a Craig-Bampton modal
/// basis, returning the reduced `(M_red, K_red)`.
///
/// `M_red = Φᵀ M_full Φ` and `K_red = Φᵀ K_full Φ`.
pub fn reduce_system(
    cb: &CraigBampton,
    m_full: DMatrix<f64>,
    k_full: DMatrix<f64>,
) -> (DMatrix<f64>, DMatrix<f64>) {
    let n = cb.boundary_dofs.len();
    let n_modes = cb.fixed_modes.ncols();
    let total = n_modes + n;

    if total == 0 {
        return (DMatrix::zeros(0, 0), DMatrix::zeros(0, 0));
    }

    let phi = assemble_phi(&cb.fixed_modes, &cb.constraint_modes, &cb.boundary_dofs);
    let mt = phi.transpose();
    let m_red = mt.clone() * m_full * phi.clone();
    let k_red = mt * k_full * phi;
    (m_red, k_red)
}

/// Solve the dense symmetric generalized eigenproblem `K x = λ M x` for the
/// smallest `num_modes` eigenpairs.
///
/// Returns the eigenvalues and the corresponding eigenvectors as columns of a
/// `DMatrix`. The eigenvectors are `M`-orthonormal (`Xᵀ M X = I`).
pub fn solve_eigenvalue_solve(
    k: DMatrix<f64>,
    m: DMatrix<f64>,
    num_modes: usize,
) -> ModalBasis {
    let n = k.nrows();
    if n == 0 || num_modes == 0 {
        return ModalBasis {
            eigenvalues: Vec::new(),
            eigenvectors: DMatrix::zeros(n, 0),
        };
    }

    // Reduce to standard form: M = L Lᵀ, K' = L⁻¹ K L⁻ᵀ, then K' y = λ y.
    let l_opt = cholesky_dense(&m);
    let (k_prime, l_inv_t) = match l_opt {
        Some(ref l) => {
            let l_inv = invert_lower(l);
            let linv_t = l_inv.transpose();
            let k_dense = mat_to_vecvec(&k);
            let linv_dense = mat_to_vecvec(&l_inv);
            let k_prime = matmul(&matmul(&linv_dense, &k_dense), &transpose(&linv_dense));
            let l_inv_t_dense = mat_to_vecvec(&linv_t);
            (k_prime, l_inv_t_dense)
        }
        None => {
            // M not PD — fall back to standard eigenproblem with M = I.
            (mat_to_vecvec(&k), identity_dense(n))
        }
    };

    let (eig, evecs) = jacobi_eig(&k_prime, 200);

    // Sort ascending by eigenvalue.
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&a, &b| {
        eig[a]
            .partial_cmp(&eig[b])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let take = num_modes.min(n);
    let mut sorted_eig = Vec::with_capacity(take);
    let mut sorted_evecs = vec![vec![0.0; n]; take];
    for (new_i, &old_i) in idx.iter().take(take).enumerate() {
        sorted_eig.push(eig[old_i]);
                    for r in 0..n {
                        sorted_evecs[new_i][r] = evecs[r][old_i];
                    }
    }

    // Back-transform eigenvectors: x = L⁻ᵀ y.
    let evecs_dense = if l_opt.is_some() {
        let mut x = vec![vec![0.0; n]; take];
        for i in 0..take {
            let y = solve_lower_dense_raw(&l_inv_t, &sorted_evecs[i]);
            x[i] = y;
        }
        x
    } else {
        sorted_evecs
    };

    // Transpose to column-major DMatrix.
    let mut result = DMatrix::zeros(n, take);
    for col in 0..take {
        for row in 0..n {
            result = with_elem(&result, row, col, evecs_dense[col][row]);
        }
    }

    ModalBasis {
        eigenvalues: sorted_eig,
        eigenvectors: result,
    }
}

/// Select a subset of modes from a full basis according to `selection`.
///
/// * Frequency cutoff: drop modes above the cutoff frequency.
/// * Participation-factor threshold: drop modes whose participation factor
///   (`||mode_shape||₂`) is below the threshold.
///
/// Modes that fail either criterion are discarded; surviving modes are returned
/// in ascending eigenvalue order.
pub fn select_modes(
    eigenvalues: &[f64],
    eigenvectors: DMatrix<f64>,
    selection: ModeSelection,
) -> ModalBasis {
    let n = eigenvectors.nrows();
    let m = eigenvectors.ncols();
    if m == 0 {
        return ModalBasis {
            eigenvalues: Vec::new(),
            eigenvectors: DMatrix::zeros(n, 0),
        };
    }

    let mut kept: Vec<usize> = Vec::new();
    for i in 0..m {
        let omega = eigenvalues[i].sqrt();
        let _freq = omega / (2.0 * std::f64::consts::PI);
        let participation = column_norm(&eigenvectors, i);

        let keep_freq = selection
            .frequency_cutoff
            .map(|cut| omega <= cut)
            .unwrap_or(true);
        let keep_part = selection
            .participation_factor_threshold
            .map(|thr| participation >= thr)
            .unwrap_or(true);

        if keep_freq && keep_part {
            kept.push(i);
        }
    }

    let mut out_eig = Vec::with_capacity(kept.len());
    let mut out_vec = DMatrix::zeros(n, kept.len());
    for (new_i, &old_i) in kept.iter().enumerate() {
        out_eig.push(eigenvalues[old_i]);
        for row in 0..n {
            out_vec = with_elem(&out_vec, row, new_i, eigenvectors[(row, old_i)]);
        }
    }

    ModalBasis {
        eigenvalues: out_eig,
        eigenvectors: out_vec,
    }
}

// ---------------------------------------------------------------------------
// Dense linear-algebra helpers (operating on Vec<Vec<f64>>)
// ---------------------------------------------------------------------------

fn identity_dense(n: usize) -> Vec<Vec<f64>> {
    (0..n).map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect()).collect()
}

fn mat_to_vecvec(m: &DMatrix<f64>) -> Vec<Vec<f64>> {
    let n = m.nrows();
    let p = m.ncols();
    (0..n).map(|i| (0..p).map(|j| m[(i, j)]).collect()).collect()
}

fn vecvec_to_mat(data: &[Vec<f64>]) -> DMatrix<f64> {
    let n = data.len();
    let p = data[0].len();
    let flat: Vec<f64> = data.iter().flat_map(|row| row.iter().cloned()).collect();
    DMatrix::from_row_slice(n, p, &flat)
}

fn transpose(a: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = a.len();
    let p = a[0].len();
    (0..p).map(|j| (0..n).map(|i| a[i][j]).collect()).collect()
}

fn matmul(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = a.len();
    let k = a[0].len();
    let m = b[0].len();
    let mut c = vec![vec![0.0; m]; n];
    for i in 0..n {
        for l in 0..k {
            let a_ik = a[i][l];
            if a_ik == 0.0 {
                continue;
            }
            for j in 0..m {
                c[i][j] += a_ik * b[l][j];
            }
        }
    }
    c
}

fn cholesky_dense(a: &DMatrix<f64>) -> Option<DMatrix<f64>> {
    let n = a.nrows();
    let mut l = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..=i {
            let mut s = a[(i, j)];
            for k in 0..j {
                s -= l[i][k] * l[j][k];
            }
            if i == j {
                if s <= 0.0 {
                    return None;
                }
                l[i][j] = s.sqrt();
            } else {
                l[i][j] = s / l[j][j];
            }
        }
    }
    Some(vecvec_to_mat(&l))
}

fn invert_lower(l: &DMatrix<f64>) -> DMatrix<f64> {
    let n = l.nrows();
    let lv = mat_to_vecvec(l);
    let mut inv = vec![vec![0.0; n]; n];
    for col in 0..n {
        let mut e = vec![0.0; n];
        e[col] = 1.0;
        let x = solve_lower_dense_raw(&lv, &e);
        for row in 0..n {
            inv[row][col] = x[row];
        }
    }
    vecvec_to_mat(&inv)
}

fn solve_lower_dense_raw(l: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = l.len();
    let mut x = vec![0.0; n];
    for i in 0..n {
        let mut s = b[i];
        for j in 0..i {
            s -= l[i][j] * x[j];
        }
        x[i] = s / l[i][i];
    }
    x
}

fn jacobi_eig(a: &[Vec<f64>], max_sweeps: usize) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = a.len();
    let mut a = a.to_vec();
    let mut v = identity_dense(n);

    for _ in 0..max_sweeps {
        let mut changed = false;
        for p in 0..n {
            for q in (p + 1)..n {
                if a[p][q].abs() > 1e-12 {
                    changed = true;
                    let app = a[p][p];
                    let aqq = a[q][q];
                    let apq = a[p][q];
                    let theta = (aqq - app) / (2.0 * apq);
                    let t = if theta >= 0.0 {
                        1.0 / (theta + (1.0 + theta * theta).sqrt())
                    } else {
                        1.0 / (theta - (1.0 + theta * theta).sqrt())
                    };
                    let c = 1.0 / (1.0 + t * t).sqrt();
                    let s = t * c;

                    for i in 0..n {
                        if i != p && i != q {
                            let api = a[p][i];
                            let aqi = a[q][i];
                            a[p][i] = c * api - s * aqi;
                            a[q][i] = s * api + c * aqi;
                            a[i][p] = a[p][i];
                            a[i][q] = a[q][i];
                        }
                    }

                    a[p][p] = c * c * app - 2.0 * c * s * apq + s * s * aqq;
                    a[q][q] = s * s * app + 2.0 * c * s * apq + c * c * aqq;
                    a[p][q] = 0.0;
                    a[q][p] = 0.0;

                    for i in 0..n {
                        let vip = v[i][p];
                        let viq = v[i][q];
                        v[i][p] = c * vip - s * viq;
                        v[i][q] = s * vip + c * viq;
                    }
                }
            }
        }
        if !changed {
            break;
        }
    }

    let eig: Vec<f64> = (0..n).map(|i| a[i][i]).collect();
    (eig, v)
}

fn column_norm(m: &DMatrix<f64>, col: usize) -> f64 {
    let n = m.nrows();
    let mut s = 0.0;
    for i in 0..n {
        s += m[(i, col)] * m[(i, col)];
    }
    s.sqrt()
}

// ---------------------------------------------------------------------------
// Block helpers
// ---------------------------------------------------------------------------

fn block(a: &DMatrix<f64>, rows: &[usize], cols: &[usize]) -> DMatrix<f64> {
    DMatrix::from_fn(rows.len(), cols.len(), |i, j| a[(rows[i], cols[j])])
}

fn expand_to_full(sub: &DMatrix<f64>, dofs: &[usize], n: usize) -> DMatrix<f64> {
    let n_i = sub.nrows();
    let n_cols = sub.ncols();
    let mut full = DMatrix::zeros(n, n_cols);
    for col in 0..n_cols {
        for i in 0..n_i {
            full = with_elem(&full, dofs[i], col, sub[(i, col)]);
        }
    }
    full
}

fn assemble_phi(
    fixed_modes: &DMatrix<f64>,
    constraint_modes: &DMatrix<f64>,
    boundary_dofs: &[usize],
) -> DMatrix<f64> {
    let n = fixed_modes.nrows();
    let n_modes = fixed_modes.ncols();
    let n_b = boundary_dofs.len();
    let total = n_modes + n_b;
    if total == 0 {
        return DMatrix::zeros(n, 0);
    }

    let mut phi = DMatrix::zeros(n, total);
    for col in 0..n_modes {
        for row in 0..n {
            phi = with_elem(&phi, row, col, fixed_modes[(row, col)]);
        }
    }
    for j in 0..n_b {
        for row in 0..n {
            phi = with_elem(&phi, row, n_modes + j, constraint_modes[(row, j)]);
        }
    }
    phi
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cms_fixed_modes_zero_at_boundary() {
        // 4-DOF spring-mass chain: [0]--[1]--[2]--[3]
        // Boundary DOFs: 0, 3 (ends fixed during mode computation)
        // Interior DOFs: 1, 2
        let k = DMatrix::from_row_slice(
            4,
            4,
            &[
                2.0, -1.0, 0.0, 0.0,
                -1.0, 2.0, -1.0, 0.0,
                0.0, -1.0, 2.0, -1.0,
                0.0, 0.0, -1.0, 2.0,
            ],
        );
        let mut mesh_builder = tpt_fem_mesh::MeshBuilder::new();
        mesh_builder.add_node(vec![0.0, 0.0, 0.0]);
        mesh_builder.add_node(vec![1.0, 0.0, 0.0]);
        mesh_builder.add_node(vec![2.0, 0.0, 0.0]);
        mesh_builder.add_node(vec![3.0, 0.0, 0.0]);
        let mesh = mesh_builder.build();

        let fixed_dofs = vec![0usize, 3];
        let cb = CraigBampton::from_fem(&mesh, k, &fixed_dofs, 1);

        // Fixed modes must be zero at boundary DOFs.
        for col in 0..cb.fixed_modes.ncols() {
            for &bdof in &cb.boundary_dofs {
                assert!(
                    cb.fixed_modes[(bdof, col)].abs() < 1e-10,
                    "fixed mode {} non-zero at boundary DOF {}: {}",
                    col,
                    bdof,
                    cb.fixed_modes[(bdof, col)]
                );
            }
        }
    }

    #[test]
    fn test_cms_modal_mass_fixed_block_diagonal() {
        // Same 4-DOF system as above.
        let k = DMatrix::from_row_slice(
            4,
            4,
            &[
                2.0, -1.0, 0.0, 0.0,
                -1.0, 2.0, -1.0, 0.0,
                0.0, -1.0, 2.0, -1.0,
                0.0, 0.0, -1.0, 2.0,
            ],
        );
        let mut mesh_builder = tpt_fem_mesh::MeshBuilder::new();
        mesh_builder.add_node(vec![0.0, 0.0, 0.0]);
        mesh_builder.add_node(vec![1.0, 0.0, 0.0]);
        mesh_builder.add_node(vec![2.0, 0.0, 0.0]);
        mesh_builder.add_node(vec![3.0, 0.0, 0.0]);
        let mesh = mesh_builder.build();

        let fixed_dofs = vec![0usize, 3];
        let cb = CraigBampton::from_fem(&mesh, k, &fixed_dofs, 1);

        let n_modes = cb.fixed_modes.ncols();
        // The top-left block (fixed-mode × fixed-mode) of modal_mass should be
        // diagonal (identity, because fixed modes are M-orthonormal with M = I).
        for i in 0..n_modes {
            for j in 0..n_modes {
                let val = cb.modal_mass[(i, j)];
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (val - expected).abs() < 1e-8,
                    "modal_mass[{},{}] = {}, expected {}",
                    i,
                    j,
                    val,
                    expected
                );
            }
        }
    }

    #[test]
    fn test_select_modes_cutoff() {
        let eig = vec![1.0, 4.0, 9.0, 16.0];
        let evec = DMatrix::from_row_slice(
            4,
            4,
            &[
                1.0, 0.0, 0.0, 0.0,
                0.0, 1.0, 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0,
                0.0, 0.0, 0.0, 1.0,
            ],
        );
        let sel = ModeSelection {
            frequency_cutoff: Some(2.5), // rad/s cutoff -> keep ω² ≤ 6.25
            participation_factor_threshold: None,
        };
        let result = select_modes(&eig, evec, sel);
        assert_eq!(result.eigenvalues.len(), 2);
        assert_eq!(result.eigenvalues, vec![1.0, 4.0]);
    }

    #[test]
    fn test_solve_eigenvalue_known_system() {
        // K = [[2, -1], [-1, 2]], M = I.  Eigenvalues: 1, 3.  Eigenvectors: [1,1]/√2, [1,-1]/√2.
        let k = DMatrix::from_row_slice(2, 2, &[2.0, -1.0, -1.0, 2.0]);
        let m = DMatrix::from_fn(2, 2, |_, _| 1.0);
        let result = solve_eigenvalue_solve(k, m, 2);
        assert_eq!(result.eigenvalues.len(), 2);
        assert!(
            (result.eigenvalues[0] - 1.0).abs() < 1e-8,
            "expected 1.0, got {}",
            result.eigenvalues[0]
        );
        assert!(
            (result.eigenvalues[1] - 3.0).abs() < 1e-8,
            "expected 3.0, got {}",
            result.eigenvalues[1]
        );
        assert_eq!(result.eigenvectors.nrows(), 2);
        assert_eq!(result.eigenvectors.ncols(), 2);
    }
}
