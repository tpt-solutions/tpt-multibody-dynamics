//! Singularity detection and workspace analysis.
//!
//! Provides:
//! - Manipulability measure (Yoshikawa's index)
//! - Jacobian condition number
//! - Distance-to-singularity via eigenvalue analysis
//! - Reachable workspace boundary tracing
//!
//! # Examples
//!
//! ```ignore
//! use tpt_mbd_kinematics::singularity::{SingularityAnalysis, distance_to_singularity};
//! use tpt_mbd_kinematics::chain::DhLink;
//!
//! let links = vec![
//!     DhLink::new(0.0, 0.0, 0.24336, 0.0),
//!     DhLink::new(0.280, 0.0, 0.0, 0.0),
//! ];
//! let analysis = SingularityAnalysis::analyze(&links, &[0.0, 0.0]);
//! let dist = distance_to_singularity(&links, &[0.0, 0.0]);
//! assert!(dist >= 0.0);
//! ```

use tpt_math_geometry::{Isometry3, UnitQuaternion};

use crate::forward::{geometric_jacobian, DhLink};

/// Singularity analysis result for a given configuration.
#[derive(Clone, Debug, PartialEq)]
pub struct SingularityAnalysis {
    pub manipulability: f64,
    pub condition_number: f64,
    pub is_singular: bool,
    pub singularity_type: SingularityType,
}

/// Type of singularity detected.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SingularityType {
    None,
    Boundary,
    Interior,
    Decoupled,
}

impl SingularityAnalysis {
    pub fn analyze(links: &[DhLink], joint_angles: &[f64]) -> Self {
        let jac = geometric_jacobian(links, joint_angles);
        let manip = crate::jacobian::manipulability(&jac);
        let cond = crate::jacobian::jacobian_condition_number(&jac);

        let is_singular = manip < 1e-6;
        let singularity_type = if is_singular {
            SingularityType::Boundary
        } else if manip < 1e-3 {
            SingularityType::Interior
        } else {
            SingularityType::None
        };

        SingularityAnalysis {
            manipulability: manip,
            condition_number: cond,
            is_singular,
            singularity_type,
        }
    }
}

/// Workspace boundary point in Cartesian space.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WorkspacePoint {
    pub position: [f64; 3],
    pub orientation: [f64; 4],
}

/// Compute the reachable workspace boundary by sampling joint space.
#[allow(clippy::needless_range_loop)]
pub fn reachable_workspace_boundary(links: &[DhLink], resolution: usize) -> Vec<WorkspacePoint> {
    let mut boundary = Vec::new();
    let n = links.len();

    if n == 0 {
        return boundary;
    }

    let step = std::f64::consts::TAU / resolution as f64;

    if n == 1 {
        for i in 0..resolution {
            let q = [i as f64 * step];
            let pose = crate::forward::forward_kinematics(links, &q);
            boundary.push(workspace_point_from_isometry(&pose));
        }
    } else if n >= 2 {
        for i in 0..resolution {
            for j in 0..resolution {
                let mut q = vec![0.0; n];
                q[0] = i as f64 * step;
                q[1] = j as f64 * step;
                let pose = crate::forward::forward_kinematics(links, &q);
                boundary.push(workspace_point_from_isometry(&pose));
            }
        }
    }

    boundary
}

fn workspace_point_from_isometry(pose: &Isometry3<f64>) -> WorkspacePoint {
    let t = pose.translation.vector;
    let position = [t.data[0], t.data[1], t.data[2]];
    let uq = UnitQuaternion::from_rotation_matrix(&pose.rotation);
    let q = uq.quaternion();
    let orientation = [
        q.coords.data[0],
        q.coords.data[1],
        q.coords.data[2],
        q.coords.data[3],
    ];
    WorkspacePoint {
        position,
        orientation,
    }
}

/// Compute the distance to the nearest singularity via eigenvalue analysis.
///
/// Computes the eigenvalues of `J·Jᵀ` and returns the square root of the
/// minimum eigenvalue, which equals the minimum singular value of `J`.
/// A distance of zero indicates a singular configuration.
#[allow(clippy::needless_range_loop)]
pub fn distance_to_singularity(links: &[DhLink], joint_angles: &[f64]) -> f64 {
    let jac = geometric_jacobian(links, joint_angles);
    let n = jac.num_joints();
    let m = 6.min(n);

    // Build J·Jᵀ as 6×6
    let mut jjt = [[0.0f64; 6]; 6];
    for row in 0..6 {
        for col in 0..6 {
            let mut sum = 0.0;
            for k in 0..m {
                let j_row_k = if row < 3 {
                    jac.angular_column(k)[row]
                } else {
                    jac.linear_column(k)[row - 3]
                };
                let j_col_k = if col < 3 {
                    jac.angular_column(k)[col]
                } else {
                    jac.linear_column(k)[col - 3]
                };
                sum += j_row_k * j_col_k;
            }
            jjt[row][col] = sum;
        }
    }

    let eigenvalues = eigenvalues_symmetric_6x6(&jjt);
    let mut min_eig = eigenvalues[0];
    for i in 1..6 {
        if eigenvalues[i] < min_eig {
            min_eig = eigenvalues[i];
        }
    }

    if min_eig < 0.0 {
        0.0
    } else {
        min_eig.sqrt()
    }
}

/// Compute all eigenvalues of a symmetric 6×6 matrix using Jacobi rotations.
#[allow(clippy::needless_range_loop)]
fn eigenvalues_symmetric_6x6(a: &[[f64; 6]; 6]) -> [f64; 6] {
    let mut a = *a;
    let mut eigenvalues = [0.0f64; 6];

    for _ in 0..100 {
        // Find largest off-diagonal element
        let mut max_val = 0.0;
        let mut p = 0usize;
        let mut q = 1usize;
        for i in 0..6 {
            for j in (i + 1)..6 {
                if a[i][j].abs() > max_val {
                    max_val = a[i][j].abs();
                    p = i;
                    q = j;
                }
            }
        }

        if max_val < 1e-12 {
            break;
        }

        // Compute rotation angle
        let app = a[p][p];
        let aqq = a[q][q];
        let apq = a[p][q];
        let two = 2.0;
        let phi = 0.5 * ((aqq - app) / apq).atan2(two);
        let c = phi.cos();
        let s = phi.sin();

        // Apply Givens rotation: A' = R^T A R
        for i in 0..6 {
            if i != p && i != q {
                let aip = a[i][p];
                let aiq = a[i][q];
                a[i][p] = c * aip - s * aiq;
                a[p][i] = a[i][p];
                a[i][q] = s * aip + c * aiq;
                a[q][i] = a[i][q];
            }
        }

        let app_new = c * c * app - two * s * c * apq + s * s * aqq;
        let aqq_new = s * s * app + two * s * c * apq + c * c * aqq;
        a[p][p] = app_new;
        a[q][q] = aqq_new;
        a[p][q] = (c * c - s * s) * apq + s * c * (aqq - app);
        a[q][p] = a[p][q];
    }

    for i in 0..6 {
        eigenvalues[i] = a[i][i];
    }

    eigenvalues
}
