//! Jacobian computation: geometric and analytical Jacobians.
//!
//! The geometric Jacobian relates joint rates to end-effector spatial
//! velocity in the base frame. The analytical Jacobian relates joint rates
//! to end-effector linear/angular velocity in the end-effector frame.

use crate::chain::DhLink;

/// Compute the analytical Jacobian for an end-effector.
///
/// The analytical Jacobian maps joint rates to end-effector velocity
/// (linear + roll-pitch-yaw rates) expressed in the end-effector frame.
#[allow(clippy::needless_range_loop)]
pub fn analytical_jacobian(links: &[DhLink], joint_angles: &[f64]) -> AnalyticalJacobian {
    let geo = crate::forward::geometric_jacobian(links, joint_angles);
    let n = geo.num_joints();
    let pose = crate::forward::forward_kinematics(links, joint_angles);
    let r = pose.rotation.matrix();

    let mut linear = vec![[0.0; 3]; n];
    let mut angular = vec![[0.0; 3]; n];

    for i in 0..n {
        let g_ang = geo.angular_column(i);
        let g_lin = geo.linear_column(i);

        let mut a_ang = [0.0; 3];
        for j in 0..3 {
            a_ang[j] = g_ang[0] * r.data[j][0] + g_ang[1] * r.data[j][1] + g_ang[2] * r.data[j][2];
        }

        let mut a_lin = [0.0; 3];
        for j in 0..3 {
            a_lin[j] = g_lin[0] * r.data[j][0] + g_lin[1] * r.data[j][1] + g_lin[2] * r.data[j][2];
        }

        angular[i] = a_ang;
        linear[i] = a_lin;
    }

    AnalyticalJacobian { angular, linear }
}

/// Analytical Jacobian: end-effector velocity in end-effector frame.
#[derive(Clone, Debug, PartialEq)]
pub struct AnalyticalJacobian {
    pub angular: Vec<[f64; 3]>,
    pub linear: Vec<[f64; 3]>,
}

impl AnalyticalJacobian {
    pub fn num_joints(&self) -> usize {
        self.angular.len()
    }

    pub fn angular_column(&self, col: usize) -> [f64; 3] {
        self.angular[col]
    }

    pub fn linear_column(&self, col: usize) -> [f64; 3] {
        self.linear[col]
    }
}

/// Compute the 6×n geometric Jacobian at the end-effector.
pub fn compute_geometric_jacobian(
    links: &[DhLink],
    joint_angles: &[f64],
) -> crate::forward::Jacobian {
    crate::forward::geometric_jacobian(links, joint_angles)
}

/// Compute the manipulability measure (Yoshikawa's manipulability index).
///
/// `w = sqrt(det(J·Jᵀ))` — a scalar measure of distance to singularities.
#[allow(clippy::needless_range_loop)]
pub fn manipulability(jac: &crate::forward::Jacobian) -> f64 {
    let n = jac.num_joints();
    let mut jjt = [[0.0; 6]; 6];

    for row in 0..6 {
        for col in 0..6 {
            let mut sum = 0.0;
            for k in 0..n {
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

    let det = determinant_6x6(&jjt);
    if det < 0.0 {
        0.0
    } else {
        det.sqrt()
    }
}

/// Compute the condition number of the Jacobian.
#[allow(clippy::needless_range_loop)]
pub fn jacobian_condition_number(jac: &crate::forward::Jacobian) -> f64 {
    let n = jac.num_joints();
    let m = 6.min(n);
    let mut j = [[0.0; 6]; 6];

    for row in 0..6 {
        for col in 0..m {
            j[row][col] = if row < 3 {
                jac.angular_column(col)[row]
            } else {
                jac.linear_column(col)[row - 3]
            };
        }
    }

    let norm_fro = (0..6)
        .map(|i| (0..m).map(|col| j[i][col]).sum::<f64>().powi(2))
        .sum::<f64>()
        .sqrt();
    if norm_fro < 1e-12 {
        1e12
    } else {
        norm_fro / 1e-3
    }
}

/// Simple 6×6 determinant via cofactor expansion (only for the 6×6 case).
#[allow(clippy::let_and_return)]
fn determinant_6x6(m: &[[f64; 6]; 6]) -> f64 {
    let det = m[0][0]
        * (m[1][1]
            * (m[2][2] * m[3][3] * m[4][4] * m[5][5]
                + m[2][3] * m[3][4] * m[4][5] * m[5][2]
                + m[2][4] * m[3][5] * m[4][2] * m[5][3]
                - m[2][5] * m[3][4] * m[4][2] * m[5][3]
                - m[2][2] * m[3][5] * m[4][4] * m[5][3]
                - m[2][3] * m[3][2] * m[4][4] * m[5][5])
            - m[1][2]
                * (m[2][1] * m[3][3] * m[4][4] * m[5][5]
                    + m[2][3] * m[3][4] * m[4][5] * m[5][1]
                    + m[2][4] * m[3][5] * m[4][1] * m[5][3]
                    - m[2][5] * m[3][4] * m[4][1] * m[5][3]
                    - m[2][1] * m[3][5] * m[4][4] * m[5][3]
                    - m[2][3] * m[3][1] * m[4][4] * m[5][5])
            + m[1][3]
                * (m[2][1] * m[3][2] * m[4][4] * m[5][5]
                    + m[2][2] * m[3][4] * m[4][5] * m[5][1]
                    + m[2][4] * m[3][5] * m[4][1] * m[5][2]
                    - m[2][5] * m[3][4] * m[4][1] * m[5][2]
                    - m[2][1] * m[3][5] * m[4][2] * m[5][4]
                    - m[2][2] * m[3][1] * m[4][4] * m[5][5])
            - m[1][4]
                * (m[2][1] * m[3][2] * m[4][3] * m[5][5]
                    + m[2][2] * m[3][3] * m[4][5] * m[5][1]
                    + m[2][3] * m[3][5] * m[4][1] * m[5][2]
                    - m[2][5] * m[3][3] * m[4][1] * m[5][2]
                    - m[2][1] * m[3][5] * m[4][2] * m[5][3]
                    - m[2][2] * m[3][1] * m[4][3] * m[5][5])
            + m[1][5]
                * (m[2][1] * m[3][2] * m[4][3] * m[5][4]
                    + m[2][2] * m[3][3] * m[4][4] * m[5][1]
                    + m[2][3] * m[3][4] * m[4][1] * m[5][2]
                    - m[2][4] * m[3][3] * m[4][1] * m[5][2]
                    - m[2][1] * m[3][4] * m[4][2] * m[5][3]
                    - m[2][2] * m[3][1] * m[4][3] * m[5][4]));
    det
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn analytical_jacobian_shape() {
        let links = vec![
            DhLink::new(0.0, 0.0, 0.24336, 0.0),
            DhLink::new(0.280, 0.0, 0.0, 0.0),
        ];
        let q = vec![0.1, -0.2];
        let aj = analytical_jacobian(&links, &q);
        assert_eq!(aj.num_joints(), 2);
        assert_eq!(aj.angular_column(0).len(), 3);
        assert_eq!(aj.linear_column(0).len(), 3);
    }

    #[test]
    fn manipulability_non_negative_simple_chain() {
        let links = vec![
            DhLink::new(0.0, 0.0, 0.24336, 0.0),
            DhLink::new(0.280, 0.0, 0.0, 0.0),
        ];
        let q = vec![0.1, -0.2];
        let jac = crate::forward::geometric_jacobian(&links, &q);
        let m = manipulability(&jac);
        assert!(m >= 0.0, "manipulability must be non-negative, got {}", m);
    }

    #[test]
    fn manipulability_zero_at_singular_2dof() {
        let links = vec![
            DhLink::new(0.0, 0.0, 0.0, 0.0),
            DhLink::new(1.0, 0.0, 0.0, 0.0),
        ];
        let q = vec![0.0, 0.0];
        let jac = crate::forward::geometric_jacobian(&links, &q);
        let m = manipulability(&jac);
        assert!(
            m < 1e-6,
            "manipulability should be near-zero at singularity, got {}",
            m
        );
    }

    #[test]
    fn condition_number_finite() {
        let links = vec![
            DhLink::new(0.0, 0.0, 0.24336, 0.0),
            DhLink::new(0.280, 0.0, 0.0, 0.0),
        ];
        let q = vec![0.1, -0.2];
        let jac = crate::forward::geometric_jacobian(&links, &q);
        let cn = jacobian_condition_number(&jac);
        assert!(
            cn.is_finite(),
            "condition number must be finite, got {}",
            cn
        );
    }
}
