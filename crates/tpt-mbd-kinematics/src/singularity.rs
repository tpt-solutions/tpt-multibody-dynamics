//! Singularity detection and workspace analysis.
//!
//! Provides:
//! - Manipulability measure (Yoshikawa's index)
//! - Jacobian condition number
//! - Distance-to-singularity via eigenvalue analysis
//! - Reachable workspace boundary tracing

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

        SingularityAnalysis { manipulability: manip, condition_number: cond, is_singular, singularity_type }
    }
}

/// Workspace boundary point in Cartesian space.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WorkspacePoint {
    pub position: [f64; 3],
    pub orientation: [f64; 4],
}

/// Compute the reachable workspace boundary by sampling joint space.
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
    let orientation = [q.coords.data[0], q.coords.data[1], q.coords.data[2], q.coords.data[3]];
    WorkspacePoint { position, orientation }
}
