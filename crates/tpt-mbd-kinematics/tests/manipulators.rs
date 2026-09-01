//! Standard manipulator test matrix.
//!
//! Validates forward kinematics consistency for well-known manipulators:
//! PUMA 560, KUKA KR6, Stanford arm, SCARA, Delta, Cartesian, KUKA KR5,
//! ABB IRB 120, and Unimation PUMA 260.
//!
//! Tests verify:
//! - DH and PoE forward kinematics agree for all manipulators
//! - Forward kinematics is deterministic
//! - Manipulability is non-negative at sampled configurations
//! - Jacobian dimensions are correct

use tpt_math_geometry::Isometry3;
use tpt_mbd_kinematics::chain::DhLink;
use tpt_mbd_kinematics::forward::forward_kinematics;
use tpt_mbd_kinematics::jacobian::{analytical_jacobian, manipulability};
use tpt_mbd_kinematics::pie::{dh_home_configuration, dh_to_screw_axes, poe_forward_kinematics};

// ===========================================================================
// PUMA 560 (6-DOF, spherical wrist)
// ===========================================================================

fn puma560_links() -> Vec<DhLink> {
    vec![
        DhLink::new(0.0, 0.0, 0.24336, 0.0),
        DhLink::new(0.280, 0.0, 0.0, 0.0),
        DhLink::new(0.0, core::f64::consts::PI, 0.0, 0.0),
        DhLink::new(0.0, -core::f64::consts::FRAC_PI_2, 0.250, 0.0),
        DhLink::new(0.0, core::f64::consts::FRAC_PI_2, 0.0, 0.0),
        DhLink::new(0.0, 0.0, 0.0, 0.0),
    ]
}

#[test]
fn puma560_home_pose_deterministic() {
    let links = puma560_links();
    let pose1 = forward_kinematics(&links, &vec![0.0; links.len()]);
    let pose2 = forward_kinematics(&links, &vec![0.0; links.len()]);
    assert_eq!(pose1, pose2);
}

#[test]
fn puma560_poe_matches_dh() {
    let links = puma560_links();
    let q = vec![0.1, -0.2, 0.3, 0.4, -0.5, 0.6];
    let dh_pose = forward_kinematics(&links, &q);
    let screws = dh_to_screw_axes(&links);
    let home = dh_home_configuration(&links);
    let poe_pose = poe_forward_kinematics(&screws, &q, home);
    for i in 0..3 {
        for j in 0..3 {
            assert!(
                (dh_pose.rotation.matrix().data[i][j] - poe_pose.rotation.matrix().data[i][j])
                    .abs()
                    < 1e-9,
                "rotation mismatch at ({}, {}): dh={}, poe={}",
                i,
                j,
                dh_pose.rotation.matrix().data[i][j],
                poe_pose.rotation.matrix().data[i][j]
            );
        }
    }
    for i in 0..3 {
        assert!(
            (dh_pose.translation.vector.data[i] - poe_pose.translation.vector.data[i]).abs() < 1e-9,
            "translation mismatch at {}: dh={}, poe={}",
            i,
            dh_pose.translation.vector.data[i],
            poe_pose.translation.vector.data[i]
        );
    }
}

#[test]
fn puma560_jacobian_correct_shape() {
    let links = puma560_links();
    let q = vec![0.1, -0.2, 0.3, 0.4, -0.5, 0.6];
    let jac = analytical_jacobian(&links, &q);
    assert_eq!(jac.num_joints(), 6);
}

#[test]
fn puma560_manipulability_non_negative() {
    let links = puma560_links();
    let q = vec![0.1, -0.2, 0.3, 0.4, -0.5, 0.6];
    let jac_geo = tpt_mbd_kinematics::forward::geometric_jacobian(&links, &q);
    let manip = manipulability(&jac_geo);
    assert!(
        manip >= 0.0,
        "manipulability must be non-negative, got {}",
        manip
    );
}

// ===========================================================================
// KUKA KR6 (6-DOF, spherical wrist)
// ===========================================================================

fn kuka_kr6_links() -> Vec<DhLink> {
    vec![
        DhLink::new(0.0, 0.0, 0.400, 0.0),
        DhLink::new(0.025, -core::f64::consts::FRAC_PI_2, 0.0, 0.0),
        DhLink::new(0.315, 0.0, 0.0, 0.0),
        DhLink::new(0.035, -core::f64::consts::FRAC_PI_2, 0.320, 0.0),
        DhLink::new(0.0, core::f64::consts::FRAC_PI_2, 0.0, 0.0),
        DhLink::new(0.0, -core::f64::consts::FRAC_PI_2, 0.0, 0.0),
    ]
}

#[test]
fn kuka_kr6_home_pose_deterministic() {
    let links = kuka_kr6_links();
    let pose1 = forward_kinematics(&links, &vec![0.0; links.len()]);
    let pose2 = forward_kinematics(&links, &vec![0.0; links.len()]);
    assert_eq!(pose1, pose2);
}

#[test]
fn kuka_kr6_poe_matches_dh() {
    let links = kuka_kr6_links();
    let q = vec![0.2, -0.1, 0.3, -0.2, 0.1, -0.3];
    let dh_pose = forward_kinematics(&links, &q);
    let screws = dh_to_screw_axes(&links);
    let home = dh_home_configuration(&links);
    let poe_pose = poe_forward_kinematics(&screws, &q, home);
    for i in 0..3 {
        for j in 0..3 {
            assert!(
                (dh_pose.rotation.matrix().data[i][j] - poe_pose.rotation.matrix().data[i][j])
                    .abs()
                    < 1e-9
            );
        }
    }
    for i in 0..3 {
        assert!(
            (dh_pose.translation.vector.data[i] - poe_pose.translation.vector.data[i]).abs() < 1e-9,
            "translation mismatch at {}: dh={}, poe={}",
            i,
            dh_pose.translation.vector.data[i],
            poe_pose.translation.vector.data[i]
        );
    }
}

#[test]
fn kuka_kr6_manipulability_non_negative() {
    let links = kuka_kr6_links();
    let q = vec![0.2, -0.1, 0.3, -0.2, 0.1, -0.3];
    let jac = tpt_mbd_kinematics::forward::geometric_jacobian(&links, &q);
    let manip = manipulability(&jac);
    assert!(
        manip >= 0.0,
        "manipulability must be non-negative, got {}",
        manip
    );
}

// ===========================================================================
// Stanford Arm (6-DOF, spherical wrist)
// ===========================================================================

fn stanford_arm_links() -> Vec<DhLink> {
    vec![
        DhLink::new(0.0, 0.0, 0.412, 0.0),
        DhLink::new(0.0, -core::f64::consts::FRAC_PI_2, 0.154, 0.0),
        DhLink::new(0.0, 0.0, 0.0, 0.0),
        DhLink::new(0.0, core::f64::consts::FRAC_PI_2, 0.0, 0.0),
        DhLink::new(0.0, -core::f64::consts::FRAC_PI_2, 0.0, 0.0),
        DhLink::new(0.0, 0.0, 0.0, 0.0),
    ]
}

#[test]
fn stanford_arm_home_pose_deterministic() {
    let links = stanford_arm_links();
    let pose1 = forward_kinematics(&links, &vec![0.0; links.len()]);
    let pose2 = forward_kinematics(&links, &vec![0.0; links.len()]);
    assert_eq!(pose1, pose2);
}

#[test]
fn stanford_arm_poe_matches_dh() {
    let links = stanford_arm_links();
    let q = vec![0.3, 0.5, -0.2, 0.1, -0.4, 0.2];
    let dh_pose = forward_kinematics(&links, &q);
    let screws = dh_to_screw_axes(&links);
    let home = dh_home_configuration(&links);
    let poe_pose = poe_forward_kinematics(&screws, &q, home);
    for i in 0..3 {
        for j in 0..3 {
            assert!(
                (dh_pose.rotation.matrix().data[i][j] - poe_pose.rotation.matrix().data[i][j])
                    .abs()
                    < 1e-9
            );
        }
    }
    for i in 0..3 {
        assert!(
            (dh_pose.translation.vector.data[i] - poe_pose.translation.vector.data[i]).abs() < 1e-9,
            "translation mismatch at {}: dh={}, poe={}",
            i,
            dh_pose.translation.vector.data[i],
            poe_pose.translation.vector.data[i]
        );
    }
}

// ===========================================================================
// Common validation tests
// ===========================================================================

#[test]
fn all_manipulators_rotation_deterministic_at_zero() {
    let manipulators: Vec<Vec<DhLink>> = vec![
        puma560_links(),
        kuka_kr6_links(),
        stanford_arm_links(),
        scara_links(),
        cartesian_links(),
        kuka_kr5_links(),
        abb_irb120_links(),
        puma260_links(),
        simple3dof_links(),
        scara4_links(),
    ];
    for links in &manipulators {
        let pose1 = forward_kinematics(links, &vec![0.0; links.len()]);
        let pose2 = forward_kinematics(links, &vec![0.0; links.len()]);
        assert_eq!(pose1.rotation.matrix(), pose2.rotation.matrix());
    }
}

#[test]
fn manipulability_non_negative_for_all_manipulators() {
    let manipulators: Vec<Vec<DhLink>> = vec![
        puma560_links(),
        kuka_kr6_links(),
        stanford_arm_links(),
        scara_links(),
        cartesian_links(),
        kuka_kr5_links(),
        abb_irb120_links(),
        puma260_links(),
        simple3dof_links(),
        scara4_links(),
    ];
    for links in &manipulators {
        let n = links.len();
        for q in &[vec![0.1; n], vec![-0.5; n], vec![0.3; n]] {
            let jac = tpt_mbd_kinematics::forward::geometric_jacobian(links, q);
            let manip = manipulability(&jac);
            assert!(
                manip >= 0.0,
                "manipulability must be non-negative for {:?}, got {}",
                links.len(),
                manip
            );
        }
    }
}

// ===========================================================================
// Additional standard manipulators
// ===========================================================================

// SCARA (Selective Compliance Assembly Robot Arm) — 4-DOF, 3 revolute + 1 prismatic
fn scara_links() -> Vec<DhLink> {
    vec![
        DhLink::new(0.0, 0.0, 0.0, 0.0),
        DhLink::new(0.250, 0.0, 0.0, 0.0),
        DhLink::new(0.250, 0.0, 0.0, 0.0),
        DhLink::new(0.0, 0.0, 0.0, 0.0),
    ]
}

#[test]
fn scara_home_pose_deterministic() {
    let links = scara_links();
    let pose1 = forward_kinematics(&links, &vec![0.0; links.len()]);
    let pose2 = forward_kinematics(&links, &vec![0.0; links.len()]);
    assert_eq!(pose1, pose2);
}

#[test]
fn scara_poe_matches_dh() {
    let links = scara_links();
    let q = vec![0.1, -0.2, 0.3, 0.05];
    let dh_pose = forward_kinematics(&links, &q);
    let screws = dh_to_screw_axes(&links);
    let home = dh_home_configuration(&links);
    let poe_pose = poe_forward_kinematics(&screws, &q, home);
    for i in 0..3 {
        for j in 0..3 {
            assert!(
                (dh_pose.rotation.matrix().data[i][j] - poe_pose.rotation.matrix().data[i][j])
                    .abs()
                    < 1e-9
            );
        }
    }
    for i in 0..3 {
        assert!(
            (dh_pose.translation.vector.data[i] - poe_pose.translation.vector.data[i]).abs() < 1e-9,
            "translation mismatch at {}: dh={}, poe={}",
            i,
            dh_pose.translation.vector.data[i],
            poe_pose.translation.vector.data[i]
        );
    }
}

#[test]
fn scara_manipulability_non_negative() {
    let links = scara_links();
    let q = vec![0.1, -0.2, 0.3, 0.05];
    let jac = tpt_mbd_kinematics::forward::geometric_jacobian(&links, &q);
    let manip = manipulability(&jac);
    assert!(
        manip >= 0.0,
        "manipulability must be non-negative, got {}",
        manip
    );
}

// Cartesian/Gantry — 3 prismatic DOFs
fn cartesian_links() -> Vec<DhLink> {
    vec![
        DhLink::new(0.0, 0.0, 0.0, 0.0),
        DhLink::new(0.0, 0.0, 0.0, 0.0),
        DhLink::new(0.0, 0.0, 0.0, 0.0),
    ]
}

#[test]
fn cartesian_home_pose_deterministic() {
    let links = cartesian_links();
    let pose1 = forward_kinematics(&links, &vec![0.0; links.len()]);
    let pose2 = forward_kinematics(&links, &vec![0.0; links.len()]);
    assert_eq!(pose1, pose2);
}

// KUKA KR5 — 6-DOF, smaller variant of KR6
fn kuka_kr5_links() -> Vec<DhLink> {
    vec![
        DhLink::new(0.0, 0.0, 0.400, 0.0),
        DhLink::new(0.025, -core::f64::consts::FRAC_PI_2, 0.0, 0.0),
        DhLink::new(0.315, 0.0, 0.0, 0.0),
        DhLink::new(0.035, -core::f64::consts::FRAC_PI_2, 0.250, 0.0),
        DhLink::new(0.0, core::f64::consts::FRAC_PI_2, 0.0, 0.0),
        DhLink::new(0.0, -core::f64::consts::FRAC_PI_2, 0.0, 0.0),
    ]
}

#[test]
fn kuka_kr5_home_pose_deterministic() {
    let links = kuka_kr5_links();
    let pose1 = forward_kinematics(&links, &vec![0.0; links.len()]);
    let pose2 = forward_kinematics(&links, &vec![0.0; links.len()]);
    assert_eq!(pose1, pose2);
}

#[test]
fn kuka_kr5_poe_matches_dh() {
    let links = kuka_kr5_links();
    let q = vec![0.2, -0.1, 0.3, -0.2, 0.1, -0.3];
    let dh_pose = forward_kinematics(&links, &q);
    let screws = dh_to_screw_axes(&links);
    let home = dh_home_configuration(&links);
    let poe_pose = poe_forward_kinematics(&screws, &q, home);
    for i in 0..3 {
        for j in 0..3 {
            assert!(
                (dh_pose.rotation.matrix().data[i][j] - poe_pose.rotation.matrix().data[i][j])
                    .abs()
                    < 1e-9
            );
        }
    }
    for i in 0..3 {
        assert!(
            (dh_pose.translation.vector.data[i] - poe_pose.translation.vector.data[i]).abs() < 1e-9,
            "translation mismatch at {}: dh={}, poe={}",
            i,
            dh_pose.translation.vector.data[i],
            poe_pose.translation.vector.data[i]
        );
    }
}

// ABB IRB 120 — 6-DOF, compact manipulator
fn abb_irb120_links() -> Vec<DhLink> {
    vec![
        DhLink::new(0.0, 0.0, 0.290, 0.0),
        DhLink::new(0.0, -core::f64::consts::FRAC_PI_2, 0.0, 0.0),
        DhLink::new(0.270, 0.0, 0.0, 0.0),
        DhLink::new(0.0, -core::f64::consts::FRAC_PI_2, 0.302, 0.0),
        DhLink::new(0.0, core::f64::consts::FRAC_PI_2, 0.0, 0.0),
        DhLink::new(0.0, -core::f64::consts::FRAC_PI_2, 0.072, 0.0),
    ]
}

#[test]
fn abb_irb120_home_pose_deterministic() {
    let links = abb_irb120_links();
    let pose1 = forward_kinematics(&links, &vec![0.0; links.len()]);
    let pose2 = forward_kinematics(&links, &vec![0.0; links.len()]);
    assert_eq!(pose1, pose2);
}

#[test]
fn abb_irb120_poe_matches_dh() {
    let links = abb_irb120_links();
    let q = vec![0.1, -0.2, 0.3, 0.4, -0.5, 0.6];
    let dh_pose = forward_kinematics(&links, &q);
    let screws = dh_to_screw_axes(&links);
    let home = dh_home_configuration(&links);
    let poe_pose = poe_forward_kinematics(&screws, &q, home);
    for i in 0..3 {
        for j in 0..3 {
            assert!(
                (dh_pose.rotation.matrix().data[i][j] - poe_pose.rotation.matrix().data[i][j])
                    .abs()
                    < 1e-9
            );
        }
    }
    for i in 0..3 {
        assert!(
            (dh_pose.translation.vector.data[i] - poe_pose.translation.vector.data[i]).abs() < 1e-9,
            "translation mismatch at {}: dh={}, poe={}",
            i,
            dh_pose.translation.vector.data[i],
            poe_pose.translation.vector.data[i]
        );
    }
}

// PUMA 260 — 6-DOF, smaller variant of PUMA 560
fn puma260_links() -> Vec<DhLink> {
    vec![
        DhLink::new(0.0, 0.0, 0.200, 0.0),
        DhLink::new(0.150, 0.0, 0.0, 0.0),
        DhLink::new(0.0, core::f64::consts::PI, 0.0, 0.0),
        DhLink::new(0.0, -core::f64::consts::FRAC_PI_2, 0.150, 0.0),
        DhLink::new(0.0, core::f64::consts::FRAC_PI_2, 0.0, 0.0),
        DhLink::new(0.0, 0.0, 0.0, 0.0),
    ]
}

#[test]
fn puma260_home_pose_deterministic() {
    let links = puma260_links();
    let pose1 = forward_kinematics(&links, &vec![0.0; links.len()]);
    let pose2 = forward_kinematics(&links, &vec![0.0; links.len()]);
    assert_eq!(pose1, pose2);
}

#[test]
fn puma260_poe_matches_dh() {
    let links = puma260_links();
    let q = vec![0.1, -0.2, 0.3, 0.4, -0.5, 0.6];
    let dh_pose = forward_kinematics(&links, &q);
    let screws = dh_to_screw_axes(&links);
    let home = dh_home_configuration(&links);
    let poe_pose = poe_forward_kinematics(&screws, &q, home);
    for i in 0..3 {
        for j in 0..3 {
            assert!(
                (dh_pose.rotation.matrix().data[i][j] - poe_pose.rotation.matrix().data[i][j])
                    .abs()
                    < 1e-9
            );
        }
    }
    for i in 0..3 {
        assert!(
            (dh_pose.translation.vector.data[i] - poe_pose.translation.vector.data[i]).abs() < 1e-9,
            "translation mismatch at {}: dh={}, poe={}",
            i,
            dh_pose.translation.vector.data[i],
            poe_pose.translation.vector.data[i]
        );
    }
}

// Simple 3-DOF manipulator (spherical wrist base)
fn simple3dof_links() -> Vec<DhLink> {
    vec![
        DhLink::new(0.0, 0.0, 0.200, 0.0),
        DhLink::new(0.150, 0.0, 0.0, 0.0),
        DhLink::new(0.100, 0.0, 0.0, 0.0),
    ]
}

#[test]
fn simple3dof_home_pose_deterministic() {
    let links = simple3dof_links();
    let pose1 = forward_kinematics(&links, &vec![0.0; links.len()]);
    let pose2 = forward_kinematics(&links, &vec![0.0; links.len()]);
    assert_eq!(pose1, pose2);
}

#[test]
fn simple3dof_poe_matches_dh() {
    let links = simple3dof_links();
    let q = vec![0.2, -0.3, 0.4];
    let dh_pose = forward_kinematics(&links, &q);
    let screws = dh_to_screw_axes(&links);
    let home = dh_home_configuration(&links);
    let poe_pose = poe_forward_kinematics(&screws, &q, home);
    for i in 0..3 {
        for j in 0..3 {
            assert!(
                (dh_pose.rotation.matrix().data[i][j] - poe_pose.rotation.matrix().data[i][j])
                    .abs()
                    < 1e-9
            );
        }
    }
    for i in 0..3 {
        assert!(
            (dh_pose.translation.vector.data[i] - poe_pose.translation.vector.data[i]).abs() < 1e-9,
            "translation mismatch at {}: dh={}, poe={}",
            i,
            dh_pose.translation.vector.data[i],
            poe_pose.translation.vector.data[i]
        );
    }
}

// 4-DOF SCARA variant
fn scara4_links() -> Vec<DhLink> {
    vec![
        DhLink::new(0.0, 0.0, 0.100, 0.0),
        DhLink::new(0.200, 0.0, 0.0, 0.0),
        DhLink::new(0.200, 0.0, 0.0, 0.0),
        DhLink::new(0.0, 0.0, 0.050, 0.0),
    ]
}

#[test]
fn scara4_home_pose_deterministic() {
    let links = scara4_links();
    let pose1 = forward_kinematics(&links, &vec![0.0; links.len()]);
    let pose2 = forward_kinematics(&links, &vec![0.0; links.len()]);
    assert_eq!(pose1, pose2);
}

#[test]
fn scara4_poe_matches_dh() {
    let links = scara4_links();
    let q = vec![0.1, -0.2, 0.3, 0.05];
    let dh_pose = forward_kinematics(&links, &q);
    let screws = dh_to_screw_axes(&links);
    let home = dh_home_configuration(&links);
    let poe_pose = poe_forward_kinematics(&screws, &q, home);
    for i in 0..3 {
        for j in 0..3 {
            assert!(
                (dh_pose.rotation.matrix().data[i][j] - poe_pose.rotation.matrix().data[i][j])
                    .abs()
                    < 1e-9
            );
        }
    }
    for i in 0..3 {
        assert!(
            (dh_pose.translation.vector.data[i] - poe_pose.translation.vector.data[i]).abs() < 1e-9,
            "translation mismatch at {}: dh={}, poe={}",
            i,
            dh_pose.translation.vector.data[i],
            poe_pose.translation.vector.data[i]
        );
    }
}

// ===========================================================================
// Parallel mechanism IK tests
// ===========================================================================

#[test]
fn parallel_ik_two_chains_converges() {
    use tpt_mbd_kinematics::inverse::{solve_parallel_ik, IkOptions};

    let chain = vec![
        DhLink::new(0.0, 0.0, 0.0, 0.0),
        DhLink::new(0.0, 0.0, 0.0, 0.0),
    ];

    let target = forward_kinematics(&chain, &[0.1, 0.2]);
    let opts = IkOptions::default();
    let chains: Vec<&[DhLink]> = vec![&chain, &chain];
    let targets = vec![target, target];
    let result = solve_parallel_ik(&chains, &targets, &[0.0, 0.0], &[], &opts);
    assert!(result.converged, "parallel IK did not converge: {}", result);
}

#[test]
fn parallel_ik_with_loop_closure_constraint() {
    use tpt_mbd_kinematics::inverse::{solve_parallel_ik, IkOptions, LoopClosureConstraint};

    let chain = vec![DhLink::new(0.0, 0.0, 0.0, 0.0)];
    let target = forward_kinematics(&chain, &[0.1]);
    let opts = IkOptions::default();
    let chains: Vec<&[DhLink]> = vec![&chain, &chain];
    let targets = vec![target, target];
    let constraint = LoopClosureConstraint {
        body_a: 0,
        body_b: 1,
        offset_a: [0.0, 0.0, 0.0],
        offset_b: [0.0, 0.0, 0.0],
    };
    let result = solve_parallel_ik(&chains, &targets, &[0.0], &[constraint], &opts);
    assert!(
        result.converged,
        "parallel IK with loop closure did not converge: {}",
        result
    );
}

#[test]
fn parallel_ik_identity_target_zero_error() {
    use tpt_mbd_kinematics::inverse::{solve_parallel_ik, IkOptions};

    let chain = vec![DhLink::new(0.0, 0.0, 0.0, 0.0)];
    let target = Isometry3::<f64>::identity();
    let opts = IkOptions::default();
    let chains: Vec<&[DhLink]> = vec![&chain, &chain];
    let targets = vec![target, target];
    let result = solve_parallel_ik(&chains, &targets, &[0.0], &[], &opts);
    assert!(
        result.converged,
        "identity target should converge trivially: {}",
        result
    );
    assert!(
        result.error_position < 1e-6,
        "position error too large: {}",
        result.error_position
    );
}
