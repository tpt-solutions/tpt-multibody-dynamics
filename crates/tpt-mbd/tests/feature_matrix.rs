//! Feature-matrix tests for the `tpt-mbd` umbrella crate.
//!
//! Verifies that each feature combination compiles independently and that the
//! expected types are available under the correct feature flags.

#[cfg(feature = "core")]
#[test]
fn core_only_exports_core_types() {
    let _ = tpt_mbd_core::frame::Frame::identity();
}

#[cfg(all(feature = "core", feature = "kinematics"))]
#[test]
fn core_plus_kinematics_exports_chain() {
    use tpt_mbd_kinematics::chain::DhLink;
    let _ = DhLink::new(0.0, 0.0, 0.0, 0.0);
}

#[cfg(all(feature = "core", feature = "joints"))]
#[test]
fn core_plus_joints_exports_joint_dof() {
    use tpt_mbd_joints::joint::JointDof;
    let _ = JointDof::Revolute;
}

#[cfg(all(feature = "core", feature = "contact"))]
#[test]
fn core_plus_contact_exports_hertzian() {
    use tpt_mbd_contact::contact::{ContactParams, HertzianContact};
    let _ = HertzianContact::new(ContactParams::default());
}

#[cfg(all(feature = "core", feature = "flexible"))]
#[test]
fn core_plus_flexible_exports_craig_bampton() {
    use std::any::TypeId;
    use tpt_mbd_flexible::cms::CraigBampton;
    let _ = TypeId::of::<CraigBampton>();
}

#[cfg(all(feature = "core", feature = "system"))]
#[test]
fn core_plus_system_exports_multibody_system() {
    use tpt_mbd_system::system::MultibodySystem;
    let _ = MultibodySystem::new();
}

#[cfg(all(
    feature = "core",
    feature = "kinematics",
    feature = "joints",
    feature = "contact",
    feature = "flexible",
    feature = "system"
))]
#[test]
fn all_features_compile_together() {
    use std::any::TypeId;
    use tpt_mbd_contact::contact::{ContactParams, HertzianContact};
    use tpt_mbd_core::frame::Frame;
    use tpt_mbd_flexible::cms::CraigBampton;
    use tpt_mbd_joints::joint::JointDof;
    use tpt_mbd_kinematics::chain::DhLink;
    use tpt_mbd_system::system::MultibodySystem;

    let _ = Frame::identity();
    let _ = DhLink::new(0.0, 0.0, 0.0, 0.0);
    let _ = JointDof::Revolute;
    let _ = HertzianContact::new(ContactParams::default());
    let _ = TypeId::of::<CraigBampton>();
    let _ = MultibodySystem::new();
}
