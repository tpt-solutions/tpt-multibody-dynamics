#![doc = include_str!("../../../README.md")]
#![warn(missing_docs)]
#![forbid(unsafe_code)]

//! Engineering-grade multibody dynamics for robotics, vehicle dynamics,
//! biomechanics, and mechanism simulation.
//!
//! This umbrella crate re-exports capabilities from the constituent crates.
//! Compiling with no features yields a minimal crate with the unified error
//! type and builder/API stubs:
//!
//! ```rust
//! # #[cfg(feature = "core")] {
//! use tpt_mbd::tpt_mbd_core;
//! # }
//! ```
//!
//! Enable individual feature flags to pull in additional solvers:
//!
//! ```toml
//! [dependencies]
//! tpt-mbd = { version = "0.1", features = ["kinematics", "joints", "contact", "system"] }
//! ```
//!
//! Available features:
//! - `core` — spatial algebra, frames, inertia (always available)
//! - `kinematics` — forward/inverse kinematics, Jacobians
//! - `joints` — joint types, constraint formulation, stabilization
//! - `contact` — collision detection, Hertzian contact, friction, impact
//! - `flexible` — Craig-Bampton CMS, modal superposition
//! - `system` — system assembly, time integration, actuators

/// Unified error type and result alias for the `tpt-mbd` ecosystem.
pub mod error;
/// Builder pattern for constructing [`MultibodySystem`][tpt_mbd_system::system::MultibodySystem].
#[cfg(feature = "system")]
pub mod builder;
/// High-level convenience API wrappers around the constituent crates.
pub mod api;

pub use error::{MbdError, Result};

#[cfg(feature = "core")]
pub use tpt_mbd_core;

#[cfg(feature = "kinematics")]
pub use tpt_mbd_kinematics;

#[cfg(feature = "joints")]
pub use tpt_mbd_joints;

#[cfg(feature = "contact")]
pub use tpt_mbd_contact;

#[cfg(feature = "flexible")]
pub use tpt_mbd_flexible;

#[cfg(feature = "system")]
pub use tpt_mbd_system;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mbderror_display_is_human_readable() {
        let err = MbdError::KinematicsError {
            message: "bad angle".to_string(),
            kind: crate::error::KinematicsErrorKind::SingularConfiguration,
            source: None,
        };
        let s = format!("{}", err);
        assert!(s.contains("KinematicsError"));
        assert!(s.contains("singular configuration"));
        assert!(s.contains("bad angle"));
    }

    #[test]
    fn mbderror_system_error_display() {
        let err: MbdError = "oops".into();
        let s = format!("{}", err);
        assert!(s.contains("SystemError"));
        assert!(s.contains("oops"));
    }

    #[test]
    fn mbderror_source_chain() {
        use std::error::Error;
        let inner = MbdError::KinematicsError {
            message: "inner".to_string(),
            kind: crate::error::KinematicsErrorKind::IkNotConverged,
            source: None,
        };
        let outer = MbdError::KinematicsError {
            message: "outer".to_string(),
            kind: crate::error::KinematicsErrorKind::SingularConfiguration,
            source: Some(Box::new(inner)),
        };
        assert!(outer.source().is_some());
    }

    #[cfg(feature = "system")]
    #[test]
    fn builder_produces_valid_system() {
        let system = MultibodySystemBuilder::new("test_system")
            .add_body("ground")
            .add_body("link1")
            .add_joint("revolute", "ground", "link1")
            .build()
            .expect("builder should succeed with bodies");
        let _ = system;
    }

    #[cfg(feature = "system")]
    #[test]
    fn builder_rejects_empty_system() {
        let err = MultibodySystemBuilder::new("empty")
            .build()
            .expect_err("builder should fail with no bodies");
        match err {
            MbdError::SystemError { kind, .. } => {
                assert_eq!(kind, crate::error::SystemErrorKind::InvalidAssembly);
            }
            _ => panic!("expected SystemError"),
        }
    }

    #[cfg(feature = "kinematics")]
    #[test]
    fn forward_kinematics_matches_chain_forward() {
        use crate::api::forward_kinematics;
        use tpt_mbd_kinematics::chain::{DhChain, DhLink};
        let links = vec![
            DhLink::new(0.0, 0.0, 0.0, 0.0),
            DhLink::new(0.0, 0.0, 0.0, 0.0),
        ];
        let chain = DhChain::new(links);
        let angles = [0.1, -0.2];
        let expected = chain.forward(&angles);
        let got = forward_kinematics(&chain, &angles);
        assert_eq!(expected, got);
    }

    #[cfg(feature = "kinematics")]
    #[test]
    fn forward_kinematics_identity_at_zero() {
        use crate::api::forward_kinematics;
        use tpt_mbd_kinematics::chain::{DhChain, DhLink};
        let links = vec![DhLink::new(0.0, 0.0, 0.0, 0.0)];
        let chain = DhChain::new(links);
        let pose = forward_kinematics(&chain, &[]);
        assert_eq!(pose, DhLink::new(0.0, 0.0, 0.0, 0.0).transform());
    }

    #[cfg(all(feature = "kinematics", feature = "system"))]
    #[test]
    fn simulate_returns_empty_result() {
        use tpt_mbd_system::system::MultibodySystem;
        let system = MultibodySystemBuilder::new("sim_test")
            .add_body("b1")
            .build()
            .unwrap();
        let result = simulate(&system, 1.0, 0.01);
        assert!(result.times.is_empty());
        assert!(result.states.is_empty());
    }

    #[cfg(feature = "core")]
    #[test]
    fn core_feature_compiles() {
        let _ = tpt_mbd_core::frame::Frame::identity();
    }

    #[cfg(feature = "kinematics")]
    #[test]
    fn kinematics_feature_compiles() {
        use tpt_mbd_kinematics::chain::DhLink;
        let _ = DhLink::new(0.0, 0.0, 0.0, 0.0);
    }

    #[cfg(feature = "joints")]
    #[test]
    fn joints_feature_compiles() {
        use tpt_mbd_joints::joint::JointDof;
        let _ = JointDof::Revolute;
    }

    #[cfg(feature = "contact")]
    #[test]
    fn contact_feature_compiles() {
        use tpt_mbd_contact::contact::Contact;
        let _ = Contact;
    }

    #[cfg(feature = "flexible")]
    #[test]
    fn flexible_feature_compiles() {
        use tpt_mbd_flexible::cms::Cms;
        let _ = Cms;
    }

    #[cfg(feature = "system")]
    #[test]
    fn system_feature_compiles() {
        use tpt_mbd_system::system::MultibodySystem;
        let _ = MultibodySystem;
    }
}
