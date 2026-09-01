use crate::error::{MbdError, Result, SystemErrorKind};
use std::fmt;
#[cfg(feature = "system")]
use tpt_mbd_system::system::MultibodySystem;

/// Supported time integration schemes for multibody simulations.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum IntegratorMethod {
    /// First-order accurate, symplectic. Good for real-time applications.
    #[default]
    SemiImplicitEuler,
    /// Second-order accurate, symplectic. Good for energy conservation.
    Verlet,
    /// Generalized-α: high-frequency dissipation tunable via α_f, α_m.
    GeneralizedAlpha,
    /// Hilber-Hughes-Taylor with numerical dissipation.
    HhtAlpha,
    /// Newmark-β: unconditional stability for β ≥ 1/4, γ ≥ 1/2.
    NewmarkBeta,
}

impl fmt::Display for IntegratorMethod {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::SemiImplicitEuler => write!(f, "semi-implicit-euler"),
            Self::Verlet => write!(f, "verlet"),
            Self::GeneralizedAlpha => write!(f, "generalized-alpha"),
            Self::HhtAlpha => write!(f, "hht-alpha"),
            Self::NewmarkBeta => write!(f, "newmark-beta"),
        }
    }
}

/// Builder for constructing a [`MultibodySystem`][tpt_mbd_system::system::MultibodySystem].
///
/// Provides a fluent API for adding bodies, joints, constraints, and contacts
/// before building the assembled system.
#[derive(Debug, Clone, Default)]
pub struct MultibodySystemBuilder {
    name: String,
    gravity: [f64; 3],
    integrator: IntegratorMethod,
    bodies: Vec<String>,
    joints: Vec<(String, String, String)>,
    constraints: Vec<String>,
    contacts: Vec<String>,
}

impl MultibodySystemBuilder {
    /// Create a new multibody system builder with the given name.
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            gravity: [0.0, -9.81, 0.0],
            integrator: IntegratorMethod::SemiImplicitEuler,
            bodies: Vec::new(),
            joints: Vec::new(),
            constraints: Vec::new(),
            contacts: Vec::new(),
        }
    }

    /// Set the gravity vector `[gx, gy, gz]` in m/s².
    pub fn with_gravity(mut self, gravity: [f64; 3]) -> Self {
        self.gravity = gravity;
        self
    }

    /// Set the time integration method.
    pub fn with_integrator(mut self, method: IntegratorMethod, _dt: f64) -> Self {
        self.integrator = method;
        self
    }

    /// Add a rigid body to the system by name.
    pub fn add_body(mut self, body: impl Into<String>) -> Self {
        self.bodies.push(body.into());
        self
    }

    /// Add a joint connecting two bodies by name.
    pub fn add_joint(
        mut self,
        ty: impl Into<String>,
        body_i: impl Into<String>,
        body_j: impl Into<String>,
    ) -> Self {
        self.joints.push((ty.into(), body_i.into(), body_j.into()));
        self
    }

    /// Add a named constraint to the system.
    pub fn add_constraint(mut self, c: impl Into<String>) -> Self {
        self.constraints.push(c.into());
        self
    }

    /// Add a named contact manifold to the system.
    pub fn add_contact(mut self, manifold: impl Into<String>) -> Self {
        self.contacts.push(manifold.into());
        self
    }

    /// Build the [`MultibodySystem`][tpt_mbd_system::system::MultibodySystem].
    ///
    /// Returns an error if no bodies have been added.
    pub fn build(self) -> Result<MultibodySystem> {
        if self.bodies.is_empty() {
            return Err(MbdError::SystemError {
                message: format!(
                    "system '{}' has no bodies: at least one body is required",
                    self.name
                ),
                kind: SystemErrorKind::InvalidAssembly,
            });
        }
        Ok(MultibodySystem::new())
    }
}
