use std::fmt;

/// Error category for kinematics operations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KinematicsErrorKind {
    /// The Jacobian lost rank at the current configuration.
    SingularConfiguration,
    /// The IK solver failed to converge within the iteration limit.
    IkNotConverged,
    /// The kinematic chain is structurally invalid.
    InvalidChain,
}

/// Error category for dynamics operations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DynamicsErrorKind {
    /// The equations of motion did not converge.
    SolverNotConverged,
    /// The mass matrix is singular.
    SingularMassMatrix,
    /// Energy drift exceeds the tolerance.
    EnergyDrift,
}

/// Error category for contact operations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContactErrorKind {
    /// Contact penetration exceeds the allowable threshold.
    PenetrationTooLarge,
    /// The contact solver did not converge.
    NoConvergence,
}

/// Error category for flexible-body operations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FlexibleErrorKind {
    /// Modal truncation error exceeds the tolerance.
    ModeTruncation,
    /// The flexible-body mesh is invalid.
    InvalidMesh,
}

/// Error category for system-level operations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SystemErrorKind {
    /// The system has no constraints and cannot be solved.
    Unconstrained,
    /// The system assembly is invalid.
    InvalidAssembly,
}

/// Unified error type for the `tpt-mbd` ecosystem.
///
/// Each top-level variant corresponds to a subsystem and carries a human-readable
/// message together with a typed `kind` discriminator.  `KinematicsError`
/// additionally stores an optional chained source error.
#[derive(Debug, Clone)]
pub enum MbdError {
    /// A kinematics-related failure (FK, IK, Jacobian).
    KinematicsError {
        /// Human-readable description of the failure.
        message: String,
        /// Specific kinematics failure mode.
        kind: KinematicsErrorKind,
        /// Optional chained source error.
        source: Option<Box<MbdError>>,
    },
    /// A dynamics-related failure (forward/inverse dynamics).
    DynamicsError {
        /// Human-readable description of the failure.
        message: String,
        /// Specific dynamics failure mode.
        kind: DynamicsErrorKind,
    },
    /// A contact-mechanics failure.
    ContactError {
        /// Human-readable description of the failure.
        message: String,
        /// Specific contact failure mode.
        kind: ContactErrorKind,
    },
    /// A flexible-body failure.
    FlexibleError {
        /// Human-readable description of the failure.
        message: String,
        /// Specific flexible-body failure mode.
        kind: FlexibleErrorKind,
    },
    /// A system-level failure (assembly, integration).
    SystemError {
        /// Human-readable description of the failure.
        message: String,
        /// Specific system failure mode.
        kind: SystemErrorKind,
    },
}

impl fmt::Display for KinematicsErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::SingularConfiguration => write!(f, "singular configuration"),
            Self::IkNotConverged => write!(f, "IK did not converge"),
            Self::InvalidChain => write!(f, "invalid kinematic chain"),
        }
    }
}

impl fmt::Display for DynamicsErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::SolverNotConverged => write!(f, "dynamics solver did not converge"),
            Self::SingularMassMatrix => write!(f, "singular mass matrix"),
            Self::EnergyDrift => write!(f, "excessive energy drift"),
        }
    }
}

impl fmt::Display for ContactErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::PenetrationTooLarge => write!(f, "penetration too large"),
            Self::NoConvergence => write!(f, "contact solver did not converge"),
        }
    }
}

impl fmt::Display for FlexibleErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ModeTruncation => write!(f, "mode truncation error too large"),
            Self::InvalidMesh => write!(f, "invalid flexible mesh"),
        }
    }
}

impl fmt::Display for SystemErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Unconstrained => write!(f, "system is unconstrained"),
            Self::InvalidAssembly => write!(f, "invalid system assembly"),
        }
    }
}

impl fmt::Display for MbdError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::KinematicsError { message, kind, .. } => {
                write!(f, "KinematicsError({}): {}", kind, message)
            }
            Self::DynamicsError { message, kind } => {
                write!(f, "DynamicsError({}): {}", kind, message)
            }
            Self::ContactError { message, kind } => {
                write!(f, "ContactError({}): {}", kind, message)
            }
            Self::FlexibleError { message, kind } => {
                write!(f, "FlexibleError({}): {}", kind, message)
            }
            Self::SystemError { message, kind } => {
                write!(f, "SystemError({}): {}", kind, message)
            }
        }
    }
}

impl std::error::Error for MbdError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::KinematicsError { source: Some(s), .. } => {
                Some(s.as_ref() as &(dyn std::error::Error + 'static))
            }
            _ => None,
        }
    }
}

impl From<String> for MbdError {
    fn from(message: String) -> Self {
        Self::SystemError {
            message,
            kind: SystemErrorKind::InvalidAssembly,
        }
    }
}

impl From<&'static str> for MbdError {
    fn from(message: &'static str) -> Self {
        Self::SystemError {
            message: message.to_string(),
            kind: SystemErrorKind::InvalidAssembly,
        }
    }
}

/// Convenience result alias used throughout `tpt-mbd`.
pub type Result<T> = std::result::Result<T, MbdError>;
