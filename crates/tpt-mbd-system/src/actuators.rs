//! Actuator models: ideal, DC motor, hydraulic, Hill-type muscle.
//!
//! Provides force-producing elements that convert control inputs (voltage,
//! pressure, activation signal) into generalized forces.

// ===========================================================================
// Actuator type enumeration
// ===========================================================================

/// Supported actuator families.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ActuatorType {
    /// Applies a prescribed force/torque with no dynamics.
    Ideal,
    /// brushed/brushless DC motor with electrical dynamics.
    DcMotor,
    /// Hydraulic actuator with bulk modulus compliance.
    Hydraulic,
    /// Hill-type skeletal muscle with activation dynamics.
    HillMuscle,
}

// ===========================================================================
// Ideal actuator
// ===========================================================================

/// An ideal actuator that applies a prescribed force or torque directly.
#[derive(Clone, Copy, Debug, Default)]
pub struct IdealActuator {
    /// Prescribed force magnitude (N) or torque (N·m).
    pub force: f64,
}

impl IdealActuator {
    /// Create an ideal actuator applying the given force/torque.
    pub fn new(force: f64) -> Self {
        Self { force }
    }
}

// ===========================================================================
// DC motor
// ===========================================================================

/// Brushed / brushless DC motor with armature electrical dynamics.
#[derive(Clone, Copy, Debug, Default)]
pub struct DcMotor {
    /// Torque constant (N·m/A).
    pub torque_constant: f64,
    /// Armature resistance (Ω).
    pub resistance: f64,
    /// Armature inductance (H).
    pub inductance: f64,
    /// Rotor inertia (kg·m²).
    pub inertia: f64,
    /// Viscous friction coefficient (N·m·s/rad).
    pub viscous_friction: f64,
}

impl DcMotor {
    /// Create a new DC motor with the given parameters.
    pub fn new(
        torque_constant: f64,
        resistance: f64,
        inductance: f64,
        inertia: f64,
        viscous_friction: f64,
    ) -> Self {
        Self {
            torque_constant,
            resistance,
            inductance,
            inertia,
            viscous_friction,
        }
    }

    /// Advance the electrical and mechanical dynamics by one time step.
    ///
    /// Returns `(torque, current)` after applying `voltage` with the motor
    /// spinning at `angular_velocity`.
    pub fn dynamics(&self, voltage: f64, angular_velocity: f64, dt: f64) -> (f64, f64) {
        let back_emf = self.torque_constant * angular_velocity;
        let di_dt = (voltage - back_emf - self.torque_constant * 0.0) / self.inductance;
        let current = (di_dt * dt).max(0.0);
        let torque = self.torque_constant * current - self.viscous_friction * angular_velocity;
        (torque.max(0.0), current)
    }
}

// ===========================================================================
// Hydraulic actuator (stub)
// ===========================================================================

/// Hydraulic actuator with fluid compressibility (stub — returns zero force).
#[derive(Clone, Copy, Debug, Default)]
pub struct Hydraulic {
    /// Peak force capacity (N).
    pub max_force: f64,
    /// Bulk modulus of the fluid (Pa).
    pub bulk_modulus: f64,
    /// Piston area (m²).
    pub piston_area: f64,
}

impl Hydraulic {
    /// Create a new hydraulic actuator.
    pub fn new(max_force: f64, bulk_modulus: f64, piston_area: f64) -> Self {
        Self {
            max_force,
            bulk_modulus,
            piston_area,
        }
    }

    /// Evaluate force (stub — returns zero).
    pub fn force(&self, _pressure: f64) -> f64 {
        0.0
    }
}

// ===========================================================================
// Hill-type muscle
// ===========================================================================

/// Hill-type skeletal muscle with force-length-velocity properties.
#[derive(Clone, Copy, Debug, Default)]
pub struct HillMuscle {
    /// Maximum isometric force (N).
    pub max_isometric_force: f64,
    /// Optimal fiber length (m).
    pub optimal_fiber_length: f64,
    /// Maximum contraction velocity (m/s).
    pub contraction_velocity: f64,
    /// Current activation level [0, 1].
    pub activation: f64,
}

impl HillMuscle {
    /// Create a new Hill muscle model.
    pub fn new(
        max_isometric_force: f64,
        optimal_fiber_length: f64,
        contraction_velocity: f64,
    ) -> Self {
        Self {
            max_isometric_force,
            optimal_fiber_length,
            contraction_velocity,
            activation: 0.0,
        }
    }

    /// Set the activation level (0 = relaxed, 1 = fully activated).
    pub fn set_activation(&mut self, activation: f64) {
        self.activation = activation.clamp(0.0, 1.0);
    }

    /// Evaluate active muscle force for the given fiber length and velocity.
    ///
    /// Uses the classic Hill-type force-length-velocity relationship:
    /// `F = F_max * activation * f_l(l) * f_v(v)`.
    pub fn force(&self, length: f64, velocity: f64) -> f64 {
        let fl = self.force_length(length);
        let fv = self.force_velocity(velocity);
        self.max_isometric_force * self.activation * fl * fv
    }

    /// Force-length scaling factor (Gaussian-like bell curve centred on optimal length).
    fn force_length(&self, length: f64) -> f64 {
        let l_rel = length / self.optimal_fiber_length.max(1e-12);
        let width = 0.5;
        let diff = (l_rel - 1.0) / width;
        (-diff * diff).exp()
    }

    /// Force-velocity scaling factor (Hill's hyperbolic decay).
    fn force_velocity(&self, velocity: f64) -> f64 {
        let v_rel = velocity / self.contraction_velocity.max(1e-12);
        if v_rel >= 0.0 {
            1.0 - v_rel
        } else {
            (1.0 + v_rel).max(0.0)
        }
    }
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dc_motor_step_response() {
        let motor = DcMotor::new(0.1, 1.0, 0.01, 0.001, 0.001);
        let dt = 0.001f64;
        let mut torque = 0.0f64;
        let mut current = 0.0f64;

        for _ in 0..1000 {
            let (t, c) = motor.dynamics(12.0, 0.0, dt);
            torque = t;
            current = c;
        }

        assert!(
            torque > 0.0,
            "torque should be positive after step, got {}",
            torque
        );
        assert!(
            current > 0.0,
            "current should be positive after step, got {}",
            current
        );
        assert!(
            current < 12.0 / 1.0 + 1.0,
            "current should be below steady-state, got {}",
            current,
        );
    }

    #[test]
    fn hill_muscle_force_at_optimal_length() {
        let mut muscle = HillMuscle::new(1000.0, 0.1, 0.5);
        muscle.set_activation(1.0);
        let f = muscle.force(0.1, 0.0);
        assert!(
            (f - 1000.0).abs() < 1.0,
            "force at optimal length should be ~F_max, got {}",
            f,
        );
    }

    #[test]
    fn ideal_actuator_force() {
        let act = IdealActuator::new(50.0);
        assert_eq!(act.force, 50.0);
    }
}
