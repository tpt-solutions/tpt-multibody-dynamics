//! Drift detection and re-projection for constrained systems.
//!
//! Provides:
//! - [`DriftMonitor`] — monitors constraint drift over time
//! - [`DriftReport`] — drift detection report
//! - [`check`] — evaluate drift and determine if re-projection is needed

#![allow(missing_docs)]

/// Drift detection report.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DriftReport {
    /// Maximum constraint violation across all equations.
    pub max_violation: f64,
    /// Whether re-projection is needed (violation exceeds threshold).
    pub needs_reprojection: bool,
    /// Number of constraint equations checked.
    pub num_constraints: usize,
}

impl DriftReport {
    /// Create a new drift report.
    pub fn new(max_violation: f64, needs_reprojection: bool, num_constraints: usize) -> Self {
        Self {
            max_violation,
            needs_reprojection,
            num_constraints,
        }
    }
}

/// Constraint drift monitor.
///
/// Tracks constraint violation magnitude over time and triggers
/// re-projection when drift exceeds a configurable tolerance.
///
/// # Default tolerance
///
/// Default tolerance is `1e-6`.
#[derive(Clone, Debug, PartialEq)]
pub struct DriftMonitor {
    tolerance: f64,
    history: alloc::vec::Vec<f64>,
    max_history: usize,
}

impl DriftMonitor {
    /// Create a new drift monitor with the given tolerance.
    pub fn new(tolerance: f64) -> Self {
        Self {
            tolerance,
            history: alloc::vec::Vec::new(),
            max_history: 100,
        }
    }

    /// Create a new drift monitor with the default tolerance (1e-6).
    pub fn default_tolerance() -> Self {
        Self::new(1e-6)
    }

    /// Check constraint violations and return a drift report.
    ///
    /// # Arguments
    ///
    /// * `phi` — constraint violation values Φ(q)
    ///
    /// # Returns
    ///
    /// A [`DriftReport`] with the max violation and re-projection flag.
    pub fn check(&mut self, phi: &[f64]) -> DriftReport {
        let max_violation = if phi.is_empty() {
            0.0
        } else {
            phi.iter()
                .map(|x| x.abs())
                .fold(f64::NEG_INFINITY, f64::max)
        };

        let needs_reprojection = max_violation > self.tolerance;
        self.history.push(max_violation);
        if self.history.len() > self.max_history {
            self.history.remove(0);
        }

        DriftReport::new(max_violation, needs_reprojection, phi.len())
    }

    /// Check without recording to history.
    pub fn peek(&self, phi: &[f64]) -> DriftReport {
        let max_violation = if phi.is_empty() {
            0.0
        } else {
            phi.iter()
                .map(|x| x.abs())
                .fold(f64::NEG_INFINITY, f64::max)
        };
        let needs_reprojection = max_violation > self.tolerance;
        DriftReport::new(max_violation, needs_reprojection, phi.len())
    }

    /// Get the tolerance threshold.
    pub fn tolerance(&self) -> f64 {
        self.tolerance
    }

    /// Set the tolerance threshold.
    pub fn set_tolerance(&mut self, tolerance: f64) {
        self.tolerance = tolerance;
    }

    /// Get the maximum historical violation.
    pub fn max_history_violation(&self) -> f64 {
        self.history
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Clear the history.
    pub fn clear_history(&mut self) {
        self.history.clear();
    }
}

impl Default for DriftMonitor {
    fn default() -> Self {
        Self::default_tolerance()
    }
}

/// Check drift for a set of constraint violations.
///
/// Convenience function that creates a [`DriftMonitor`] with the given
/// tolerance and checks the violations.
pub fn check(phi: &[f64], tolerance: f64) -> DriftReport {
    let mut monitor = DriftMonitor::new(tolerance);
    monitor.check(phi)
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;

    #[test]
    fn test_drift_within_tolerance() {
        let mut monitor = DriftMonitor::new(1e-6);
        let phi = vec![1e-7, -5e-7, 1e-8];
        let report = monitor.check(&phi);
        assert!(!report.needs_reprojection);
        assert_eq!(report.max_violation, 5e-7);
    }

    #[test]
    fn test_drift_exceeds_tolerance() {
        let mut monitor = DriftMonitor::new(1e-6);
        let phi = vec![1e-7, 2e-5, -1e-8];
        let report = monitor.check(&phi);
        assert!(report.needs_reprojection);
        assert_eq!(report.max_violation, 2e-5);
    }

    #[test]
    fn test_drift_empty() {
        let mut monitor = DriftMonitor::new(1e-6);
        let phi = vec![];
        let report = monitor.check(&phi);
        assert!(!report.needs_reprojection);
        assert_eq!(report.max_violation, 0.0);
    }

    #[test]
    fn test_drift_default_monitor() {
        let mut monitor = DriftMonitor::default();
        let phi = vec![1e-5];
        let report = monitor.check(&phi);
        assert!(report.needs_reprojection);
    }

    #[test]
    fn test_drift_history() {
        let mut monitor = DriftMonitor::new(1e-6);
        monitor.check(&vec![1e-7]);
        monitor.check(&vec![2e-5]);
        assert_eq!(monitor.max_history_violation(), 2e-5);
    }

    #[test]
    fn test_drift_peek() {
        let monitor = DriftMonitor::new(1e-6);
        let phi = vec![1e-5];
        let report = monitor.peek(&phi);
        assert!(report.needs_reprojection);
    }
}
