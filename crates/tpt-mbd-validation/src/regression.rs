//! Lightweight regression-tracking record.
//!
//! Each `RegressionEntry` records a single benchmark name, its pass/fail
//! status, a numeric metric (typically an error or drift), and a free-form
//! description. The full set of entries across all Phase 8 modules forms
//! the workspace regression log.

extern crate alloc;

/// One row in the regression table.
#[derive(Debug, Clone)]
pub struct RegressionEntry {
    /// Benchmark identifier.
    pub name: &'static str,
    /// Pass/fail flag.
    pub passed: bool,
    /// Numeric metric (drift, relative error, compute time in seconds, …).
    pub metric: f64,
    /// Free-form text describing the test.
    pub description: &'static str,
}

impl RegressionEntry {
    /// Construct a new entry.
    pub const fn new(
        name: &'static str,
        passed: bool,
        metric: f64,
        description: &'static str,
    ) -> Self {
        Self {
            name,
            passed,
            metric,
            description,
        }
    }
}

/// Convert seconds to milliseconds for printing.
///
/// (Retained as a small helper; not currently called by the validation
/// harness itself but useful for future log formatting.)
#[allow(dead_code)]
#[inline]
pub fn ms(s: f64) -> f64 {
    s * 1.0e3
}
