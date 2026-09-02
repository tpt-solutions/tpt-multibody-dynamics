#![doc = include_str!("../../../README.md")]
#![allow(missing_docs)]
#![allow(clippy::needless_range_loop)]
#![forbid(unsafe_code)]

//! Cross-cutting validation harness for `tpt-mbd` (Phase 8).
//!
//! Exercises the complete workspace against analytical benchmarks and
//! regression targets defined in `spec.txt` §6:
//!
//! - Dynamics vs. analytical solutions for 20+ benchmark problems
//!   (pendulum, double pendulum, acrobot, cart-pole, gyroscope).
//! - Constraint satisfaction `||Φ|| < 1e-6` over 10 000 steps for 15+
//!   constrained systems (pendulum, four-bar, slider-crank, Stewart platform).
//! - Craig-Bampton vs. full FE for 10+ flexible-body benchmarks
//!   (cantilever beam, rotating plate, flexible manipulator).
//! - Contact force vs. Hertzian analytical solutions (sphere-sphere,
//!   sphere-plane, cylinder-cylinder).
//! - Friction vs. analytical solutions for block-on-plane, rolling wheel,
//!   brake pad; correct stick-slip transition.
//! - Performance at 10/100/1000 DOFs targeting 1 kHz real-time for
//!   <100 DOFs.
//! - Numerical stability for stiff systems over 100 000+ time steps.
//! - Regression tracking of energy drift, constraint violation, and
//!   compute time across code changes.
//!
//! Each test records its result into a small regression summary printed at
//! the end of the test binary, enabling `cargo test` to act as the
//! continuous-regression harness.

extern crate alloc;

mod benchmarks;
mod constraint_satisfaction;
mod contact_validation;
mod flexible_validation;
mod friction_validation;
mod performance;
mod regression;
mod stability;

pub use benchmarks::run_all as run_dynamics;

/// Render the full Phase 8 regression summary to a string.
pub fn full_summary() -> String {
    let mut s = benchmarks::run_all().render();
    s.push('\n');
    s.push_str(&constraint_satisfaction::run_all().render());
    s.push('\n');
    s.push_str(&contact_validation::run_all().render());
    s.push('\n');
    s.push_str(&friction_validation::run_all().render());
    s.push('\n');
    s.push_str(&flexible_validation::run_all().render());
    s.push('\n');
    s.push_str(&performance::run_all().render());
    s.push('\n');
    s.push_str(&stability::run_all().render());
    s
}
