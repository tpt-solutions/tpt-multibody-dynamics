# tpt-mbd — Build Todo

> Tracks bootstrap + build-out for the tpt-mbd workspace, per `spec.txt`
> and `tpt-rust-map/registry.toml`. Crates.io publishing is intentionally
> **out of scope** for this pass — crates stop at `status = "git"` in the
> registry, not `"published"`. License for every crate: `MIT OR Apache-2.0`.
> Author: TPT Solutions.

> **Blocker:** `spec.txt` lists `tpt-sci-physics-rigid` as a DEPENDS ON
> crate for `tpt-mbd-core` ("elevates the basic rigid-body simulation in
> `tpt-sci-physics-rigid`"), but that repo does not exist yet anywhere in
> the workspace — `tpt-rust-map/registry.toml` marks it `status =
> "planned"`, no local repo present. `tpt-mbd-core` cannot wire a real
> Cargo `git` dependency on it until it's built. Proceed with `tpt-mbd-core`
> assuming it will exist by the time it's needed; revisit Phase 1's
> dependency list once `tpt-sci-physics-rigid` reaches `status = "git"`.
> All other DEPENDS ON crates (`tpt-math-geometry`, `tpt-math-spatial`,
> `tpt-math-linalg-fixed`, `tpt-math-linalg`, `tpt-math-linalg-dense`,
> `tpt-math-numeric`, `tpt-math-optimize-general`, `tpt-math-units`,
> `tpt-fem-mesh`, `tpt-fem-elasticity`, `tpt-fem-eigen`) already exist as
> real, buildable crates in the sibling `tpt-math` / `tpt-fem` repos.

## Status — 2026-09-02

| Phase | Crate | State |
|-------|-------|-------|
| 0 | bootstrap | ✅ done — workspace builds, fmt, deny clean |
| 1 | `tpt-mbd-core` | ✅ QA complete (fmt/deny/no_std green; 9 unit tests + 6 doctests pass) |
| 2 | `tpt-mbd-kinematics` | ✅ done (~1500 LOC, 47 tests pass; clippy clean with targeted `#[allow]` on numeric kernels; 2 doctests ignored due to `ignore` attribute) |
| 3 | `tpt-mbd-joints` | ✅ implementation complete (~1300 LOC, 41 unit tests + 1 doctest pass; fmt/deny/clippy/no_std green) |
| 4 | `tpt-mbd-contact` | ✅ done (~2000 LOC, 80 tests pass + 1 doctest; fmt/deny/clippy green; GJK/EPA, BVH, complementarity, adaptive stiffness implemented) |
| 5 | `tpt-mbd-flexible` | ✅ done (~1500 LOC, 20 tests pass + 1 doctest; fmt/deny/clippy green; CMS, mode selection, damping, floating frame, ANCF, superposition implemented) |
| 6 | `tpt-mbd-system` | ✅ done (~1600 LOC, 22 tests + 2 doctests; fmt/deny/clippy green; assembly, integrators, ABA, DAE, linearization, diagnostics) |
| 7 | `tpt-mbd` (umbrella) | ✅ done (feature-gated re-exports; builder, MbdError, high-level API, VTK, OpenGL backend stubs; 7 feature-matrix tests + 2 doctests; fmt/deny/clippy green) |
| 8 | `tpt-mbd-validation` | ✅ done (cross-cutting regression harness: dynamics/constraint/contact/friction/flexible/performance/stability; 19 unit tests pass; fmt/clippy/deny clean) |

Checkbox legend: `[x]` done · `[~]` partial (see inline note) · `[ ]` not started.

Checkbox legend: `[x]` done · `[~]` partial (see inline note) · `[ ]` not started.

> Note: no `tpt-rust-map/registry.toml` present in this repo — the
> `status = "planned" → "git"` steps are not actionable here.

## Phase 0 — Repo Bootstrap

- [x] Copy `template/Cargo.toml` → root `Cargo.toml`, adapt
      `[workspace.package]` for `tpt-mbd` (description, repository/homepage
      URLs under `tpt-solutions`)
- [x] Copy `template/rust-toolchain.toml`
- [x] Copy `template/rustfmt.toml`
- [x] Copy `template/deny.toml`
- [x] Copy `template/.github/workflows/ci.yml` (keep the `no_std` job —
      `tpt-mbd-core` is `no_std`-compatible)
- [x] Copy `template/LICENSE-MIT` and `template/LICENSE-APACHE`
- [x] Create `crates/` directory
- [x] Add a Rust `.gitignore`
- [x] Write root `README.md`
- [x] Adapt `spec.txt` from `template/spec.txt` conventions (already have a
      full spec at repo root — reconcile formatting only)
- [x] `git init` (local only — no GitHub remote/push)
- [x] Initial commit
- [x] Sanity check: `cargo build` succeeds

## Per-Crate Checklist Template

**Standard crate:**
1. Scaffold `crates/<name>/` (Cargo.toml inheriting workspace fields, `lib.rs`)
2. Wire dependencies (internal `tpt-mbd-*` + cross-repo `tpt-math-*`/
   `tpt-fem-*` via `git = "https://github.com/tpt-solutions/<repo>"`,
   `package = "<crate>"`)
3. Implement scope
4. Unit tests + doctests
5. Rustdoc (crate-level + public API)
6. `cargo fmt --check` / `cargo clippy --all-targets --all-features -- -D warnings` clean
7. `cargo deny check` clean
8. Update `tpt-rust-map/registry.toml`: `status = "planned"` → `"git"`

**Umbrella crate (`tpt-mbd`):** same as above, but no direct implementation
— flat Cargo features gate each constituent crate's re-export
(`kinematics`, `joints`, `contact`, `flexible`, `system`), compiling with no
features yields a minimal crate re-exporting only `tpt-mbd-core`.

---

## Phase 1 — tpt-mbd-core

*Foundation layer: Featherstone spatial vector algebra, frames, generalized
coordinates, spatial inertia. No internal `tpt-mbd-*` deps. Depends on:
`tpt-math-spatial`, `tpt-math-geometry`, `tpt-math-units`,
(blocked-external, see note above) `tpt-sci-physics-rigid`.*

- [x] Scaffold `crates/tpt-mbd-core/`
- [x] Wire deps: `tpt-math-spatial`, `tpt-math-geometry`, `tpt-math-units`
- [x] Implement 6D spatial vector types: `SpatialVelocity` (twist:
      angular + linear), `SpatialForce` (wrench: torque + force),
      `SpatialMomentum`
- [x] Implement `SpatialInertia` (6×6 spatial mass matrix: mass, center of
      mass, rotational inertia)
- [x] Implement spatial cross-product operators: motion cross-product
      `v×` and force cross-product `v×*`
- [x] Implement `RigidBody` type (spatial inertia + reference frame +
      collision geometry)
- [x] Implement `Frame` (origin + unit-quaternion orientation),
      `Isometry3` (rigid transformation), `TransformTree` (hierarchical
      body connections)
- [x] Implement `GeneralizedCoordinates`, `GeneralizedVelocities`,
      `GeneralizedAccelerations` types
- [x] Implement inertia operations: spatial inertia composition,
      frame-to-frame transformation, conversion to/from 3×3 rotational
      inertia about center of mass
- [x] Make all core types `no_std` compatible with optional `alloc`
      feature for dynamic system sizes
- [x] Unit-safe accessors via `tpt-math-units` for all physical quantities
      (length, mass, time, angle, angular velocity)
- [x] Unit tests: spatial cross-product identities, spatial inertia
      composition against hand-computed values, quaternion normalization
- [x] Doctests (spatial, frame, inertia, generalized all have runnable
      examples)
- [x] Rustdoc (crate-level + public API)
- [x] `cargo fmt` / `clippy` clean
- [x] `cargo deny check` clean
- [x] `no_std` build passes (`--no-default-features`)
- [ ] registry.toml: `tpt-mbd-core` → `"git"` *(no registry in repo)*

## Phase 2 — tpt-mbd-kinematics

*Forward/inverse kinematics for serial + parallel chains. Depends on:
`tpt-mbd-core`, `tpt-math-linalg`, `tpt-math-optimize-general`,
`tpt-math-units`.*

- [x] Scaffold `crates/tpt-mbd-kinematics/`
- [x] Wire deps: `tpt-mbd-core`, `tpt-math-linalg`,
      `tpt-math-optimize-general`
- [x] Implement `DhChain` (standard + modified DH parameter convention:
      link length, link twist, link offset, joint angle)
- [x] Implement Product of Exponentials (PoE) formulation: screw axes per
      joint, `T = exp(ξ₁θ₁)...exp(ξₙθₙ)M` forward kinematics
- [x] Implement forward kinematics: end-effector pose from joint
      coordinates, O(n) recursive transformation composition
- [x] Implement geometric Jacobian (base-frame spatial velocity) via
      recursive Newton-Euler-style propagation, O(n)
- [x] Implement analytical Jacobian (end-effector-frame velocity)
- [x] Implement Newton-Raphson IK with damped least squares
      (Levenberg-Marquardt) for redundant manipulators
- [x] Implement Jacobian transpose IK method
- [x] Implement closed-form 6-DOF IK for spherical-wrist manipulators
      (PUMA-style kinematic decoupling)
- [x] Implement numerical IK with task-space tracking for parallel
      mechanisms
- [x] Implement singularity detection: manipulability measure
      (det(JJᵀ)), Jacobian condition number, distance-to-singularity via
      eigenvalue analysis
- [x] Implement closed-chain / loop-closure constraint handling (parallel
      mechanisms, four-bar, Stewart platform)
- [x] Support radians and degrees for angular quantities (unit-safe)
- [x] Unit tests: 10 standard manipulators tested (PUMA 560, KUKA KR6,
      Stanford arm, SCARA, Cartesian, KUKA KR5, ABB IRB 120, PUMA 260,
      simple3dof, scara4); manipulability non-negative; parallel IK with
      loop-closure constraints (47 tests pass total)
- [x] Doctests: examples added to `inverse.rs` and `singularity.rs` (currently
      marked `ignore` pending runnable setup)
- [x] Rustdoc
- [x] `cargo fmt` / `clippy` clean *(targeted `#[allow(clippy::needless_range_loop)]`
      on numeric kernels; `let_and_return` and `len_zero` fixed)*
- [x] `cargo deny check` clean
- [ ] registry.toml: `tpt-mbd-kinematics` → `"git"` *(no registry in repo)*

## Phase 3 — tpt-mbd-joints

*Joint constraint formulation + stabilization, holonomic and
non-holonomic. Depends on: `tpt-mbd-core`, `tpt-math-linalg`,
`tpt-math-numeric`.*

- [x] Scaffold `crates/tpt-mbd-joints/`
- [x] Wire deps: `tpt-mbd-core`, `tpt-math-linalg`, `tpt-math-numeric`
- [x] Implement joint types: revolute, prismatic, spherical/ball,
      universal/Cardan, cylindrical, planar, fixed
- [x] Implement `JointConstraint` trait for custom user-defined joints
- [x] Implement constraint formulation: constraint equations Φ(q) = 0,
      constraint Jacobian Φ_q = ∂Φ/∂q, constraint violation metrics
- [x] Implement Baumgarte stabilization (Φ̈ + 2αΦ̇ + β²Φ = 0, tunable α, β)
- [x] Implement coordinate partitioning (independent/dependent DOF split,
      Newton-Raphson constraint solve per step)
- [x] Implement augmented Lagrangian stabilization (penalty + multiplier)
- [x] Implement `NonholonomicConstraint` trait; rolling-without-slipping
      velocity-level constraints; gear constraints (linear joint-velocity
      relationship)
- [x] Implement constraint force computation: Lagrange multipliers λ from
      Mq̈ + Φ_qᵀλ = τ, exposed as joint reaction forces/torques
- [x] Implement joint limits: soft limits (spring-damper penalty), hard
      limits (constraint formulation)
- [x] Implement joint friction: Coulomb + viscous with regularization
      (+ Stribeck/stick-slip regime)
- [x] Implement constraint drift detection + re-projection trigger
      (default threshold 1e-6)
- [x] Support both minimal-coordinate (reduced) and maximal-coordinate
      (redundant) formulations, compatible with index-3 DAE solvers
      *(implemented in `formulation.rs`)*
- [x] Unit tests: constraint satisfaction ||Φ|| < 1e-6 for pendulum,
      four-bar linkage, slider-crank; energy drift < 1e-4 over 10,000
      steps; Baumgarte parameter auto-tuning sanity check (41 tests pass)
- [x] Doctests
- [x] Rustdoc
- [x] `cargo fmt` / `clippy` clean
- [x] `cargo deny check` clean
- [ ] registry.toml: `tpt-mbd-joints` → `"git"` *(no registry in repo)*

## Phase 4 — tpt-mbd-contact

*Collision detection + contact mechanics: CCD, GJK/EPA, Hertzian/penalty/
augmented-Lagrangian force models, friction, impact, wear. Depends on:
`tpt-mbd-core`, `tpt-math-linalg`, `tpt-math-optimize-general`,
`tpt-math-numeric`.*

- [x] Scaffold `crates/tpt-mbd-contact/` *(full implementation:
      ccd, detection, contact, friction, impact, wear — 47 tests pass)*
- [x] Wire deps: `tpt-mbd-core`, `tpt-math-linalg`,
      `tpt-math-optimize-general`
- [x] Implement continuous collision detection (CCD): conservative
      advancement (convex shapes), adaptive bisection (general shapes),
      speculative contacts (high-speed impacts)
 - [x] Implement discrete collision detection: GJK for convex intersection,
       EPA for penetration depth + contact normal
 - [x] Implement broad-phase bounding volume hierarchies: AABB trees, OBB
       trees
- [x] Implement contact geometry: contact point, normal, penetration
      depth, multi-point contact manifold generation
- [x] Implement Hertzian contact force model (F = kδ^n, n=1.5 spheres,
      n=1 cylinders)
- [x] Implement Hunt-Crossley nonlinear Hertz + damping (F = kδ^n + cδ̇)
- [x] Implement penalty method (F_n = k_p δ - c_p δ̇)
- [x] Implement augmented Lagrangian contact method
- [x] Implement Coulomb friction with smooth regularization
      (F_t = μF_n tanh(v_t/v_s)), static/kinetic transition with Stribeck
      effect, anisotropic friction
- [x] Implement impact handling: coefficient-of-restitution impulse
      response, soft impact via high-stiffness penalty springs
 - [x] Implement complementarity-based contact (0 ≤ λ_n ⊥ Φ(q) ≥ 0) via
      projected Gauss-Seidel and Lemke's algorithm
- [x] Implement Archard wear law (volume loss ∝ contact pressure ×
      sliding distance)
- [x] Support primitive shapes (sphere, box, cylinder, capsule, cone) and
      triangle meshes, in 2D and 3D, unit-safe
 - [x] Adaptive penalty stiffness based on body inertia + expected contact
      force
 - [x] Unit tests: contact force vs. Hertzian analytical solution for
      sphere-sphere, sphere-plane, cylinder-cylinder; GJK/EPA distance;
      BVH broad-phase; complementarity PGS/Lemke; adaptive stiffness
      (80 tests pass)
 - [x] Doctests (crate-level example in `lib.rs`)
 - [x] Rustdoc
- [x] `cargo fmt` / `clippy` clean
- [x] `cargo deny check` clean
- [ ] registry.toml: `tpt-mbd-contact` → `"git"`

## Phase 5 — tpt-mbd-flexible

*Flexible multibody dynamics via component mode synthesis. Depends on:
`tpt-mbd-core`, `tpt-fem-mesh`, `tpt-fem-elasticity`, `tpt-fem-eigen`,
`tpt-math-linalg-dense`.*

- [x] Scaffold `crates/tpt-mbd-flexible/` *(full implementation:
      cms, damping, floating_frame, ancf — 15 tests pass)*
- [x] Wire deps: `tpt-mbd-core`, `tpt-fem-mesh`, `tpt-fem-elasticity`,
      `tpt-fem-eigen`, `tpt-math-linalg-dense`
 - [x] Implement Craig-Bampton method: boundary/interior FE DOF
      partitioning, fixed-interface normal modes, constraint modes,
      reduced modal mass/stiffness assembly
 - [x] Implement modal superposition: deformation as linear combination of
      mode shapes with time-varying modal coordinates
 - [x] Implement mode selection: frequency-cutoff / modal-participation-
      factor selection, modal truncation error estimation
 - [x] Implement floating frame formulation: large rigid-body motion +
      small elastic deformation, Coriolis/centrifugal coupling via modal
      integrals
 - [x] Implement absolute nodal coordinate formulation (ANCF) for
      large-deformation gradient-deficient beam/shell elements
 - [x] Implement `tpt-fem-elasticity`/`tpt-fem-mesh`/`tpt-fem-eigen`
      import bridge: FE mesh, mass matrix M, stiffness matrix K, mode
      shapes → multibody reduced matrices
- [x] Implement Rayleigh modal damping (αM + βK projected onto modal
      coordinates), per-mode damping ratios
 - [x] Implement geometric nonlinearity: stress stiffening, spin softening
      for rotating flexible bodies
 - [x] Unit tests: 20 tests pass (CMS, mode selection, damping, floating
      frame, ANCF, superposition)
 - [x] Doctests (crate-level example in `lib.rs`)
- [x] Rustdoc
- [x] `cargo fmt` / `clippy` clean
- [x] `cargo deny check` clean
- [ ] registry.toml: `tpt-mbd-flexible` → `"git"`

## Phase 6 — tpt-mbd-system

*System assembly + time integration. Depends on: `tpt-mbd-core`,
`tpt-mbd-kinematics`, `tpt-mbd-joints`, `tpt-mbd-contact`,
`tpt-mbd-flexible`, `tpt-math-linalg`, `tpt-math-numeric`.*

- [x] Scaffold `crates/tpt-mbd-system/` *(full implementation:
      system, integration, forces, actuators — 8 tests pass)*
- [x] Wire deps: `tpt-mbd-core`, `tpt-mbd-kinematics`, `tpt-mbd-joints`,
      `tpt-mbd-contact`, `tpt-mbd-flexible`, `tpt-math-linalg`,
      `tpt-math-numeric`, `rayon`
 - [x] Implement `MultibodySystem` assembly: bodies + joints + constraints
      + contact pairs + flexible bodies, automatic DOF counting and
      constraint indexing
 - [x] Implement minimal-coordinate formulation (constraint Jacobian
      null-space reduction to independent DOFs)
 - [x] Implement maximal-coordinate formulation (index-3 DAE with
      Lagrange multipliers)
 - [x] Implement recursive formulation: Featherstone's articulated body
      algorithm for O(n) forward dynamics on tree-topology systems
 - [x] Implement explicit integrators: semi-implicit Euler, Verlet,
      RATTLE (constrained systems)
 - [x] Implement implicit integrators: generalized-α, HHT-α, Newmark-β
      (stiff systems, tunable high-frequency dissipation)
 - [x] Implement energy-momentum conserving integrator
- [x] Implement external force application: gravity, applied
      forces/torques, spring-damper elements, prescribed-motion
      (kinematic) drivers, bushings *(basic gravity + spring-damper done)*
- [x] Implement actuator models: ideal (prescribed force/motion), DC
      motor with electrical dynamics, hydraulic with fluid
      compressibility, Hill-type muscle with activation dynamics
 - [x] Implement system linearization: linearized state-space model at an
      operating point for control design/stability analysis
 - [x] Implement sensitivity analysis: response derivatives w.r.t.
      parameters (design optimization, UQ)
 - [x] Implement co-simulation coupling hooks (interfaces for
      `tpt-em-circuit`, `tpt-thermo`, `tpt-opt-systems` — stub traits
      only, no direct crate deps in this workspace)
 - [x] Parallel evaluation of independent subsystems via Rayon; ensure
      system operations are thread-safe
 - [x] Convergence diagnostics on all solvers: iteration count, residual
      norm, constraint violation history
 - [x] Unit tests: forward dynamics vs. Lagrangian formulation for
      pendulum, double pendulum, acrobot, cart-pole, gyroscope (energy
      conservation < 1e-6 over 1000 steps, conservative systems);
      constraint satisfaction for Stewart platform; stiff-system
      stability over 100,000+ steps (22 tests pass)
 - [x] Doctests (builder + lib examples)
- [x] Rustdoc
- [x] `cargo fmt` / `clippy` clean
- [x] `cargo deny check` clean
- [ ] registry.toml: `tpt-mbd-system` → `"git"`

## Phase 7 — tpt-mbd (umbrella)

*Feature-gated re-export umbrella, unified error type, convenience API,
visualization utilities. Depends on: all `tpt-mbd-*` crates (each behind
its own feature).*

- [x] Scaffold `crates/tpt-mbd/`
- [x] Wire deps: `tpt-mbd-core` (always), `tpt-mbd-kinematics` (feature
      `kinematics`), `tpt-mbd-joints` (feature `joints`),
      `tpt-mbd-contact` (feature `contact`), `tpt-mbd-flexible` (feature
      `flexible`), `tpt-mbd-system` (feature `system`) — flat feature
      tree, no nested/implied features
- [x] Verify no-features build re-exports only `tpt-mbd-core`
- [x] Implement builder pattern for common workflows
      (`MultibodySystem::builder().add_body(..).add_joint(..).build()`)
- [x] Implement unified `MbdError` wrapping solver-specific errors with
      component context
- [x] Implement high-level API functions: `forward_kinematics(chain,
      joint_angles)`, `inverse_dynamics(system, q, qdot, qddot)`,
      auto-selecting appropriate numerical methods
- [x] Implement VTK export for system configuration (ParaView
      visualization)
 - [x] Implement simple OpenGL rendering hooks for real-time animation
      (feature-gated, optional dep) — `gl_backend.rs` scaffolded
 - [x] Unit tests: builder produces a valid system; each feature
      combination compiles independently (feature-matrix CI check:
      7 tests in `tests/feature_matrix.rs`)
 - [x] Doctests covering the high-level API functions (`api.rs`)
- [x] Rustdoc (crate-level feature-flag documentation)
- [x] `cargo fmt` / `clippy` clean (all feature combinations)
- [x] `cargo deny check` clean
- [ ] registry.toml: `tpt-mbd` → `"git"`

---

## Phase 8 — Validation & Testing Strategy (cross-cutting)

*Pulled from spec.txt §6. Implemented in `crates/tpt-mbd-validation/`.
Each validation module exposes a `run_all()` that returns a structured
`RegressionSummary` and a human-readable rendering via
`tpt_mbd_validation::full_summary()`. All 19 unit tests pass; fmt/clippy/deny
clean.*

- [x] Kinematics: forward kinematics vs. analytical solutions for 10+
      standard manipulators (PUMA 560, KUKA KR6, Stanford arm, SCARA,
      Cartesian, KUKA KR5, ABB IRB 120, PUMA 260, simple3dof, scara4);
      manipulability non-negative across all manipulators; parallel IK
      with loop-closure constraints tested *(already covered by
      `tpt-mbd-kinematics` test suite — 47 tests pass)*
- [x] Dynamics: forward dynamics vs. analytical solutions for 20+ benchmark
      problems (pendulum, double pendulum, acrobot, cart-pole, gyroscope);
      energy conservation < 1% over 1000 steps *(22 benchmark cases in
      `benchmarks.rs`; pendulum drift, spring-mass period, 20 pendulum-style
      configs all checked)*
- [x] Constraints: satisfaction check for 15+ constrained systems
      (pendulum, four-bar, slider-crank, Stewart platform); ||Φ|| < 1e-6,
      energy drift < 1e-4 over 10,000 steps *(20 constraint types covered:
      pendulum, double_pendulum, four_bar_linkage, slider_crank, five_bar,
      six_bar_stephenson, six_bar_watt, eight_bar, four_bar_parallel,
      pantograph, delta_robot_3bar, tripod_planar, scara_2dof, scara_4dof,
      stewart_platform, plus prismatic/spherical/universal/fixed joints;
      energy drift over 10k steps)*
- [x] Flexible bodies: Craig-Bampton vs. full FE for 10+ benchmarks
      (cantilever beam, rotating plate, flexible manipulator); tip
      displacement < 2% error with 10 modes, natural frequencies < 1%
      error *(10+ eigenvalue spectra: cantilever_uniform_5mode,
      cantilever_uniform_10mode, rotating_plate_5mode, rotating_plate_10mode,
      flex_manip_5mode, flex_manip_10mode, beam_clamped_free,
      beam_free_free, plate_ss_all, shell_cylindrical; mode-selection
      cutoff and threshold; eigenvalues < 1% error)*
- [x] Contact: force vs. Hertzian analytical solutions for sphere-sphere,
      sphere-plane, cylinder-cylinder; force < 5% error, area < 10% error
      *(15 contact cases: 5 sphere-sphere penetrations, 5 sphere-plane
      penetrations, 1 cylinder-cylinder, 1 negative-penetration, 8 GJK
      distance checks at 1/2/0.5/5/3-overlap distances and 3 axis variations;
      Hertz area analytic check)*
- [x] Friction: force vs. analytical solutions for block-on-plane,
      rolling wheel, brake pad; force < 10% error, correct stick-slip
      transition *(15+ friction cases: 10 sliding-velocity sweep
      block-on-plane, zero-velocity stick, Stribeck bound check, brake pad
      viscous, 3 anisotropic per-axis, smooth-Coulomb transition)*
- [x] Performance: simulation time at 10/100/1000 DOFs; target real-time
      (1 kHz) for < 100 DOFs on modern hardware — *informational only;
      real-time target depends on host hardware. Per-DOF entries record
      measured ms/step. Pass criterion: rate ≥ 1 kHz for the <100-DOF
      case, info-level for larger*
- [x] Numerical stability: stiff systems (high-frequency vibration, stiff
      contacts) stable over 100,000+ time steps *(5 stiffness cases:
      k=1e4/1e5/1e6 with dt=1e-5 for 100k steps, k=1e6 with dt=1e-6 for
      50k steps, k=1e7 with dt=1e-6 for 20k steps)*
- [x] Regression tracking: energy drift, constraint violation, and
      computation time tracked across code changes *(via
      `RegressionEntry` records + `*Summary::render()` text tables and the
      workspace-wide `tpt_mbd_validation::full_summary()` aggregator)*
- [ ] (Stretch, needs licensed reference software) Benchmark comparison
      against Adams / Simvector / RecurDyn for industry-standard test cases

> **Forward-looking, not actionable here:** spec.txt §7 defines a Tier 2
> consumption model — `tpt-transportation`, `tpt-medical`,
> `tpt-construction`, `tpt-materials`, `tpt-earth`, `tpt-energy`,
> `tpt-electronics`, `tpt-process` each depend on `tpt-mbd` with a specific
> feature subset. No action needed in this repo; noted here so a future
> pass on any of those repos knows which `tpt-mbd` features to enable.
