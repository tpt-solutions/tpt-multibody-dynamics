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

## Phase 0 — Repo Bootstrap

- [ ] Copy `template/Cargo.toml` → root `Cargo.toml`, adapt
      `[workspace.package]` for `tpt-mbd` (description, repository/homepage
      URLs under `tpt-solutions`)
- [ ] Copy `template/rust-toolchain.toml`
- [ ] Copy `template/rustfmt.toml`
- [ ] Copy `template/deny.toml`
- [ ] Copy `template/.github/workflows/ci.yml` (keep the `no_std` job —
      `tpt-mbd-core` is `no_std`-compatible)
- [ ] Copy `template/LICENSE-MIT` and `template/LICENSE-APACHE`
- [ ] Create `crates/` directory
- [ ] Add a Rust `.gitignore`
- [ ] Write root `README.md`
- [ ] Adapt `spec.txt` from `template/spec.txt` conventions (already have a
      full spec at repo root — reconcile formatting only)
- [ ] `git init` (local only — no GitHub remote/push)
- [ ] Initial commit
- [ ] Sanity check: `cargo build` succeeds

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

- [ ] Scaffold `crates/tpt-mbd-core/`
- [ ] Wire deps: `tpt-math-spatial`, `tpt-math-geometry`, `tpt-math-units`
- [ ] Implement 6D spatial vector types: `SpatialVelocity` (twist:
      angular + linear), `SpatialForce` (wrench: torque + force),
      `SpatialMomentum`
- [ ] Implement `SpatialInertia` (6×6 spatial mass matrix: mass, center of
      mass, rotational inertia)
- [ ] Implement spatial cross-product operators: motion cross-product
      `v×` and force cross-product `v×*`
- [ ] Implement `RigidBody` type (spatial inertia + reference frame +
      collision geometry)
- [ ] Implement `Frame` (origin + unit-quaternion orientation),
      `Isometry3` (rigid transformation), `TransformTree` (hierarchical
      body connections)
- [ ] Implement `GeneralizedCoordinates`, `GeneralizedVelocities`,
      `GeneralizedAccelerations` types
- [ ] Implement inertia operations: spatial inertia composition,
      frame-to-frame transformation, conversion to/from 3×3 rotational
      inertia about center of mass
- [ ] Make all core types `no_std` compatible with optional `alloc`
      feature for dynamic system sizes
- [ ] Unit-safe accessors via `tpt-math-units` for all physical quantities
      (length, mass, time, angle, angular velocity)
- [ ] Unit tests: spatial cross-product identities, spatial inertia
      composition against hand-computed values, quaternion normalization
- [ ] Doctests
- [ ] Rustdoc (crate-level + public API)
- [ ] `cargo fmt` / `clippy` clean
- [ ] `cargo deny check` clean
- [ ] `no_std` CI job passes (default-`alloc`-off build)
- [ ] registry.toml: `tpt-mbd-core` → `"git"`

## Phase 2 — tpt-mbd-kinematics

*Forward/inverse kinematics for serial + parallel chains. Depends on:
`tpt-mbd-core`, `tpt-math-linalg`, `tpt-math-optimize-general`,
`tpt-math-units`.*

- [ ] Scaffold `crates/tpt-mbd-kinematics/`
- [ ] Wire deps: `tpt-mbd-core`, `tpt-math-linalg`,
      `tpt-math-optimize-general`
- [ ] Implement `DhChain` (standard + modified DH parameter convention:
      link length, link twist, link offset, joint angle)
- [ ] Implement Product of Exponentials (PoE) formulation: screw axes per
      joint, `T = exp(ξ₁θ₁)...exp(ξₙθₙ)M` forward kinematics
- [ ] Implement forward kinematics: end-effector pose from joint
      coordinates, O(n) recursive transformation composition
- [ ] Implement geometric Jacobian (base-frame spatial velocity) via
      recursive Newton-Euler-style propagation, O(n)
- [ ] Implement analytical Jacobian (end-effector-frame velocity)
- [ ] Implement Newton-Raphson IK with damped least squares
      (Levenberg-Marquardt) for redundant manipulators
- [ ] Implement Jacobian transpose IK method
- [ ] Implement closed-form 6-DOF IK for spherical-wrist manipulators
      (PUMA-style kinematic decoupling)
- [ ] Implement numerical IK with task-space tracking for parallel
      mechanisms
- [ ] Implement singularity detection: manipulability measure
      (det(JJᵀ)), Jacobian condition number, distance-to-singularity via
      eigenvalue analysis
- [ ] Implement workspace analysis: reachable-workspace boundary tracing,
      dexterous workspace identification
- [ ] Implement closed-chain / loop-closure constraint handling (parallel
      mechanisms, four-bar, Stewart platform)
- [ ] Support radians and degrees for angular quantities (unit-safe)
- [ ] Unit tests: forward kinematics vs. analytical solutions for PUMA
      560, KUKA KR6, Stanford arm (position < 1e-9 m, orientation
      < 1e-9 rad); IK convergence < 1e-6 m against closed-form targets;
      manipulability at known singular configurations
- [ ] Doctests
- [ ] Rustdoc
- [ ] `cargo fmt` / `clippy` clean
- [ ] `cargo deny check` clean
- [ ] registry.toml: `tpt-mbd-kinematics` → `"git"`

## Phase 3 — tpt-mbd-joints

*Joint constraint formulation + stabilization, holonomic and
non-holonomic. Depends on: `tpt-mbd-core`, `tpt-math-linalg`,
`tpt-math-numeric`.*

- [ ] Scaffold `crates/tpt-mbd-joints/`
- [ ] Wire deps: `tpt-mbd-core`, `tpt-math-linalg`, `tpt-math-numeric`
- [ ] Implement joint types: revolute, prismatic, spherical/ball,
      universal/Cardan, cylindrical, planar, fixed
- [ ] Implement `JointConstraint` trait for custom user-defined joints
- [ ] Implement constraint formulation: constraint equations Φ(q) = 0,
      constraint Jacobian Φ_q = ∂Φ/∂q, constraint violation metrics
- [ ] Implement Baumgarte stabilization (Φ̈ + 2αΦ̇ + β²Φ = 0, tunable α, β)
- [ ] Implement coordinate partitioning (independent/dependent DOF split,
      Newton-Raphson constraint solve per step)
- [ ] Implement augmented Lagrangian stabilization (penalty + multiplier)
- [ ] Implement `NonholonomicConstraint` trait; rolling-without-slipping
      velocity-level constraints; gear constraints (linear joint-velocity
      relationship)
- [ ] Implement constraint force computation: Lagrange multipliers λ from
      Mq̈ + Φ_qᵀλ = τ, exposed as joint reaction forces/torques
- [ ] Implement joint limits: soft limits (spring-damper penalty), hard
      limits (constraint formulation)
- [ ] Implement joint friction: Coulomb + viscous with regularization
- [ ] Implement constraint drift detection + re-projection trigger
      (default threshold 1e-6)
- [ ] Support both minimal-coordinate (reduced) and maximal-coordinate
      (redundant) formulations, compatible with index-3 DAE solvers
- [ ] Unit tests: constraint satisfaction ||Φ|| < 1e-6 for pendulum,
      four-bar linkage, slider-crank; energy drift < 1e-4 over 10,000
      steps; Baumgarte parameter auto-tuning sanity check
- [ ] Doctests
- [ ] Rustdoc
- [ ] `cargo fmt` / `clippy` clean
- [ ] `cargo deny check` clean
- [ ] registry.toml: `tpt-mbd-joints` → `"git"`

## Phase 4 — tpt-mbd-contact

*Collision detection + contact mechanics: CCD, GJK/EPA, Hertzian/penalty/
augmented-Lagrangian force models, friction, impact, wear. Depends on:
`tpt-mbd-core`, `tpt-math-linalg`, `tpt-math-optimize-general`,
`tpt-math-numeric`.*

- [ ] Scaffold `crates/tpt-mbd-contact/`
- [ ] Wire deps: `tpt-mbd-core`, `tpt-math-linalg`,
      `tpt-math-optimize-general`
- [ ] Implement continuous collision detection (CCD): conservative
      advancement (convex shapes), adaptive bisection (general shapes),
      speculative contacts (high-speed impacts)
- [ ] Implement discrete collision detection: GJK for convex intersection,
      EPA for penetration depth + contact normal
- [ ] Implement broad-phase bounding volume hierarchies: AABB trees, OBB
      trees
- [ ] Implement contact geometry: contact point, normal, penetration
      depth, multi-point contact manifold generation
- [ ] Implement Hertzian contact force model (F = kδ^n, n=1.5 spheres,
      n=1 cylinders)
- [ ] Implement Hunt-Crossley nonlinear Hertz + damping (F = kδ^n + cδ̇)
- [ ] Implement penalty method (F_n = k_p δ - c_p δ̇)
- [ ] Implement augmented Lagrangian contact method
- [ ] Implement Coulomb friction with smooth regularization
      (F_t = μF_n tanh(v_t/v_s)), static/kinetic transition with Stribeck
      effect, anisotropic friction
- [ ] Implement impact handling: coefficient-of-restitution impulse
      response, soft impact via high-stiffness penalty springs
- [ ] Implement complementarity-based contact (0 ≤ λ_n ⊥ Φ(q) ≥ 0) via
      projected Gauss-Seidel and Lemke's algorithm
- [ ] Implement Archard wear law (volume loss ∝ contact pressure ×
      sliding distance)
- [ ] Support primitive shapes (sphere, box, cylinder, capsule, cone) and
      triangle meshes, in 2D and 3D, unit-safe
- [ ] Adaptive penalty stiffness based on body inertia + expected contact
      force
- [ ] Unit tests: contact force vs. Hertzian analytical solution for
      sphere-sphere, sphere-plane, cylinder-cylinder (force < 5% error,
      area < 10% error); friction force vs. analytical solutions for
      block-on-plane, rolling wheel, brake pad (< 10% error, correct
      stick-slip transition)
- [ ] Doctests
- [ ] Rustdoc
- [ ] `cargo fmt` / `clippy` clean
- [ ] `cargo deny check` clean
- [ ] registry.toml: `tpt-mbd-contact` → `"git"`

## Phase 5 — tpt-mbd-flexible

*Flexible multibody dynamics via component mode synthesis. Depends on:
`tpt-mbd-core`, `tpt-fem-mesh`, `tpt-fem-elasticity`, `tpt-fem-eigen`,
`tpt-math-linalg-dense`.*

- [ ] Scaffold `crates/tpt-mbd-flexible/`
- [ ] Wire deps: `tpt-mbd-core`, `tpt-fem-mesh`, `tpt-fem-elasticity`,
      `tpt-fem-eigen`, `tpt-math-linalg-dense`
- [ ] Implement Craig-Bampton method: boundary/interior FE DOF
      partitioning, fixed-interface normal modes (eigenvectors of K_ii
      with boundary fixed), constraint modes (static shapes from unit
      boundary displacement), reduced modal mass/stiffness assembly
- [ ] Implement modal superposition: deformation as linear combination of
      mode shapes with time-varying modal coordinates
- [ ] Implement mode selection: frequency-cutoff / modal-participation-
      factor selection, modal truncation error estimation
- [ ] Implement floating frame formulation: large rigid-body motion +
      small elastic deformation, Coriolis/centrifugal coupling via modal
      integrals (consistent mass matrix projection)
- [ ] Implement absolute nodal coordinate formulation (ANCF) for
      large-deformation gradient-deficient beam/shell elements
- [ ] Implement `tpt-fem-elasticity`/`tpt-fem-mesh`/`tpt-fem-eigen`
      import bridge: FE mesh, mass matrix M, stiffness matrix K, mode
      shapes → multibody reduced matrices
- [ ] Implement Rayleigh modal damping (αM + βK projected onto modal
      coordinates), per-mode damping ratios
- [ ] Implement geometric nonlinearity: stress stiffening (stress-
      dependent stiffness updates), spin softening (centrifugal effects
      on stiffness) for rotating flexible bodies
- [ ] Unit tests: Craig-Bampton vs. full FE simulation for cantilever
      beam, rotating plate, flexible manipulator (tip displacement < 2%
      error with 10 modes, natural frequencies < 1% error)
- [ ] Doctests
- [ ] Rustdoc
- [ ] `cargo fmt` / `clippy` clean
- [ ] `cargo deny check` clean
- [ ] registry.toml: `tpt-mbd-flexible` → `"git"`

## Phase 6 — tpt-mbd-system

*System assembly + time integration. Depends on: `tpt-mbd-core`,
`tpt-mbd-kinematics`, `tpt-mbd-joints`, `tpt-mbd-contact`,
`tpt-mbd-flexible`, `tpt-math-linalg`, `tpt-math-numeric`.*

- [ ] Scaffold `crates/tpt-mbd-system/`
- [ ] Wire deps: `tpt-mbd-core`, `tpt-mbd-kinematics`, `tpt-mbd-joints`,
      `tpt-mbd-contact`, `tpt-mbd-flexible`, `tpt-math-linalg`,
      `tpt-math-numeric`, `rayon`
- [ ] Implement `MultibodySystem` assembly: bodies + joints + constraints
      + contact pairs + flexible bodies, automatic DOF counting and
      constraint indexing
- [ ] Implement minimal-coordinate formulation (constraint Jacobian
      null-space reduction to independent DOFs)
- [ ] Implement maximal-coordinate formulation (index-3 DAE with
      Lagrange multipliers)
- [ ] Implement recursive formulation: Featherstone's articulated body
      algorithm for O(n) forward dynamics on tree-topology systems
- [ ] Implement explicit integrators: semi-implicit Euler, Verlet,
      RATTLE (constrained systems)
- [ ] Implement implicit integrators: generalized-α, HHT-α, Newmark-β
      (stiff systems, tunable high-frequency dissipation)
- [ ] Implement energy-momentum conserving integrator
- [ ] Implement external force application: gravity, applied
      forces/torques, spring-damper elements, prescribed-motion
      (kinematic) drivers, bushings
- [ ] Implement actuator models: ideal (prescribed force/motion), DC
      motor with electrical dynamics, hydraulic with fluid
      compressibility, Hill-type muscle with activation dynamics
- [ ] Implement system linearization: linearized state-space model at an
      operating point for control design/stability analysis
- [ ] Implement sensitivity analysis: response derivatives w.r.t.
      parameters (design optimization, UQ)
- [ ] Implement co-simulation coupling hooks (interfaces for
      `tpt-em-circuit`, `tpt-thermo`, `tpt-opt-systems` — stub traits
      only, no direct crate deps in this workspace)
- [ ] Parallel evaluation of independent subsystems via Rayon; ensure
      system operations are thread-safe
- [ ] Convergence diagnostics on all solvers: iteration count, residual
      norm, constraint violation history
- [ ] Unit tests: forward dynamics vs. Lagrangian formulation for
      pendulum, double pendulum, acrobot, cart-pole, gyroscope (energy
      conservation < 1e-6 over 1000 steps, conservative systems);
      constraint satisfaction for Stewart platform; stiff-system
      stability over 100,000+ steps
- [ ] Doctests
- [ ] Rustdoc
- [ ] `cargo fmt` / `clippy` clean
- [ ] `cargo deny check` clean
- [ ] registry.toml: `tpt-mbd-system` → `"git"`

## Phase 7 — tpt-mbd (umbrella)

*Feature-gated re-export umbrella, unified error type, convenience API,
visualization utilities. Depends on: all `tpt-mbd-*` crates (each behind
its own feature).*

- [ ] Scaffold `crates/tpt-mbd/`
- [ ] Wire deps: `tpt-mbd-core` (always), `tpt-mbd-kinematics` (feature
      `kinematics`), `tpt-mbd-joints` (feature `joints`),
      `tpt-mbd-contact` (feature `contact`), `tpt-mbd-flexible` (feature
      `flexible`), `tpt-mbd-system` (feature `system`) — flat feature
      tree, no nested/implied features
- [ ] Verify no-features build re-exports only `tpt-mbd-core`
- [ ] Implement builder pattern for common workflows
      (`MultibodySystem::builder().add_body(..).add_joint(..).build()`)
- [ ] Implement unified `MbdError` wrapping solver-specific errors with
      component context (`KinematicsError::SingularConfiguration`,
      `ContactError::PenetrationTooLarge`, etc., per spec §4 Singularity
      Handling error enum shape)
- [ ] Implement high-level API functions: `forward_kinematics(chain,
      joint_angles)`, `inverse_dynamics(system, q, qdot, qddot)`,
      auto-selecting appropriate numerical methods
- [ ] Implement VTK export for system configuration (ParaView
      visualization)
- [ ] Implement simple OpenGL rendering hooks for real-time animation
      (feature-gated, optional dep)
- [ ] Unit tests: builder produces a valid system; each feature
      combination compiles independently (feature-matrix CI check)
- [ ] Doctests covering the high-level API functions
- [ ] Rustdoc (crate-level feature-flag documentation)
- [ ] `cargo fmt` / `clippy` clean (all feature combinations)
- [ ] `cargo deny check` clean
- [ ] registry.toml: `tpt-mbd` → `"git"`

---

## Phase 8 — Validation & Testing Strategy (cross-cutting)

*Pulled from spec.txt §6. These are workspace-level validation passes
that exercise multiple crates together, run after Phases 1-7 land.*

- [ ] Kinematics: forward kinematics vs. analytical solutions for 10+
      standard manipulators (PUMA 560, KUKA, Stanford); inverse
      kinematics vs. closed-form where available (position < 1e-9 m,
      orientation < 1e-9 rad FK; IK convergence < 1e-6 m)
- [ ] Dynamics: forward dynamics vs. Lagrangian formulation for 20+
      benchmark problems (pendulum, double pendulum, acrobot, cart-pole,
      gyroscope); energy conservation < 1e-6 over 1000 steps
- [ ] Constraints: satisfaction check for 15+ constrained systems
      (pendulum, four-bar, slider-crank, Stewart platform); ||Φ|| < 1e-6,
      energy drift < 1e-4 over 10,000 steps
- [ ] Flexible bodies: Craig-Bampton vs. full FE for 10+ benchmarks
      (cantilever beam, rotating plate, flexible manipulator); tip
      displacement < 2% error with 10 modes, natural frequencies < 1%
      error
- [ ] Contact: force vs. Hertzian analytical solutions for sphere-sphere,
      sphere-plane, cylinder-cylinder; force < 5% error, area < 10% error
- [ ] Friction: force vs. analytical solutions for block-on-plane,
      rolling wheel, brake pad; force < 10% error, correct stick-slip
      transition
- [ ] Performance: simulation time at 10/100/1000 DOFs; target real-time
      (1 kHz) for < 100 DOFs on modern hardware
- [ ] Numerical stability: stiff systems (high-frequency vibration, stiff
      contacts) stable over 100,000+ time steps
- [ ] Regression tracking: energy drift, constraint violation, and
      computation time tracked across code changes
- [ ] (Stretch, needs licensed reference software) Benchmark comparison
      against Adams / Simpack / RecurDyn for industry-standard test cases

> **Forward-looking, not actionable here:** spec.txt §7 defines a Tier 2
> consumption model — `tpt-transportation`, `tpt-medical`,
> `tpt-construction`, `tpt-materials`, `tpt-earth`, `tpt-energy`,
> `tpt-electronics`, `tpt-process` each depend on `tpt-mbd` with a specific
> feature subset. No action needed in this repo; noted here so a future
> pass on any of those repos knows which `tpt-mbd` features to enable.
