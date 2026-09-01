# tpt-multibody-dynamics

[![CI](https://github.com/tpt-solutions/tpt-multibody-dynamics/actions/workflows/ci.yml/badge.svg)](https://github.com/tpt-solutions/tpt-multibody-dynamics/actions/workflows/ci.yml)
[![license: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](LICENSE-MIT)

Engineering-grade multibody dynamics (MBD) for robotics, vehicle dynamics, biomechanics, and mechanism simulation. Provides forward/inverse kinematics, holonomic/non-holonomic joint constraints, flexible body reduction via Craig-Bampton CMS, and contact mechanics (Hertzian contact, Coulomb friction, impact). All implementations are from scratch in pure Rust with no external multibody library dependencies, ensuring full numerical control, deterministic reproducibility, and WASM compatibility.

## Crates

| Crate | Purpose |
|-------|---------|
| `tpt-mbd-core` | Spatial vector algebra, frames, inertia, generalized coordinates (`no_std`-compatible) |
| `tpt-mbd-kinematics` | DH parameters, PoE, forward/inverse kinematics, Jacobians, singularity detection |
| `tpt-mbd-joints` | Joint types, constraint formulation, Baumgarte stabilization, friction |
| `tpt-mbd-contact` | Collision detection (CCD, GJK/EPA), Hertzian contact, friction, impact |
| `tpt-mbd-flexible` | Craig-Bampton CMS, modal superposition, ANCF, Rayleigh damping |
| `tpt-mbd-system` | System assembly, recursive dynamics, time integration, actuators |
| `tpt-mbd` | Feature-gated umbrella crate with builder pattern and unified error type |

## Build Order

1. `tpt-mbd-core`
2. `tpt-mbd-kinematics`
3. `tpt-mbd-joints`
4. `tpt-mbd-contact`
5. `tpt-mbd-flexible`
6. `tpt-mbd-system`
7. `tpt-mbd` (umbrella)

## Building

This workspace depends on the [`tpt-math`](https://github.com/tpt-solutions/tpt-math) and [`tpt-fem`](https://github.com/tpt-solutions/tpt-fem) substrates, resolved as git dependencies:

```sh
git clone https://github.com/tpt-solutions/tpt-multibody-dynamics.git
cd tpt-multibody-dynamics
cargo build --workspace
cargo test --workspace
```

- **Edition:** `2021`.
- **MSRV:** `1.85` (pinned via `rust-version` in `[workspace.package]`).
- **Toolchain:** fixed by `rust-toolchain.toml` via `rustup`.
- **License:** MIT OR Apache-2.0.

## Status

Scaffolded from `tpt-rust-map/template/`. Build order and crate scope follow `todo.md` and `spec.txt`.

Dual-licensed under [MIT](LICENSE-MIT) or [Apache-2.0](LICENSE-APACHE) at your option.

Copyright (c) 2026 TPT Solutions.
