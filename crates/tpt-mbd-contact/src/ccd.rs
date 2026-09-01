//! Continuous collision detection: conservative advancement and speculative contacts.

use crate::{
    detection::{self, BroadPhaseAABB, ContactPoint},
    Vector3,
};

/// Result of a continuous collision detection query.
#[derive(Clone, Debug, PartialEq)]
pub enum CcdResult {
    /// No collision within the time step.
    NoHit,
    /// Collision detected at time of impact `toi` with contact information.
    Hit(f64, ContactPoint),
}

fn binary_search_toi(
    t0: f64,
    t1: f64,
    rel_vel: Vector3,
    aabb_i: &BroadPhaseAABB,
    aabb_j: &BroadPhaseAABB,
) -> f64 {
    let mut lo = t0;
    let mut hi = t1;
    for _ in 0..24 {
        let mid = (lo + hi) * 0.5;
        let offset = rel_vel * (mid - t0);
        let moved_a = BroadPhaseAABB {
            min: aabb_i.min + offset,
            max: aabb_i.max + offset,
        };
        if detection::aabb_overlap(&moved_a, aabb_j) {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    hi
}

/// Performs conservative advancement between two moving AABBs.
///
/// Returns the earliest time of impact within `dt`, or `NoHit` if no
/// collision occurs.
pub fn conservative_advancement(
    vel_i: Vector3,
    vel_j: Vector3,
    aabb_i: &BroadPhaseAABB,
    aabb_j: &BroadPhaseAABB,
    dt: f64,
) -> CcdResult {
    let rel_vel = vel_i - vel_j;
    let mut t = 0.0;
    let mut current_a = *aabb_i;

    while t < dt {
        if detection::aabb_overlap(&current_a, aabb_j) {
            let toi = binary_search_toi(0.0, t, rel_vel, aabb_i, aabb_j);
            return CcdResult::Hit(
                toi.max(0.0),
                ContactPoint {
                    position: Vector3::zero(),
                    normal: Vector3::new(0.0, 0.0, 0.0),
                    penetration_depth: 0.0,
                },
            );
        }
        let step = dt * 0.1;
        t += step;
        let offset = rel_vel * step;
        current_a.min = current_a.min + offset;
        current_a.max = current_a.max + offset;
    }
    CcdResult::NoHit
}

/// Returns a speculative contact penetration depth if bodies are closing,
/// or None if they are separating.
///
/// A negative `separation_velocity` indicates the bodies are approaching.
pub fn speculative_contact(dt: f64, separation_velocity: f64) -> Option<f64> {
    if separation_velocity < 0.0 {
        Some((-separation_velocity) * dt)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::detection::BroadPhaseAABB;

    #[test]
    fn test_ccd_no_hit_when_separating() {
        let aabb_i = BroadPhaseAABB::new(Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 1.0, 1.0));
        let aabb_j = BroadPhaseAABB::new(Vector3::new(3.0, 3.0, 3.0), Vector3::new(4.0, 4.0, 4.0));
        let vel_i = Vector3::new(1.0, 0.0, 0.0);
        let vel_j = Vector3::new(0.0, 0.0, 0.0);
        let result = conservative_advancement(vel_i, vel_j, &aabb_i, &aabb_j, 1.0);
        assert_eq!(result, CcdResult::NoHit);
    }

    #[test]
    fn test_ccd_hit_construction() {
        let cp = ContactPoint {
            position: Vector3::new(1.0, 0.0, 0.0),
            normal: Vector3::new(0.0, 1.0, 0.0),
            penetration_depth: 0.1,
        };
        let result = CcdResult::Hit(0.005, cp);
        match result {
            CcdResult::Hit(toi, _) => assert!(toi > 0.0 && toi < 1.0),
            CcdResult::NoHit => panic!("expected hit"),
        }
    }

    #[test]
    fn test_ccd_hit_returns_valid_toi() {
        let aabb_i = BroadPhaseAABB::new(Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 1.0, 1.0));
        let aabb_j = BroadPhaseAABB::new(Vector3::new(2.0, 0.0, 0.0), Vector3::new(3.0, 1.0, 1.0));
        let vel_i = Vector3::new(2.0, 0.0, 0.0);
        let vel_j = Vector3::new(0.0, 0.0, 0.0);
        let result = conservative_advancement(vel_i, vel_j, &aabb_i, &aabb_j, 1.0);
        match result {
            CcdResult::Hit(toi, _) => assert!(toi > 0.0 && toi < 1.0),
            CcdResult::NoHit => panic!("expected hit"),
        }
    }

    #[test]
    fn test_speculative_contact_closing() {
        let depth = speculative_contact(0.016, -1.0);
        assert!(depth.is_some());
        assert!((depth.unwrap() - 0.016).abs() < 1e-10);
    }

    #[test]
    fn test_speculative_contact_separating() {
        let depth = speculative_contact(0.016, 1.0);
        assert!(depth.is_none());
    }

    #[test]
    fn test_speculative_contact_zero_velocity() {
        let depth = speculative_contact(0.016, 0.0);
        assert!(depth.is_none());
    }
}
