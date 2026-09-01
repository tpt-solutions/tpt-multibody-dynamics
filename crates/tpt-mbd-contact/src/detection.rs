//! Discrete collision detection primitives: AABB overlap tests,
//! contact points, and contact manifolds for broad/narrow phase.

extern crate alloc;

use alloc::vec::Vec;

use crate::Vector3;

/// A single contact point between two bodies.
#[derive(Clone, Debug, PartialEq)]
pub struct ContactPoint {
    /// World-space position of the contact point.
    pub position: Vector3,
    /// Contact normal pointing from body_j toward body_i.
    pub normal: Vector3,
    /// Penetration depth (positive means interpenetration).
    pub penetration_depth: f64,
}

/// A set of contact points between two bodies.
#[derive(Clone, Debug, Default)]
pub struct ContactManifold {
    /// Contact points in this manifold.
    pub points: Vec<ContactPoint>,
    /// Index of the first body.
    pub body_i: usize,
    /// Index of the second body.
    pub body_j: usize,
}

/// Axis-aligned bounding box for broad-phase collision detection.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BroadPhaseAABB {
    /// Minimum corner of the bounding box.
    pub min: Vector3,
    /// Maximum corner of the bounding box.
    pub max: Vector3,
}

impl BroadPhaseAABB {
    /// Create a new AABB from min and max corners.
    pub fn new(min: Vector3, max: Vector3) -> Self {
        Self { min, max }
    }

    /// Center of the AABB.
    pub fn center(&self) -> Vector3 {
        Vector3::new(
            (self.min.0[0] + self.max.0[0]) * 0.5,
            (self.min.0[1] + self.max.0[1]) * 0.5,
            (self.min.0[2] + self.max.0[2]) * 0.5,
        )
    }

    /// Half extents of the AABB.
    pub fn half_extents(&self) -> Vector3 {
        Vector3::new(
            (self.max.0[0] - self.min.0[0]) * 0.5,
            (self.max.0[1] - self.min.0[1]) * 0.5,
            (self.max.0[2] - self.min.0[2]) * 0.5,
        )
    }
}

/// Returns true if two AABBs overlap.
pub fn aabb_overlap(a: &BroadPhaseAABB, b: &BroadPhaseAABB) -> bool {
    for i in 0..3 {
        if a.max.0[i] < b.min.0[i] || a.min.0[i] > b.max.0[i] {
            return false;
        }
    }
    true
}

/// Returns the intersection AABB if two AABBs overlap, or None.
pub fn aabb_intersect(a: &BroadPhaseAABB, b: &BroadPhaseAABB) -> Option<BroadPhaseAABB> {
    if !aabb_overlap(a, b) {
        return None;
    }
    let min = Vector3::new(
        a.min.0[0].max(b.min.0[0]),
        a.min.0[1].max(b.min.0[1]),
        a.min.0[2].max(b.min.0[2]),
    );
    let max = Vector3::new(
        a.max.0[0].min(b.max.0[0]),
        a.max.0[1].min(b.max.0[1]),
        a.max.0[2].min(b.max.0[2]),
    );
    Some(BroadPhaseAABB::new(min, max))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_aabb_overlap_disjoint() {
        let a = BroadPhaseAABB::new(Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 1.0, 1.0));
        let b = BroadPhaseAABB::new(Vector3::new(2.0, 2.0, 2.0), Vector3::new(3.0, 3.0, 3.0));
        assert!(!aabb_overlap(&a, &b));
    }

    #[test]
    fn test_aabb_overlap_contained() {
        let a = BroadPhaseAABB::new(Vector3::new(0.0, 0.0, 0.0), Vector3::new(2.0, 2.0, 2.0));
        let b = BroadPhaseAABB::new(Vector3::new(0.5, 0.5, 0.5), Vector3::new(1.5, 1.5, 1.5));
        assert!(aabb_overlap(&a, &b));
    }

    #[test]
    fn test_aabb_overlap_partial() {
        let a = BroadPhaseAABB::new(Vector3::new(0.0, 0.0, 0.0), Vector3::new(2.0, 2.0, 2.0));
        let b = BroadPhaseAABB::new(Vector3::new(1.0, 1.0, 1.0), Vector3::new(3.0, 3.0, 3.0));
        assert!(aabb_overlap(&a, &b));
    }

    #[test]
    fn test_aabb_intersect_returns_intersection() {
        let a = BroadPhaseAABB::new(Vector3::new(0.0, 0.0, 0.0), Vector3::new(2.0, 2.0, 2.0));
        let b = BroadPhaseAABB::new(Vector3::new(1.0, 1.0, 1.0), Vector3::new(3.0, 3.0, 3.0));
        let inter = aabb_intersect(&a, &b).unwrap();
        assert_eq!(inter.min, Vector3::new(1.0, 1.0, 1.0));
        assert_eq!(inter.max, Vector3::new(2.0, 2.0, 2.0));
    }

    #[test]
    fn test_aabb_intersect_no_overlap() {
        let a = BroadPhaseAABB::new(Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 1.0, 1.0));
        let b = BroadPhaseAABB::new(Vector3::new(2.0, 2.0, 2.0), Vector3::new(3.0, 3.0, 3.0));
        assert!(aabb_intersect(&a, &b).is_none());
    }

    #[test]
    fn test_aabb_intersect_identical() {
        let a = BroadPhaseAABB::new(Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 1.0, 1.0));
        let inter = aabb_intersect(&a, &a).unwrap();
        assert_eq!(inter.min, Vector3::new(0.0, 0.0, 0.0));
        assert_eq!(inter.max, Vector3::new(1.0, 1.0, 1.0));
    }

    #[test]
    fn test_contact_point_creation() {
        let cp = ContactPoint {
            position: Vector3::new(1.0, 0.5, 0.0),
            normal: Vector3::new(0.0, 1.0, 0.0),
            penetration_depth: 0.1,
        };
        assert_eq!(cp.penetration_depth, 0.1);
        assert_eq!(cp.normal.0[1], 1.0);
        assert_eq!(cp.position.0[0], 1.0);
    }

    #[test]
    fn test_contact_manifold_default() {
        let manifold = ContactManifold::default();
        assert!(manifold.points.is_empty());
        assert_eq!(manifold.body_i, 0);
        assert_eq!(manifold.body_j, 0);
    }

    #[test]
    fn test_contact_manifold_add_point() {
        let mut manifold = ContactManifold::default();
        manifold.points.push(ContactPoint {
            position: Vector3::new(0.5, 0.5, 0.0),
            normal: Vector3::new(0.0, 1.0, 0.0),
            penetration_depth: 0.1,
        });
        assert_eq!(manifold.points.len(), 1);
        assert_eq!(manifold.points[0].penetration_depth, 0.1);
    }

    #[test]
    fn test_broad_phase_aabb_center() {
        let aabb = BroadPhaseAABB::new(Vector3::new(0.0, 0.0, 0.0), Vector3::new(2.0, 2.0, 2.0));
        assert_eq!(aabb.center(), Vector3::new(1.0, 1.0, 1.0));
    }

    #[test]
    fn test_broad_phase_aabb_half_extents() {
        let aabb = BroadPhaseAABB::new(Vector3::new(-1.0, -1.0, -1.0), Vector3::new(1.0, 1.0, 1.0));
        assert_eq!(aabb.half_extents(), Vector3::new(1.0, 1.0, 1.0));
    }
}
