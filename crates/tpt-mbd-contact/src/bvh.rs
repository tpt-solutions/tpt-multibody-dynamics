//! Bounding Volume Hierarchy (BVH) for broad-phase collision detection.
//!
//! Provides an AABB tree for efficiently finding potentially colliding pairs
//! of objects. The BVH is a binary tree where each node stores an AABB that
//! tightly encloses all its children's AABBs.

extern crate alloc;

use alloc::vec::Vec;

use crate::{
    detection::{aabb_overlap, BroadPhaseAABB},
    gjk::{gjk_distance, SupportMap},
    Vector3,
};

// ===========================================================================
// BVH Node and Tree
// ===========================================================================

/// A node in the bounding volume hierarchy.
#[derive(Clone, Debug)]
pub struct BvhNode {
    /// Axis-aligned bounding box enclosing this node and all children.
    pub aabb: BroadPhaseAABB,
    /// Object index (leaf nodes) or None (internal nodes).
    pub object_index: Option<usize>,
    /// Left child index (internal nodes), or None (leaf nodes).
    pub left: Option<usize>,
    /// Right child index (internal nodes), or None (leaf nodes).
    pub right: Option<usize>,
    /// Parent index, or None for the root.
    pub parent: Option<usize>,
}

impl BvhNode {
    /// Create a new leaf node.
    pub fn leaf(object_index: usize, aabb: BroadPhaseAABB) -> Self {
        BvhNode {
            aabb,
            object_index: Some(object_index),
            left: None,
            right: None,
            parent: None,
        }
    }

    /// Create a new internal node.
    pub fn internal(aabb: BroadPhaseAABB, left: usize, right: usize) -> Self {
        BvhNode {
            aabb,
            object_index: None,
            left: Some(left),
            right: Some(right),
            parent: None,
        }
    }

    /// Returns true if this is a leaf node.
    pub fn is_leaf(&self) -> bool {
        self.object_index.is_some()
    }
}

/// Bounding volume hierarchy over AABBs.
#[derive(Clone, Debug, Default)]
pub struct BvhTree {
    /// Nodes of the tree.
    pub nodes: Vec<BvhNode>,
    /// Root node index.
    pub root: Option<usize>,
    /// Object AABBs indexed by object index.
    pub object_aabbs: Vec<BroadPhaseAABB>,
}

impl BvhTree {
    /// Create an empty BVH tree.
    pub fn new() -> Self {
        Self::default()
    }

    /// Build a BVH tree from object AABBs using a top-down SAH (surface area
    /// heuristic) construction.
    ///
    /// `aabbs` is a list of AABBs, one per object. The tree is rebuilt from
    /// scratch.
    pub fn build(aabbs: &[BroadPhaseAABB]) -> Self {
        let mut tree = BvhTree {
            nodes: Vec::new(),
            root: None,
            object_aabbs: aabbs.to_vec(),
        };

        if aabbs.is_empty() {
            return tree;
        }

        let mut indices: Vec<usize> = (0..aabbs.len()).collect();
        let root = tree.build_recursive(&mut indices, aabbs, None);
        tree.root = Some(root);
        tree
    }

    /// Recursively build BVH nodes using median split along the widest axis.
    fn build_recursive(
        &mut self,
        indices: &mut [usize],
        aabbs: &[BroadPhaseAABB],
        parent: Option<usize>,
    ) -> usize {
        if indices.len() == 1 {
            let idx = indices[0];
            let aabb = aabbs[idx];
            let node_idx = self.nodes.len();
            self.nodes.push(BvhNode::leaf(idx, aabb));
            self.nodes[node_idx].parent = parent;
            return node_idx;
        }

        let mut centroid_aabb_min = [f64::MAX; 3];
        let mut centroid_aabb_max = [f64::MIN; 3];
        for &idx in indices.iter() {
            let c = aabbs[idx].center();
            for i in 0..3 {
                centroid_aabb_min[i] = centroid_aabb_min[i].min(c.0[i]);
                centroid_aabb_max[i] = centroid_aabb_max[i].max(c.0[i]);
            }
        }

        let mut best_axis = 0usize;
        let mut best_extent = 0.0f64;
        for axis in 0..3 {
            let extent = centroid_aabb_max[axis] - centroid_aabb_min[axis];
            if extent > best_extent {
                best_extent = extent;
                best_axis = axis;
            }
        }

        indices.sort_by(|&a, &b| {
            let ca = aabbs[a].center().0[best_axis];
            let cb = aabbs[b].center().0[best_axis];
            ca.partial_cmp(&cb).unwrap_or(core::cmp::Ordering::Equal)
        });

        let mid = indices.len() / 2;
        let (left_indices, right_indices) = indices.split_at_mut(mid);
        let left = self.build_recursive(left_indices, aabbs, None);
        let right = self.build_recursive(right_indices, aabbs, None);

        let left_aabb = &self.nodes[left].aabb;
        let right_aabb = &self.nodes[right].aabb;
        let parent_aabb = BroadPhaseAABB::new(
            Vector3::new(
                left_aabb.min.0[0].min(right_aabb.min.0[0]),
                left_aabb.min.0[1].min(right_aabb.min.0[1]),
                left_aabb.min.0[2].min(right_aabb.min.0[2]),
            ),
            Vector3::new(
                left_aabb.max.0[0].max(right_aabb.max.0[0]),
                left_aabb.max.0[1].max(right_aabb.max.0[1]),
                left_aabb.max.0[2].max(right_aabb.max.0[2]),
            ),
        );

        let node_idx = self.nodes.len();
        self.nodes.push(BvhNode::internal(parent_aabb, left, right));
        self.nodes[node_idx].parent = parent;
        self.nodes[left].parent = Some(node_idx);
        self.nodes[right].parent = Some(node_idx);
        node_idx
    }

    /// Compute the surface area of an AABB.
    pub fn surface_area(aabb: &BroadPhaseAABB) -> f64 {
        let dx = aabb.max.0[0] - aabb.min.0[0];
        let dy = aabb.max.0[1] - aabb.min.0[1];
        let dz = aabb.max.0[2] - aabb.min.0[2];
        2.0 * (dx * dy + dy * dz + dz * dx)
    }

    /// Find all potentially overlapping pairs by traversing the BVH.
    ///
    /// Returns a list of `(i, j)` pairs where `i < j` and the AABBs may overlap.
    /// Further narrow-phase collision detection is needed to confirm contact.
    pub fn find_potential_pairs(&self) -> Vec<(usize, usize)> {
        let mut pairs = Vec::new();
        if let Some(root) = self.root {
            self.find_pairs_recursive(root, &mut pairs);
        }
        pairs
    }

    /// Recursively find overlapping pairs.
    fn find_pairs_recursive(&self, node_idx: usize, pairs: &mut Vec<(usize, usize)>) {
        let node = &self.nodes[node_idx];

        if let (Some(left_idx), Some(right_idx)) = (node.left, node.right) {
            let left = &self.nodes[left_idx];
            let right = &self.nodes[right_idx];

            if aabb_overlap(&left.aabb, &right.aabb) {
                if left.is_leaf() && right.is_leaf() {
                    if let (Some(li), Some(ri)) = (left.object_index, right.object_index) {
                        pairs.push((li.min(ri), li.max(ri)));
                    }
                } else if left.is_leaf() {
                    self.leaf_vs_subtree(right_idx, left.object_index.unwrap(), pairs);
                } else if right.is_leaf() {
                    self.leaf_vs_subtree(left_idx, right.object_index.unwrap(), pairs);
                } else {
                    self.find_pairs_recursive(left_idx, pairs);
                    self.find_pairs_recursive(right_idx, pairs);
                }
            }

            self.find_pairs_recursive(left_idx, pairs);
            self.find_pairs_recursive(right_idx, pairs);
        }
    }

    /// Test a leaf object against all leaves in a subtree.
    fn leaf_vs_subtree(&self, node_idx: usize, obj_idx: usize, pairs: &mut Vec<(usize, usize)>) {
        let node = &self.nodes[node_idx];

        if node.is_leaf() {
            if let Some(other_idx) = node.object_index {
                if other_idx != obj_idx {
                    pairs.push((obj_idx.min(other_idx), obj_idx.max(other_idx)));
                }
            }
            return;
        }

        if let Some(left_idx) = node.left {
            let left = &self.nodes[left_idx];
            if aabb_overlap(&left.aabb, &self.object_aabbs[obj_idx]) {
                self.leaf_vs_subtree(left_idx, obj_idx, pairs);
            }
        }

        if let Some(right_idx) = node.right {
            let right = &self.nodes[right_idx];
            if aabb_overlap(&right.aabb, &self.object_aabbs[obj_idx]) {
                self.leaf_vs_subtree(right_idx, obj_idx, pairs);
            }
        }
    }

    /// Rebuild the BVH after object AABBs have changed.
    pub fn refit(&mut self) {
        if let Some(root) = self.root {
            self.refit_recursive(root);
        }
    }

    /// Refit a node's AABB from its children.
    fn refit_recursive(&mut self, node_idx: usize) -> BroadPhaseAABB {
        let (left_idx, right_idx, is_leaf, obj_idx) = {
            let node = &self.nodes[node_idx];
            (node.left, node.right, node.is_leaf(), node.object_index)
        };

        if is_leaf {
            if let Some(obj_idx) = obj_idx {
                self.nodes[node_idx].aabb = self.object_aabbs[obj_idx];
            }
            return self.nodes[node_idx].aabb;
        }

        let mut left_aabb = BroadPhaseAABB::new(Vector3::zero(), Vector3::zero());
        let mut right_aabb = BroadPhaseAABB::new(Vector3::zero(), Vector3::zero());

        if let Some(left_idx) = left_idx {
            left_aabb = self.refit_recursive(left_idx);
        }
        if let Some(right_idx) = right_idx {
            right_aabb = self.refit_recursive(right_idx);
        }

        self.nodes[node_idx].aabb = BroadPhaseAABB::new(
            Vector3::new(
                left_aabb.min.0[0].min(right_aabb.min.0[0]),
                left_aabb.min.0[1].min(right_aabb.min.0[1]),
                left_aabb.min.0[2].min(right_aabb.min.0[2]),
            ),
            Vector3::new(
                left_aabb.max.0[0].max(right_aabb.max.0[0]),
                left_aabb.max.0[1].max(right_aabb.max.0[1]),
                left_aabb.max.0[2].max(right_aabb.max.0[2]),
            ),
        );

        self.nodes[node_idx].aabb
    }
}

// ===========================================================================
// Collision query helpers
// ===========================================================================

/// Query the BVH for all object pairs whose AABBs overlap.
///
/// This is the broad-phase result. Use narrow-phase (GJK/EPA) to confirm
/// actual contact.
pub fn bvh_broad_phase(aabbs: &[BroadPhaseAABB]) -> Vec<(usize, usize)> {
    let tree = BvhTree::build(aabbs);
    tree.find_potential_pairs()
}

/// Narrow-phase collision test between two objects with given AABBs and
/// support maps.
///
/// Returns `Some(GjkResult)` if the shapes are detected, `None` otherwise.
pub fn narrow_phase_test(
    a: &dyn SupportMap,
    b: &dyn SupportMap,
    aabb_a: &BroadPhaseAABB,
    aabb_b: &BroadPhaseAABB,
) -> Option<crate::gjk::GjkResult> {
    if !aabb_overlap(aabb_a, aabb_b) {
        return None;
    }
    let sep = aabb_a.center() - aabb_b.center();
    if sep.norm() < 1e-12 {
        return Some(gjk_distance(a, b, &Vector3::new(1.0, 0.0, 0.0)));
    }
    Some(gjk_distance(a, b, &sep))
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bvh_empty() {
        let tree = BvhTree::new();
        assert!(tree.root.is_none());
        assert!(tree.find_potential_pairs().is_empty());
    }

    #[test]
    fn test_bvh_single_object() {
        let aabbs = vec![BroadPhaseAABB::new(
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 1.0),
        )];
        let tree = BvhTree::build(&aabbs);
        assert!(tree.root.is_some());
        assert!(tree.nodes[tree.root.unwrap()].is_leaf());
        assert!(tree.find_potential_pairs().is_empty());
    }

    #[test]
    fn test_bvh_two_non_overlapping() {
        let aabbs = vec![
            BroadPhaseAABB::new(Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 1.0, 1.0)),
            BroadPhaseAABB::new(Vector3::new(5.0, 5.0, 5.0), Vector3::new(6.0, 6.0, 6.0)),
        ];
        let tree = BvhTree::build(&aabbs);
        assert!(tree.find_potential_pairs().is_empty());
    }

    #[test]
    fn test_bvh_two_overlapping() {
        let aabbs = vec![
            BroadPhaseAABB::new(Vector3::new(0.0, 0.0, 0.0), Vector3::new(2.0, 2.0, 2.0)),
            BroadPhaseAABB::new(Vector3::new(1.0, 1.0, 1.0), Vector3::new(3.0, 3.0, 3.0)),
        ];
        let pairs = bvh_broad_phase(&aabbs);
        assert!(pairs.contains(&(0, 1)));
    }

    #[test]
    fn test_bvh_multiple_objects() {
        let aabbs = vec![
            BroadPhaseAABB::new(Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 1.0, 1.0)),
            BroadPhaseAABB::new(Vector3::new(0.5, 0.5, 0.5), Vector3::new(1.5, 1.5, 1.5)),
            BroadPhaseAABB::new(Vector3::new(2.0, 2.0, 2.0), Vector3::new(3.0, 3.0, 3.0)),
            BroadPhaseAABB::new(Vector3::new(5.0, 5.0, 5.0), Vector3::new(6.0, 6.0, 6.0)),
        ];
        let pairs = bvh_broad_phase(&aabbs);
        assert!(pairs.contains(&(0, 1)));
        assert!(!pairs.contains(&(0, 3)));
        assert!(!pairs.contains(&(2, 3)));
    }

    #[test]
    fn test_bvh_refit() {
        let aabbs = vec![
            BroadPhaseAABB::new(Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 1.0, 1.0)),
            BroadPhaseAABB::new(Vector3::new(2.0, 2.0, 2.0), Vector3::new(3.0, 3.0, 3.0)),
        ];
        let mut tree = BvhTree::build(&aabbs);
        assert!(tree.root.is_some());
        tree.refit();
        assert!(tree.root.is_some());
    }
}
