//! Rigid-body reference frame: origin + unit-quaternion orientation.
//!
//! Re-exports the core geometric primitives from [`tpt_math_geometry`] and
//! [`tpt_math_spatial`] at the [`tpt_mbd_core`] namespace level, and provides
//! [`TransformTree`] for hierarchical body connections.

#[cfg(feature = "alloc")]
pub use tpt_math_geometry::{Isometry3, Quaternion, Rotation3, Translation, UnitQuaternion};
#[cfg(feature = "alloc")]
pub use tpt_math_spatial::{DualQuaternion, Screw};

#[cfg(not(feature = "alloc"))]
pub use tpt_math_geometry::{Isometry3, Quaternion, Rotation3, Translation, UnitQuaternion};
#[cfg(not(feature = "alloc"))]
pub use tpt_math_spatial::{DualQuaternion, Screw};

/// A rigid-body reference frame: origin + unit-quaternion orientation.
pub type Frame = Isometry3<f64>;

// ===========================================================================
// alloc feature: heap-backed TransformTree
// ===========================================================================

#[cfg(feature = "alloc")]
pub struct TransformTree {
    nodes: Vec<TreeNode>,
}

#[cfg(feature = "alloc")]
pub struct TreeNode {
    pub name: &'static str,
    pub transform_from_parent: Isometry3<f64>,
    pub parent: Option<usize>,
    pub children: Vec<usize>,
}

#[cfg(feature = "alloc")]
impl TransformTree {
    pub fn new() -> Self {
        Self { nodes: Vec::new() }
    }

    pub fn add_node(
        &mut self,
        name: &'static str,
        transform_from_parent: Isometry3<f64>,
        parent: Option<usize>,
    ) -> usize {
        let index = self.nodes.len();
        self.nodes.push(TreeNode {
            name,
            transform_from_parent,
            parent,
            children: Vec::new(),
        });
        if let Some(p) = parent {
            self.nodes[p].children.push(index);
        }
        index
    }

    pub fn children(&self, node: usize) -> &[usize] {
        &self.nodes[node].children
    }

    pub fn transform(&self, from_node: usize, to_node: usize) -> Isometry3<f64> {
        let w_from = self.world_transform(from_node);
        let w_to = self.world_transform(to_node);
        w_to * w_from.inverse()
    }

    fn world_transform(&self, node: usize) -> Isometry3<f64> {
        let mut result = Isometry3::identity();
        let mut current = node;
        loop {
            let t = self.nodes[current].transform_from_parent;
            result = t * result;
            if let Some(p) = self.nodes[current].parent {
                current = p;
            } else {
                break;
            }
        }
        result
    }
}

// ===========================================================================
// no_alloc: fixed-capacity TransformTree
// ===========================================================================

#[cfg(not(feature = "alloc"))]
#[derive(Copy, Clone, Debug)]
pub struct TreeNode {
    pub name: &'static str,
    pub transform_from_parent: Isometry3<f64>,
    pub parent: Option<usize>,
    pub children: [usize; 8],
    pub child_count: usize,
}

#[cfg(not(feature = "alloc"))]
#[derive(Copy, Clone, Debug)]
pub struct TransformTree {
    nodes: [TreeNode; 16],
    active: [bool; 16],
    count: usize,
}

#[cfg(not(feature = "alloc"))]
impl TransformTree {
    pub fn new() -> Self {
        Self {
            nodes: [TreeNode {
                name: "",
                transform_from_parent: Isometry3::identity(),
                parent: None,
                children: [0; 8],
                child_count: 0,
            }; 16],
            active: [false; 16],
            count: 0,
        }
    }

    pub fn add_node(
        &mut self,
        name: &'static str,
        transform_from_parent: Isometry3<f64>,
        parent: Option<usize>,
    ) -> usize {
        let index = self.count;
        self.nodes[index] = TreeNode {
            name,
            transform_from_parent,
            parent,
            children: [0; 8],
            child_count: 0,
        };
        self.active[index] = true;
        self.count += 1;
        if let Some(p) = parent {
            let parent_node = &mut self.nodes[p];
            let ci = parent_node.child_count;
            parent_node.children[ci] = index;
            parent_node.child_count += 1;
        }
        index
    }

    pub fn children(&self, node: usize) -> &[usize] {
        let n = &self.nodes[node];
        &n.children[..n.child_count]
    }

    pub fn transform(&self, from_node: usize, to_node: usize) -> Isometry3<f64> {
        let w_from = self.world_transform(from_node);
        let w_to = self.world_transform(to_node);
        w_to * w_from.inverse()
    }

    fn world_transform(&self, node: usize) -> Isometry3<f64> {
        let mut result = Isometry3::identity();
        let mut current = node;
        loop {
            let node_data = &self.nodes[current];
            let t = node_data.transform_from_parent;
            result = t * result;
            if let Some(p) = node_data.parent {
                current = p;
            } else {
                break;
            }
        }
        result
    }
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_math_geometry::{Rotation3, Translation};
    use tpt_math_linalg_fixed::Vector3;

    #[test]
    fn serial_chain_transform_composition() {
        let mut tree = TransformTree::new();

        let rot = Rotation3::<f64>::from_axis_angle(
            &Vector3::new([0.0, 0.0, 1.0]),
            core::f64::consts::FRAC_PI_2,
        );
        let translation = Translation::new(Vector3::new([1.0, 0.0, 0.0]));
        let iso = Isometry3::new(translation, rot);

        let root = tree.add_node("root", Isometry3::identity(), None);
        let link1 = tree.add_node("link1", iso, Some(root));
        let link2 = tree.add_node("link2", iso, Some(link1));
        let link3 = tree.add_node("link3", iso, Some(link2));

        let expected = iso * iso * iso;
        let computed = tree.transform(root, link3);

        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (computed.rotation.matrix().data[i][j] - expected.rotation.matrix().data[i][j])
                        .abs()
                        < 1e-12,
                    "rotation {i},{j}"
                );
            }
            assert!(
                (computed.translation.vector.data[i] - expected.translation.vector.data[i]).abs()
                    < 1e-12,
                "translation {i}"
            );
        }

        assert_eq!(tree.children(link1), &[link2]);
    }
}
