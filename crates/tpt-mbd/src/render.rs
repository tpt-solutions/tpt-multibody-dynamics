//! OpenGL rendering hooks for real-time multibody animation.
//!
//! This module provides a minimal OpenGL rendering interface for visualizing
//! multibody systems in real time. It is intentionally kept as a lightweight
//! abstraction so it can be replaced with a more sophisticated renderer
//! (e.g., `wgpu`, `three-rs`) without affecting the simulation core.
//!
//! # Examples
//!
//! ```ignore
//! use tpt_mbd::render::{RenderContext, RenderableBody};
//!
//! let ctx = RenderContext::new();
//! ```

#[cfg(feature = "gl")]
pub mod gl_backend {
    //! OpenGL backend for multibody rendering.
}

/// A renderable body in the scene.
#[derive(Clone, Debug, Default)]
pub struct RenderableBody {
    /// Body name.
    pub name: String,
    /// World-space position.
    pub position: [f64; 3],
    /// World-space rotation (quaternion: x, y, z, w).
    pub rotation: [f64; 4],
    /// Geometry type identifier.
    pub geometry_type: GeometryType,
    /// Geometry parameters (radius for sphere, half-extents for box, etc.).
    pub geometry_params: Vec<f64>,
}

/// Supported geometry types for rendering.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum GeometryType {
    /// Sphere with given radius.
    #[default]
    Sphere,
    /// Box with given half-extents.
    Box,
    /// Cylinder with given radius and half-length.
    Cylinder,
    /// Capsule with given radius and half-length.
    Capsule,
    /// Plane (infinite).
    Plane,
}

impl GeometryType {
    /// Number of geometry parameters required for this type.
    pub fn param_count(&self) -> usize {
        match self {
            GeometryType::Sphere => 1,
            GeometryType::Box => 3,
            GeometryType::Cylinder => 2,
            GeometryType::Capsule => 2,
            GeometryType::Plane => 0,
        }
    }
}

/// Minimal renderer context for multibody systems.
#[derive(Clone, Debug, Default)]
pub struct RenderContext {
    /// List of bodies to render.
    pub bodies: Vec<RenderableBody>,
    /// Background color (r, g, b).
    pub background: [f64; 3],
    /// Camera position.
    pub camera_position: [f64; 3],
    /// Camera target.
    pub camera_target: [f64; 3],
}

impl RenderContext {
    /// Create a new empty render context.
    pub fn new() -> Self {
        Self {
            bodies: Vec::new(),
            background: [0.1, 0.1, 0.1],
            camera_position: [5.0, 5.0, 5.0],
            camera_target: [0.0, 0.0, 0.0],
        }
    }

    /// Add a body to the render context.
    pub fn add_body(&mut self, body: RenderableBody) {
        self.bodies.push(body);
    }

    /// Clear all bodies from the render context.
    pub fn clear_bodies(&mut self) {
        self.bodies.clear();
    }

    /// Set the camera view.
    pub fn set_camera(&mut self, position: [f64; 3], target: [f64; 3]) {
        self.camera_position = position;
        self.camera_target = target;
    }

    /// Set the background color.
    pub fn set_background(&mut self, r: f64, g: f64, b: f64) {
        self.background = [r, g, b];
    }

    /// Update body poses from simulation state.
    pub fn update_body_pose(&mut self, index: usize, position: [f64; 3], rotation: [f64; 4]) {
        if let Some(body) = self.bodies.get_mut(index) {
            body.position = position;
            body.rotation = rotation;
        }
    }
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_render_context_new() {
        let ctx = RenderContext::new();
        assert!(ctx.bodies.is_empty());
        assert_eq!(ctx.background, [0.1, 0.1, 0.1]);
    }

    #[test]
    fn test_render_context_add_body() {
        let mut ctx = RenderContext::new();
        let body = RenderableBody {
            name: "link0".to_string(),
            position: [0.0, 0.0, 0.0],
            rotation: [0.0, 0.0, 0.0, 1.0],
            geometry_type: GeometryType::Box,
            geometry_params: vec![0.5, 0.5, 0.5],
        };
        ctx.add_body(body);
        assert_eq!(ctx.bodies.len(), 1);
        assert_eq!(ctx.bodies[0].name, "link0");
    }

    #[test]
    fn test_render_context_update_pose() {
        let mut ctx = RenderContext::new();
        let body = RenderableBody {
            name: "link0".to_string(),
            position: [0.0; 3],
            rotation: [0.0; 4],
            geometry_type: GeometryType::Sphere,
            geometry_params: vec![0.5],
        };
        ctx.add_body(body);
        ctx.update_body_pose(0, [1.0, 2.0, 3.0], [0.0, 0.0, 0.0, 1.0]);
        assert_eq!(ctx.bodies[0].position, [1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_geometry_type_param_count() {
        assert_eq!(GeometryType::Sphere.param_count(), 1);
        assert_eq!(GeometryType::Box.param_count(), 3);
        assert_eq!(GeometryType::Cylinder.param_count(), 2);
        assert_eq!(GeometryType::Plane.param_count(), 0);
    }

    #[test]
    fn test_render_context_set_camera() {
        let mut ctx = RenderContext::new();
        ctx.set_camera([1.0, 2.0, 3.0], [0.0, 0.0, 0.0]);
        assert_eq!(ctx.camera_position, [1.0, 2.0, 3.0]);
        assert_eq!(ctx.camera_target, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_render_context_set_background() {
        let mut ctx = RenderContext::new();
        ctx.set_background(0.5, 0.5, 0.5);
        assert_eq!(ctx.background, [0.5, 0.5, 0.5]);
    }

    #[test]
    fn test_render_context_clear() {
        let mut ctx = RenderContext::new();
        ctx.add_body(RenderableBody {
            name: "b1".into(),
            position: [0.0; 3],
            rotation: [0.0, 0.0, 0.0, 1.0],
            geometry_type: GeometryType::Sphere,
            geometry_params: vec![0.5],
        });
        ctx.clear_bodies();
        assert!(ctx.bodies.is_empty());
    }
}
