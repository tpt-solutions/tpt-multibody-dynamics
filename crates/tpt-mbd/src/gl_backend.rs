//! OpenGL backend for multibody system rendering.
//!
//! Provides minimal OpenGL 3.3+ core-profile helpers for rendering
//! `RenderableBody` instances.  This module is intentionally thin:
//! it issues GL draw calls but does not manage a windowing context.
//!
//! The caller is responsible for:
//! - Creating an OpenGL context (e.g., via `glutin`, `glfw`, or `sdl2`)
//! - Setting up projection / view matrices
//! - Running the render loop

use gl::types::*;

use crate::render::{GeometryType, RenderContext, RenderableBody};

// ===========================================================================
// Rendering helpers
// ===========================================================================

/// Initialize the OpenGL backend (call once after context creation).
pub fn init() {
    unsafe {
        gl::Enable(gl::DEPTH_TEST);
        gl::Enable(gl::CULL_FACE);
        gl::CullFace(gl::BACK);
        gl::FrontFace(gl::CCW);
        gl::ClearColor(0.1, 0.1, 0.1, 1.0);
    }
}

/// Render all bodies in the context.
pub fn render(ctx: &RenderContext) {
    unsafe {
        gl::Clear(gl::COLOR_BUFFER_BIT | gl::DEPTH_BUFFER_BIT);
    }

    for body in &ctx.bodies {
        render_body(body);
    }
}

/// Render a single body.
fn render_body(body: &RenderableBody) {
    unsafe {
        gl::PushMatrix();
        gl::Translatef(
            body.position[0] as GLfloat,
            body.position[1] as GLfloat,
            body.position[2] as GLfloat,
        );
        gl::Rotatef(
            body.rotation[3] as GLfloat,
            body.rotation[0] as GLfloat,
            body.rotation[1] as GLfloat,
            body.rotation[2] as GLfloat,
        );

        match body.geometry_type {
            GeometryType::Sphere => draw_sphere(body.geometry_params.get(0).copied().unwrap_or(0.5)),
            GeometryType::Box => draw_box(
                body.geometry_params.get(0).copied().unwrap_or(0.5),
                body.geometry_params.get(1).copied().unwrap_or(0.5),
                body.geometry_params.get(2).copied().unwrap_or(0.5),
            ),
            GeometryType::Cylinder => draw_cylinder(
                body.geometry_params.get(0).copied().unwrap_or(0.5),
                body.geometry_params.get(1).copied().unwrap_or(1.0),
            ),
            GeometryType::Capsule => draw_capsule(
                body.geometry_params.get(0).copied().unwrap_or(0.5),
                body.geometry_params.get(1).copied().unwrap_or(1.0),
            ),
            GeometryType::Plane => draw_plane(),
        }

        gl::PopMatrix();
    }
}

// ===========================================================================
// Primitive draw calls
// ===========================================================================

fn draw_sphere(radius: f64) {
    unsafe {
        let r = radius as GLfloat;
        let slices = 24;
        let stacks = 16;
        for i in 0..stacks {
            let phi1 = std::f64::consts::PI * i as f64 / stacks as f64;
            let phi2 = std::f64::consts::PI * (i + 1) as f64 / stacks as f64;
            gl::Begin(gl::QUAD_STRIP);
            for j in 0..=slices {
                let theta = 2.0 * std::f64::consts::PI * j as f64 / slices as f64;
                let (s1, c1) = phi1.sin_cos();
                let (s2, c2) = phi2.sin_cos();
                let (st, ct) = theta.sin_cos();
                let (st2, ct2) = (theta + 2.0 * std::f64::consts::PI / slices as f64).sin_cos();
                gl::Normal3f(st * c1, ct * c1, s1);
                gl::Vertex3f(r * st * c1, r * ct * c1, r * s1);
                gl::Normal3f(st * c2, ct * c2, s2);
                gl::Vertex3f(r * st * c2, r * ct * c2, r * s2);
            }
            gl::End();
        }
    }
}

fn draw_box(hx: f64, hy: f64, hz: f64) {
    unsafe {
        let x = hx as GLfloat;
        let y = hy as GLfloat;
        let z = hz as GLfloat;
        gl::Begin(gl::QUADS);
        // front
        gl::Normal3f(0.0, 0.0, 1.0);
        gl::Vertex3f(-x, -y, z);
        gl::Vertex3f(x, -y, z);
        gl::Vertex3f(x, y, z);
        gl::Vertex3f(-x, y, z);
        // back
        gl::Normal3f(0.0, 0.0, -1.0);
        gl::Vertex3f(-x, -y, -z);
        gl::Vertex3f(-x, y, -z);
        gl::Vertex3f(x, y, -z);
        gl::Vertex3f(x, -y, -z);
        // top
        gl::Normal3f(0.0, 1.0, 0.0);
        gl::Vertex3f(-x, y, -z);
        gl::Vertex3f(-x, y, z);
        gl::Vertex3f(x, y, z);
        gl::Vertex3f(x, y, -z);
        // bottom
        gl::Normal3f(0.0, -1.0, 0.0);
        gl::Vertex3f(-x, -y, -z);
        gl::Vertex3f(x, -y, -z);
        gl::Vertex3f(x, -y, z);
        gl::Vertex3f(-x, -y, z);
        // right
        gl::Normal3f(1.0, 0.0, 0.0);
        gl::Vertex3f(x, -y, -z);
        gl::Vertex3f(x, y, -z);
        gl::Vertex3f(x, y, z);
        gl::Vertex3f(x, -y, z);
        // left
        gl::Normal3f(-1.0, 0.0, 0.0);
        gl::Vertex3f(-x, -y, -z);
        gl::Vertex3f(-x, -y, z);
        gl::Vertex3f(-x, y, z);
        gl::Vertex3f(-x, y, -z);
        gl::End();
    }
}

fn draw_cylinder(radius: f64, half_length: f64) {
    unsafe {
        let r = radius as GLfloat;
        let h = half_length as GLfloat;
        let slices = 24;
        gl::Begin(gl::QUAD_STRIP);
        for i in 0..=slices {
            let theta = 2.0 * std::f64::consts::PI * i as f64 / slices as f64;
            let (st, ct) = theta.sin_cos();
            gl::Normal3f(st, ct, 0.0);
            gl::Vertex3f(r * st, r * ct, -h);
            gl::Vertex3f(r * st, r * ct, h);
        }
        gl::End();
    }
}

fn draw_capsule(radius: f64, half_length: f64) {
    draw_cylinder(radius, half_length);
    draw_sphere(radius);
}

fn draw_plane() {
    unsafe {
        gl::Begin(gl::QUADS);
        gl::Normal3f(0.0, 1.0, 0.0);
        gl::Vertex3f(-10.0, 0.0, -10.0);
        gl::Vertex3f(10.0, 0.0, -10.0);
        gl::Vertex3f(10.0, 0.0, 10.0);
        gl::Vertex3f(-10.0, 0.0, 10.0);
        gl::End();
    }
}
