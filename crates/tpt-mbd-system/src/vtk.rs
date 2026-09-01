//! VTK legacy-format export for multibody system visualization.
//!
//! Writes body poses as an unstructured grid that ParaView and other VTK
//! consumers can render.

extern crate alloc;

use alloc::string::String;

use crate::system::MultibodySystem;

/// Export a [`MultibodySystem`] as a VTK legacy unstructured-grid file.
///
/// Each rigid body is written as a vertex at its reference-frame origin.
/// Flexible bodies are written as vertices at their reference-frame origin
/// with an additional scalar field `modal_dofs` indicating how many modal
/// coordinates they carry.
///
/// Contact manifolds are written as line cells connecting the two body
/// origins involved in each contact.
///
/// Returns the VTK file content as a UTF-8 string.
pub fn to_vtk(system: &MultibodySystem) -> String {
    let mut out = String::new();

    out.push_str("# vtk DataFile Version 3.0\n");
    out.push_str("tpt-mbd system export\n");
    out.push_str("ASCII\n");
    out.push_str("DATASET UNSTRUCTURED_GRID\n");

    let num_rigid = system.bodies.len();
    let num_flex = system.flexible.len();
    let total_points = num_rigid + num_flex;
    let num_contact_cells = system.contacts.len();

    let total_cells = num_contact_cells;
    let cell_capacity = total_cells * 3; // 2 pts + size prefix per line cell

    out.push_str(&format!("POINTS {} float\n", total_points));

    for body in &system.bodies {
        let t = &body.spatial_inertia.matrix;
        let x = t.data[0][3];
        let y = t.data[1][3];
        let z = t.data[2][3];
        out.push_str(&format!("{} {} {}\n", x, y, z));
    }

    for flex in &system.flexible {
        let t = &flex.rigid_body.spatial_inertia.matrix;
        let x = t.data[0][3];
        let y = t.data[1][3];
        let z = t.data[2][3];
        out.push_str(&format!("{} {} {}\n", x, y, z));
    }

    out.push_str(&format!("CELLS {} {}\n", total_cells, cell_capacity));

    for (idx, contact) in system.contacts.iter().enumerate() {
        let p0 = contact.body_i;
        let p1 = contact.body_j;
        out.push_str(&format!("2 {} {}\n", p0, p1));
        let _ = idx;
    }

    out.push_str(&format!("CELL_TYPES {}\n", total_cells));
    for _ in 0..total_cells {
        out.push_str("3\n");
    }

    out.push_str(&format!("POINT_DATA {}\n", total_points));
    out.push_str("SCALARS body_id int 1\n");
    out.push_str("LOOKUP_TABLE default\n");
    for i in 0..total_points {
        out.push_str(&format!("{}\n", i));
    }

    out.push_str("SCALARS modal_dofs int 1\n");
    out.push_str("LOOKUP_TABLE default\n");
    for _i in 0..num_rigid {
        out.push_str("0\n");
    }
    for flex in &system.flexible {
        out.push_str(&format!("{}\n", flex.num_modes));
    }

    out
}
