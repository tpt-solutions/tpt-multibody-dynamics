//! Absolute Nodal Coordinate Formulation (ANCF): large-deformation beam/shell
//! elements with gradient-deficient shape functions.
//!
//! ANCF uses absolute nodal coordinates (positions and slopes) so the
//! kinematic description remains valid for arbitrarily large rotations.

use tpt_math_linalg_dense::DMatrix;

/// 2-node ANCF beam element in 2-D (each node carries the transverse
/// displacement and its slope).
#[derive(Clone, Debug)]
pub struct AncfBeam {
    /// Number of elements in the body.
    pub num_elements: usize,
    /// Reference length of a single element.
    pub element_length: f64,
    /// Cross-section parameters: `(cross_section_area, moment_of_inertia)`.
    pub cross_section: (f64, f64),
}

impl AncfBeam {
    /// Build a new ANCF beam description.
    pub fn new(num_elements: usize, element_length: f64, cross_section: (f64, f64)) -> Self {
        AncfBeam {
            num_elements,
            element_length,
            cross_section,
        }
    }

    /// Consistent mass matrix for a single 2-node ANCF beam element.
    ///
    /// Each node has 2 DOFs (transverse displacement and slope), giving a
    /// `4 × 4` element matrix.  The matrix is derived from the consistent
    /// mass definition `M_e = ∫ ρ A Nᵀ N dX` with linear shape functions.
    pub fn element_mass_matrix(&self) -> DMatrix<f64> {
        let l = self.element_length;
        let rho_a = 1.0; // unit mass per unit length for the generic element matrix
        let coeff = rho_a * l / 30.0;
        DMatrix::from_row_slice(
            4,
            4,
            &[
                156.0 * coeff,
                22.0 * l * coeff,
                54.0 * coeff,
                -13.0 * l * coeff,
                22.0 * l * coeff,
                4.0 * l * l * coeff,
                13.0 * l * coeff,
                -3.0 * l * l * coeff,
                54.0 * coeff,
                13.0 * l * coeff,
                156.0 * coeff,
                -22.0 * l * coeff,
                -13.0 * l * coeff,
                -3.0 * l * l * coeff,
                -22.0 * l * coeff,
                4.0 * l * l * coeff,
            ],
        )
    }

    /// Linear stiffness matrix for a single 2-node ANCF beam element.
    ///
    /// `e` is Young's modulus and `i` is the second moment of area of the
    /// cross-section.  The `4 × 4` element matrix corresponds to the Euler-
    /// Bernoulli beam stiffness `K_e = ∫ EI Bᵀ B dX`.
    pub fn element_stiffness_matrix(&self, e: f64, i: f64) -> DMatrix<f64> {
        let l = self.element_length;
        let ei = e * i;
        let coeff = ei / (l * l * l);
        DMatrix::from_row_slice(
            4,
            4,
            &[
                12.0 * coeff,
                6.0 * l * coeff,
                -12.0 * coeff,
                6.0 * l * coeff,
                6.0 * l * coeff,
                4.0 * l * l * coeff,
                -6.0 * l * coeff,
                2.0 * l * l * coeff,
                -12.0 * coeff,
                -6.0 * l * coeff,
                12.0 * coeff,
                -6.0 * l * coeff,
                6.0 * l * coeff,
                2.0 * l * coeff,
                -6.0 * l * coeff,
                4.0 * l * l * coeff,
            ],
        )
    }
}

/// 4-node ANCF shell element (quadrilateral, 2-D surface with 3 DOFs per node).
#[derive(Clone, Debug)]
pub struct AncfShell {
    /// Shell thickness.
    pub thickness: f64,
    /// Young's modulus.
    pub young_modulus: f64,
    /// Poisson's ratio.
    pub poisson_ratio: f64,
}

impl AncfShell {
    /// Build a new ANCF shell description.
    pub fn new(thickness: f64, young_modulus: f64, poisson_ratio: f64) -> Self {
        AncfShell {
            thickness,
            young_modulus,
            poisson_ratio,
        }
    }
}

/// Supported ANCF element topologies.
#[derive(Clone, Debug)]
pub enum AncfElement {
    /// 2-node beam element.
    Beam(AncfBeam),
    /// 4-node shell (membrane) element.
    Plate(AncfShell),
    /// 4-node shell element with transverse shear.
    Shell(AncfShell),
}

impl AncfElement {
    /// Return the number of nodes for this element type.
    pub fn node_count(&self) -> usize {
        match self {
            AncfElement::Beam(_) => 2,
            AncfElement::Plate(_) => 4,
            AncfElement::Shell(_) => 4,
        }
    }

    /// Return the degrees of freedom per node for this element type.
    pub fn dofs_per_node(&self) -> usize {
        match self {
            AncfElement::Beam(_) => 2,
            AncfElement::Plate(_) => 3,
            AncfElement::Shell(_) => 3,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_beam_element_mass_matrix_symmetric() {
        let beam = AncfBeam::new(4, 1.0, (1.0, 1.0));
        let m = beam.element_mass_matrix();
        assert_eq!(m.nrows(), 4);
        assert_eq!(m.ncols(), 4);
        for i in 0..4 {
            for j in 0..4 {
                assert!(
                    (m[(i, j)] - m[(j, i)]).abs() < 1e-12,
                    "mass matrix not symmetric at ({},{})",
                    i,
                    j
                );
            }
        }
    }

    #[test]
    fn test_beam_element_stiffness_symmetric() {
        let beam = AncfBeam::new(4, 1.0, (1.0, 1.0));
        let k = beam.element_stiffness_matrix(1.0, 1.0);
        assert_eq!(k.nrows(), 4);
        assert_eq!(k.ncols(), 4);
        for i in 0..4 {
            for j in 0..4 {
                assert!(
                    (k[(i, j)] - k[(j, i)]).abs() < 1e-12,
                    "stiffness matrix not symmetric at ({},{})",
                    i,
                    j
                );
            }
        }
    }

    #[test]
    fn test_ancf_element_node_count() {
        let beam = AncfElement::Beam(AncfBeam::new(1, 1.0, (1.0, 1.0)));
        let plate = AncfElement::Plate(AncfShell::new(0.1, 1.0, 0.3));
        let shell = AncfElement::Shell(AncfShell::new(0.1, 1.0, 0.3));
        assert_eq!(beam.node_count(), 2);
        assert_eq!(plate.node_count(), 4);
        assert_eq!(shell.node_count(), 4);
    }

    #[test]
    fn test_beam_mass_positive_diagonal() {
        let beam = AncfBeam::new(4, 1.0, (1.0, 1.0));
        let m = beam.element_mass_matrix();
        for i in 0..4 {
            assert!(m[(i, i)] > 0.0, "mass diagonal must be positive");
        }
    }

    #[test]
    fn test_beam_stiffness_positive_diagonal() {
        let beam = AncfBeam::new(4, 1.0, (1.0, 1.0));
        let k = beam.element_stiffness_matrix(210e9, 1e-4);
        for i in [0, 2] {
            assert!(k[(i, i)] > 0.0, "stiffness diagonal must be positive");
        }
    }
}
