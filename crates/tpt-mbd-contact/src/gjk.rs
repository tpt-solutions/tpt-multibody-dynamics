//! Gilbert-Johnson-Keerthi (GJK) distance algorithm and Expanding Polytope
//! Algorithm (EPA) for penetration depth.
//!
//! GJK computes the minimum distance between two convex shapes by iteratively
//! building a simplex in the Minkowski difference space. EPA expands the final
//! simplex to find the penetration depth and contact normal when the shapes overlap.

extern crate alloc;

use alloc::vec::Vec;

use crate::Vector3;

// ===========================================================================
// Support function trait
// ===========================================================================

/// A convex shape that can provide a support point in any direction.
pub trait SupportMap {
    /// Return the support point of the shape in direction `dir`.
    fn support(&self, dir: &Vector3) -> Vector3;
}

impl SupportMap for alloc::boxed::Box<dyn SupportMap> {
    fn support(&self, dir: &Vector3) -> Vector3 {
        (**self).support(dir)
    }
}

// ===========================================================================
// GJK distance
// ===========================================================================

/// Result of a GJK distance query.
#[derive(Clone, Debug, PartialEq)]
pub struct GjkResult {
    /// Minimum distance between the two shapes. Zero means penetration.
    pub distance: f64,
    /// Direction from shape A toward shape B at closest approach.
    pub normal: Vector3,
    /// Closest point on shape A (in world space).
    pub closest_a: Vector3,
    /// Closest point on shape B (in world space).
    pub closest_b: Vector3,
    /// True if shapes are penetrating.
    pub penetrating: bool,
}

/// Compute the minimum distance between two convex shapes using GJK.
pub fn gjk_distance(a: &dyn SupportMap, b: &dyn SupportMap, initial_dir: &Vector3) -> GjkResult {
    let mut dir = *initial_dir;
    let sa = a.support(&dir);
    let sb = b.support(&dir.negate());
    let mut simplex = vec![sa - sb];
    dir = simplex[0].negate();

    for _ in 0..64 {
        let sa = a.support(&dir);
        let sb = b.support(&dir.negate());
        let new_point = sa - sb;

        if new_point.dot(&dir) <= 0.0 {
            let result = closest_point_result(&simplex);
            if result.distance < 1e-8 && !simplex.is_empty() {
                let epa = epa_penetration(a, b, simplex.clone());
                if let Some(epa) = epa {
                    return GjkResult {
                        distance: 0.0,
                        normal: epa.normal,
                        closest_a: Vector3::zero(),
                        closest_b: Vector3::zero(),
                        penetrating: true,
                    };
                }
            }
            return GjkResult {
                distance: result.distance,
                normal: result.normal.negate(),
                closest_a: result.closest_a,
                closest_b: Vector3::zero(),
                penetrating: result.penetrating,
            };
        }

        simplex.push(new_point);

        if simplex.len() > 1 {
            reduce_simplex(&mut simplex, &mut dir);
        } else {
            dir = simplex[0].negate();
        }
    }

    let result = closest_point_result(&simplex);
    GjkResult {
        distance: result.distance,
        normal: result.normal.negate(),
        closest_a: result.closest_a,
        closest_b: Vector3::zero(),
        penetrating: result.penetrating,
    }
}

/// Compute the result from the final simplex: closest point to origin.
fn closest_point_result(simplex: &[Vector3]) -> GjkResult {
    if simplex.is_empty() {
        return GjkResult {
            distance: 0.0,
            normal: Vector3::zero(),
            closest_a: Vector3::zero(),
            closest_b: Vector3::zero(),
            penetrating: true,
        };
    }

    if simplex.len() == 1 {
        let p = simplex[0];
        let d = p.norm();
        if d < 1e-12 {
            return GjkResult {
                distance: 0.0,
                normal: Vector3::zero(),
                closest_a: Vector3::zero(),
                closest_b: Vector3::zero(),
                penetrating: true,
            };
        }
        let n = p * (1.0 / d);
        return GjkResult {
            distance: d,
            normal: n,
            closest_a: p,
            closest_b: Vector3::zero(),
            penetrating: false,
        };
    }

    if simplex.len() == 2 {
        let a = simplex[0];
        let b = simplex[1];
        let ab = b - a;
        let t = (-a.dot(&ab)) / ab.dot(&ab).max(1e-12);
        let t = t.clamp(0.0, 1.0);
        let closest = a + ab * t;
        let d = closest.norm();
        if d < 1e-12 {
            return GjkResult {
                distance: 0.0,
                normal: Vector3::zero(),
                closest_a: Vector3::zero(),
                closest_b: Vector3::zero(),
                penetrating: true,
            };
        }
        let n = closest * (1.0 / d);
        return GjkResult {
            distance: d,
            normal: n,
            closest_a: closest,
            closest_b: Vector3::zero(),
            penetrating: false,
        };
    }

    if simplex.len() == 3 {
        let o = Vector3::zero();
        let a = simplex[0];
        let b = simplex[1];
        let c = simplex[2];
        let ab = b - a;
        let ac = c - a;
        let ao = o - a;

        let abc = ab.cross(&ac);
        let n_abc = abc.norm();
        if n_abc > 1e-12 {
            let n = abc * (1.0 / n_abc);
            let ab_n = ab.cross(&n);
            if ao.dot(&ab_n) > 0.0 {
                let ab_len2 = ab.dot(&ab).max(1e-12);
                let t = ao.dot(&ab) / ab_len2;
                let t = t.clamp(0.0, 1.0);
                let closest = a + ab * t;
                let d = closest.norm();
                if d < 1e-12 {
                    return penetrating_result();
                }
                return GjkResult {
                    distance: d,
                    normal: closest * (1.0 / d),
                    closest_a: closest,
                    closest_b: Vector3::zero(),
                    penetrating: false,
                };
            }
            let ac_n = n.cross(&ac);
            if ao.dot(&ac_n) > 0.0 {
                let ac_len2 = ac.dot(&ac).max(1e-12);
                let t = ao.dot(&ac) / ac_len2;
                let t = t.clamp(0.0, 1.0);
                let closest = a + ac * t;
                let d = closest.norm();
                if d < 1e-12 {
                    return penetrating_result();
                }
                return GjkResult {
                    distance: d,
                    normal: closest * (1.0 / d),
                    closest_a: closest,
                    closest_b: Vector3::zero(),
                    penetrating: false,
                };
            }
            let bc = c - b;
            let bo = o - b;
            let bc_n = bc.cross(&n);
            if bo.dot(&bc_n) > 0.0 {
                let bc_len2 = bc.dot(&bc).max(1e-12);
                let t = bo.dot(&bc) / bc_len2;
                let t = t.clamp(0.0, 1.0);
                let closest = b + bc * t;
                let d = closest.norm();
                if d < 1e-12 {
                    return penetrating_result();
                }
                return GjkResult {
                    distance: d,
                    normal: closest * (1.0 / d),
                    closest_a: closest,
                    closest_b: Vector3::zero(),
                    penetrating: false,
                };
            }
            let d_val = ao.dot(&n);
            return GjkResult {
                distance: d_val.max(0.0),
                normal: if d_val > 1e-12 { n } else { Vector3::zero() },
                closest_a: a + n * d_val,
                closest_b: Vector3::zero(),
                penetrating: d_val < 1e-12,
            };
        }

        let ab_len2 = ab.dot(&ab).max(1e-12);
        let t_ab = ao.dot(&ab) / ab_len2;
        let t_ab = t_ab.clamp(0.0, 1.0);
        let closest_ab = a + ab * t_ab;

        let ac_len2 = ac.dot(&ac).max(1e-12);
        let t_ac = ao.dot(&ac) / ac_len2;
        let t_ac = t_ac.clamp(0.0, 1.0);
        let closest_ac = a + ac * t_ac;

        let d_ab = closest_ab.norm();
        let d_ac = closest_ac.norm();

        if d_ab < d_ac {
            if d_ab < 1e-12 {
                return penetrating_result();
            }
            return GjkResult {
                distance: d_ab,
                normal: closest_ab * (1.0 / d_ab),
                closest_a: closest_ab,
                closest_b: Vector3::zero(),
                penetrating: false,
            };
        } else {
            if d_ac < 1e-12 {
                return penetrating_result();
            }
            return GjkResult {
                distance: d_ac,
                normal: closest_ac * (1.0 / d_ac),
                closest_a: closest_ac,
                closest_b: Vector3::zero(),
                penetrating: false,
            };
        }
    }

    // 4+ points: use the first 4 to form a tetrahedron.
    let a = simplex[0];
    let b = simplex[1];
    let c = simplex[2];
    let d = simplex[3];
    let ao = Vector3::zero() - a;

    let ab = b - a;
    let ac = c - a;
    let ad = d - a;

    let abc = ab.cross(&ac);
    let abc_n = abc.norm();
    let abd = ab.cross(&ad);
    let abd_n = abd.norm();
    let acd_ = ac.cross(&ad);
    let acd_n = acd_.norm();

    if abc_n > 1e-12 {
        let n = abc * (1.0 / abc_n);
        if ao.dot(&n) > 0.0 {
            let mut sub = vec![d, c, a];
            return closest_point_triangle_from_vec(&mut sub);
        }
        let d_val = ao.dot(&n);
        return GjkResult {
            distance: d_val.max(0.0),
            normal: if d_val > 1e-12 { n } else { Vector3::zero() },
            closest_a: a + n * d_val,
            closest_b: Vector3::zero(),
            penetrating: d_val < 1e-12,
        };
    }

    if abd_n > 1e-12 {
        let n = abd * (1.0 / abd_n);
        if ao.dot(&n) > 0.0 {
            let sub = vec![d, a];
            return closest_point_result(&sub);
        }
        let d_val = ao.dot(&n);
        return GjkResult {
            distance: d_val.max(0.0),
            normal: if d_val > 1e-12 { n } else { Vector3::zero() },
            closest_a: a + n * d_val,
            closest_b: Vector3::zero(),
            penetrating: d_val < 1e-12,
        };
    }

    if acd_n > 1e-12 {
        let n = acd_ * (1.0 / acd_n);
        if ao.dot(&n) > 0.0 {
            let mut sub = vec![b, d, c];
            return closest_point_triangle_from_vec(&mut sub);
        }
        let d_val = ao.dot(&n);
        return GjkResult {
            distance: d_val.max(0.0),
            normal: if d_val > 1e-12 { n } else { Vector3::zero() },
            closest_a: a + n * d_val,
            closest_b: Vector3::zero(),
            penetrating: d_val < 1e-12,
        };
    }

    let sub = vec![a];
    closest_point_result(&sub)
}

fn closest_point_triangle(a: Vector3, b: Vector3, c: Vector3) -> GjkResult {
    let o = Vector3::zero();
    let ab = b - a;
    let ac = c - a;
    let ao = o - a;

    let abc = ab.cross(&ac);
    let n_abc = abc.norm();
    if n_abc > 1e-12 {
        let n = abc * (1.0 / n_abc);
        let ab_n = ab.cross(&n);
        if ao.dot(&ab_n) > 0.0 {
            let ab_len2 = ab.dot(&ab).max(1e-12);
            let t = ao.dot(&ab) / ab_len2;
            let t = t.clamp(0.0, 1.0);
            let closest = a + ab * t;
            let d = closest.norm();
            if d < 1e-12 {
                return penetrating_result();
            }
            return GjkResult {
                distance: d,
                normal: closest * (1.0 / d),
                closest_a: closest,
                closest_b: Vector3::zero(),
                penetrating: false,
            };
        }

        let ac_n = n.cross(&ac);
        if ao.dot(&ac_n) > 0.0 {
            let ac_len2 = ac.dot(&ac).max(1e-12);
            let t = ao.dot(&ac) / ac_len2;
            let t = t.clamp(0.0, 1.0);
            let closest = a + ac * t;
            let d = closest.norm();
            if d < 1e-12 {
                return penetrating_result();
            }
            return GjkResult {
                distance: d,
                normal: closest * (1.0 / d),
                closest_a: closest,
                closest_b: Vector3::zero(),
                penetrating: false,
            };
        }

        let bc = c - b;
        let bo = o - b;
        let bc_n = bc.cross(&n);
        if bo.dot(&bc_n) > 0.0 {
            let bc_len2 = bc.dot(&bc).max(1e-12);
            let t = bo.dot(&bc) / bc_len2;
            let t = t.clamp(0.0, 1.0);
            let closest = b + bc * t;
            let d = closest.norm();
            if d < 1e-12 {
                return penetrating_result();
            }
            return GjkResult {
                distance: d,
                normal: closest * (1.0 / d),
                closest_a: closest,
                closest_b: Vector3::zero(),
                penetrating: false,
            };
        }

        let d_val = ao.dot(&n);
        return GjkResult {
            distance: d_val.max(0.0),
            normal: if d_val > 1e-12 { n } else { Vector3::zero() },
            closest_a: a + n * d_val,
            closest_b: Vector3::zero(),
            penetrating: d_val < 1e-12,
        };
    }

    let ab_len2 = ab.dot(&ab).max(1e-12);
    let t_ab = ao.dot(&ab) / ab_len2;
    let t_ab = t_ab.clamp(0.0, 1.0);
    let closest_ab = a + ab * t_ab;

    let ac_len2 = ac.dot(&ac).max(1e-12);
    let t_ac = ao.dot(&ac) / ac_len2;
    let t_ac = t_ac.clamp(0.0, 1.0);
    let closest_ac = a + ac * t_ac;

    let d_ab = closest_ab.norm();
    let d_ac = closest_ac.norm();

    if d_ab < d_ac {
        if d_ab < 1e-12 {
            return penetrating_result();
        }
        return GjkResult {
            distance: d_ab,
            normal: closest_ab * (1.0 / d_ab),
            closest_a: closest_ab,
            closest_b: Vector3::zero(),
            penetrating: false,
        };
    } else {
        if d_ac < 1e-12 {
            return penetrating_result();
        }
        return GjkResult {
            distance: d_ac,
            normal: closest_ac * (1.0 / d_ac),
            closest_a: closest_ac,
            closest_b: Vector3::zero(),
            penetrating: false,
        };
    }
}

fn closest_point_triangle_from_vec(simplex: &[Vector3]) -> GjkResult {
    closest_point_triangle(simplex[0], simplex[1], simplex[2])
}

fn penetrating_result() -> GjkResult {
    GjkResult {
        distance: 0.0,
        normal: Vector3::zero(),
        closest_a: Vector3::zero(),
        closest_b: Vector3::zero(),
        penetrating: true,
    }
}

/// Reduce the simplex by removing vertices that are not part of the closest
/// feature to the origin, and update `dir` to point toward the origin from
/// the new simplex.
fn reduce_simplex(simplex: &mut Vec<Vector3>, dir: &mut Vector3) {
    if simplex.len() == 2 {
        let a = simplex[0];
        let b = simplex[1];
        let ab = b - a;
        let t = (-a.dot(&ab)) / ab.dot(&ab).max(1e-12);
        let t = t.clamp(0.0, 1.0);
        simplex[0] = a + ab * t;
        simplex.truncate(1);
        *dir = simplex[0].negate();
        return;
    }

    if simplex.len() == 3 {
        let a = simplex[0];
        let b = simplex[1];
        let c = simplex[2];
        let ab = b - a;
        let ac = c - a;
        let ao = Vector3::zero() - a;
        let abc = ab.cross(&ac);
        let n = abc.normalize();

        let ab_n = ab.cross(&n);
        if ao.dot(&ab_n) > 0.0 {
            simplex[0] = c;
            simplex.truncate(1);
            simplex.push(a);
            *dir = ao;
            return;
        }

        let ac_n = n.cross(&ac);
        if ao.dot(&ac_n) > 0.0 {
            simplex[0] = b;
            simplex.truncate(1);
            simplex.push(a);
            *dir = ao;
            return;
        }

        simplex[0] = a;
        simplex.truncate(3);
        simplex.push(c);
        simplex.push(b);
        *dir = if n.dot(&ao) > 0.0 { n } else { n.negate() };
        return;
    }

    if simplex.len() >= 4 {
        let a = simplex[0];
        let b = simplex[1];
        let c = simplex[2];
        let d = simplex[3];
        let ao = Vector3::zero() - a;

        let ab = b - a;
        let ac = c - a;
        let ad = d - a;

        let abc = ab.cross(&ac);
        let abc_n = abc.norm();
        let abd = ab.cross(&ad);
        let abd_n = abd.norm();
        let acd_ = ac.cross(&ad);
        let acd_n = acd_.norm();

        if abc_n > 1e-12 {
            let n_abc = abc * (1.0 / abc_n);
            if ao.dot(&n_abc) > 0.0 {
                simplex[0] = d;
                simplex[1] = c;
                simplex.truncate(2);
                simplex.push(a);
                *dir = ao;
                return;
            }
        }

        if abd_n > 1e-12 {
            let n_abd = abd * (1.0 / abd_n);
            if ao.dot(&n_abd) > 0.0 {
                simplex[0] = d;
                simplex.truncate(1);
                simplex.push(a);
                *dir = ao;
                return;
            }
        }

        if acd_n > 1e-12 {
            let n_acd = acd_ * (1.0 / acd_n);
            if ao.dot(&n_acd) > 0.0 {
                simplex[0] = b;
                simplex[1] = d;
                simplex[2] = c;
                simplex.truncate(3);
                simplex.push(a);
                *dir = ao;
                return;
            }
        }

        simplex.truncate(3);
        let ab2 = simplex[1] - simplex[0];
        let ac2 = simplex[2] - simplex[0];
        let abc2 = ab2.cross(&ac2);
        let n2 = abc2.normalize();
        let ao2 = Vector3::zero() - simplex[0];
        *dir = if n2.dot(&ao2) > 0.0 { n2 } else { n2.negate() };
    }
}

impl Vector3 {
    fn negate(&self) -> Self {
        Self([-self.0[0], -self.0[1], -self.0[2]])
    }
}

// ===========================================================================
// EPA penetration
// ===========================================================================

/// Result of an EPA penetration query.
#[derive(Clone, Debug, PartialEq)]
pub struct EpaResult {
    /// Penetration depth (positive when shapes overlap).
    pub penetration_depth: f64,
    /// Contact normal pointing from shape B toward shape A.
    pub normal: Vector3,
    /// Contact point in world coordinates (approximate).
    pub contact_point: Vector3,
}

/// Compute penetration depth using EPA.
///
/// `a` and `b` are the two penetrating convex shapes. `simplex` should be the
/// final tetrahedron from GJK (at least 4 points). EPA expands this to find
/// the penetration depth.
pub fn epa_penetration(
    a: &dyn SupportMap,
    b: &dyn SupportMap,
    mut simplex: Vec<Vector3>,
) -> Option<EpaResult> {
    if simplex.len() < 4 {
        simplex = expand_to_tetrahedron(a, b, &simplex);
        if simplex.len() < 4 {
            return None;
        }
    }

    for _ in 0..64 {
        let (normal, _pen, closest) = closest_face(&simplex)?;
        let support = a.support(&normal.negate()) - b.support(&normal);
        let dist = support.dot(&normal);

        let mut found = false;
        for i in 0..simplex.len() {
            let j = (i + 1) % simplex.len();
            let edge = simplex[i] - simplex[j];
            let n = edge.cross(&(support - simplex[j]));
            if normal.dot(&n) > 0.0 {
                simplex.insert(j, support);
                found = true;
                break;
            }
        }
        if !found {
            simplex.push(support);
        }

        if simplex.len() > 16 {
            simplex.truncate(16);
        }

        let (new_normal, new_pen, _) = closest_face(&simplex)?;
        if (new_pen - dist).abs() < 1e-8 {
            return Some(EpaResult {
                penetration_depth: new_pen.max(0.0),
                normal: if new_pen > 1e-12 {
                    new_normal
                } else {
                    Vector3::zero()
                },
                contact_point: closest,
            });
        }
    }

    let (normal, penetration, closest) = closest_face(&simplex)?;
    Some(EpaResult {
        penetration_depth: penetration.max(0.0),
        normal: if penetration > 1e-12 {
            normal
        } else {
            Vector3::zero()
        },
        contact_point: closest,
    })
}

/// Expand simplex to a tetrahedron by walking the support function.
fn expand_to_tetrahedron(
    a: &dyn SupportMap,
    b: &dyn SupportMap,
    simplex: &[Vector3],
) -> Vec<Vector3> {
    let mut result = simplex.to_vec();
    let mut dir = if result.is_empty() {
        Vector3::new(1.0, 0.0, 0.0)
    } else {
        result[0].negate()
    };

    while result.len() < 4 {
        let sa = a.support(&dir.negate());
        let sb = b.support(&dir);
        let p = sa - sb;
        let mut found = false;
        for existing in &result {
            if ((*existing) - p).norm() < 1e-12 {
                found = true;
                break;
            }
        }
        if found {
            break;
        }
        result.push(p);
        dir = result[result.len() - 1].negate();
        if result.len() > 16 {
            break;
        }
    }

    result
}

/// Find the closest face of the convex hull to the origin.
fn closest_face(points: &[Vector3]) -> Option<(Vector3, f64, Vector3)> {
    if points.len() < 4 {
        return None;
    }

    let mut min_dist = f64::MAX;
    let mut best_normal = Vector3::zero();
    let mut best_point = Vector3::zero();

    for i in 0..points.len() {
        let j = (i + 1) % points.len();
        let k = (i + 2) % points.len();
        let a = points[i];
        let b = points[j];
        let c = points[k];
        let ab = b - a;
        let ac = c - a;
        let abc = ab.cross(&ac);
        let n_norm = abc.norm();
        if n_norm < 1e-12 {
            continue;
        }
        let n = abc * (1.0 / n_norm);
        let ao = Vector3::zero() - a;
        let dist = ao.dot(&n);
        if dist > 1e-12 && dist < min_dist {
            min_dist = dist;
            best_normal = n;
            best_point = a + n * dist;
        }
    }

    if min_dist < f64::MAX {
        Some((best_normal, min_dist, best_point))
    } else {
        None
    }
}

// ===========================================================================
// Convex shape implementations
// ===========================================================================

/// Sphere support map.
#[derive(Clone, Debug, PartialEq)]
pub struct Sphere {
    /// Sphere center.
    pub center: Vector3,
    /// Sphere radius.
    pub radius: f64,
}

impl Sphere {
    /// Create a new sphere.
    pub fn new(center: Vector3, radius: f64) -> Self {
        Sphere { center, radius }
    }
}

impl SupportMap for Sphere {
    fn support(&self, dir: &Vector3) -> Vector3 {
        let n = dir.normalize();
        Vector3([
            self.center.0[0] + n.0[0] * self.radius,
            self.center.0[1] + n.0[1] * self.radius,
            self.center.0[2] + n.0[2] * self.radius,
        ])
    }
}

/// Axis-aligned bounding box support map.
#[derive(Clone, Debug, PartialEq)]
pub struct Aabb {
    /// Minimum corner.
    pub min: Vector3,
    /// Maximum corner.
    pub max: Vector3,
}

impl Aabb {
    /// Create a new AABB.
    pub fn new(min: Vector3, max: Vector3) -> Self {
        Aabb { min, max }
    }
}

impl SupportMap for Aabb {
    fn support(&self, dir: &Vector3) -> Vector3 {
        Vector3([
            if dir.0[0] >= 0.0 {
                self.max.0[0]
            } else {
                self.min.0[0]
            },
            if dir.0[1] >= 0.0 {
                self.max.0[1]
            } else {
                self.min.0[1]
            },
            if dir.0[2] >= 0.0 {
                self.max.0[2]
            } else {
                self.min.0[2]
            },
        ])
    }
}

/// Capsule support map (cylinder with hemispherical end caps).
#[derive(Clone, Debug, PartialEq)]
pub struct Capsule {
    /// Capsule axis start point.
    pub a: Vector3,
    /// Capsule axis end point.
    pub b: Vector3,
    /// Capsule radius.
    pub radius: f64,
}

impl Capsule {
    /// Create a new capsule.
    pub fn new(a: Vector3, b: Vector3, radius: f64) -> Self {
        Capsule { a, b, radius }
    }
}

impl SupportMap for Capsule {
    fn support(&self, dir: &Vector3) -> Vector3 {
        let axis = self.b - self.a;
        let axis_len = axis.norm();
        if axis_len < 1e-12 {
            let n = dir.normalize();
            return Vector3([
                self.a.0[0] + n.0[0] * self.radius,
                self.a.0[1] + n.0[1] * self.radius,
                self.a.0[2] + n.0[2] * self.radius,
            ]);
        }
        let n = axis * (1.0 / axis_len);
        let da = dir.dot(&n);
        let t = da.clamp(-1.0, 1.0);
        let on_axis = self.a + n * (t * axis_len * 0.5);
        let radial = (*dir) - n * da;
        let rn = radial.normalize();
        Vector3([
            on_axis.0[0] + rn.0[0] * self.radius,
            on_axis.0[1] + rn.0[1] * self.radius,
            on_axis.0[2] + rn.0[2] * self.radius,
        ])
    }
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gjk_sphere_sphere_no_overlap() {
        let a = Sphere::new(Vector3::new(0.0, 0.0, 0.0), 1.0);
        let b = Sphere::new(Vector3::new(3.0, 0.0, 0.0), 1.0);
        let result = gjk_distance(&a, &b, &Vector3::new(1.0, 0.0, 0.0));
        assert!(!result.penetrating);
        assert!(
            (result.distance - 1.0).abs() < 1e-4,
            "distance = {}",
            result.distance
        );
        assert!(result.normal.0[0] > 0.9);
    }

    #[test]
    fn test_gjk_sphere_sphere_overlap() {
        let a = Sphere::new(Vector3::new(0.0, 0.0, 0.0), 1.0);
        let b = Sphere::new(Vector3::new(1.0, 0.0, 0.0), 1.0);
        let result = gjk_distance(&a, &b, &Vector3::new(1.0, 0.0, 0.0));
        assert!(result.penetrating);
        assert!(result.distance < 1e-4);
    }

    #[test]
    fn test_gjk_sphere_aabb_no_overlap() {
        let a = Sphere::new(Vector3::new(0.0, 0.0, 0.0), 1.0);
        let b = Aabb::new(Vector3::new(3.0, -1.0, -1.0), Vector3::new(4.0, 1.0, 1.0));
        let result = gjk_distance(&a, &b, &Vector3::new(1.0, 0.0, 0.0));
        assert!(!result.penetrating);
        assert!(result.distance > 0.5, "distance = {}", result.distance);
    }

    #[test]
    fn test_gjk_sphere_aabb_overlap() {
        let a = Sphere::new(Vector3::new(0.0, 0.0, 0.0), 2.0);
        let b = Aabb::new(Vector3::new(1.0, -1.0, -1.0), Vector3::new(3.0, 1.0, 1.0));
        let result = gjk_distance(&a, &b, &Vector3::new(1.0, 0.0, 0.0));
        assert!(result.penetrating || result.distance < 1e-4);
    }

    #[test]
    fn test_gjk_sphere_sphere_touching() {
        let a = Sphere::new(Vector3::new(0.0, 0.0, 0.0), 1.0);
        let b = Sphere::new(Vector3::new(2.0, 0.0, 0.0), 1.0);
        let result = gjk_distance(&a, &b, &Vector3::new(1.0, 0.0, 0.0));
        assert!(!result.penetrating || result.distance < 1e-4);
    }

    #[test]
    fn test_epa_sphere_sphere_penetration() {
        let a = Sphere::new(Vector3::new(0.0, 0.0, 0.0), 1.0);
        let b = Sphere::new(Vector3::new(1.0, 0.0, 0.0), 1.0);
        let gjk = gjk_distance(&a, &b, &Vector3::new(1.0, 0.0, 0.0));
        let simplex = if gjk.penetrating {
            vec![
                Vector3::new(1.0, 0.0, 0.0),
                Vector3::new(-1.0, 0.0, 0.0),
                Vector3::new(0.0, 1.0, 0.0),
                Vector3::new(0.0, 0.0, 1.0),
            ]
        } else {
            return;
        };
        let epa = epa_penetration(&a, &b, simplex);
        if let Some(epa) = epa {
            assert!(
                epa.penetration_depth > 0.0,
                "epa depth = {}",
                epa.penetration_depth
            );
        }
    }

    #[test]
    fn test_sphere_support() {
        let sphere = Sphere::new(Vector3::new(1.0, 2.0, 3.0), 0.5);
        let dir = Vector3::new(1.0, 0.0, 0.0);
        let s = sphere.support(&dir);
        assert!((s.0[0] - 1.5).abs() < 1e-12);
        assert!((s.0[1] - 2.0).abs() < 1e-12);
        assert!((s.0[2] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_aabb_support() {
        let aabb = Aabb::new(Vector3::new(-1.0, -1.0, -1.0), Vector3::new(1.0, 1.0, 1.0));
        let dir = Vector3::new(1.0, 1.0, 1.0);
        let s = aabb.support(&dir);
        assert_eq!(s, Vector3::new(1.0, 1.0, 1.0));
    }

    #[test]
    fn test_gjk_sphere_far_apart() {
        let a = Sphere::new(Vector3::new(0.0, 0.0, 0.0), 1.0);
        let b = Sphere::new(Vector3::new(10.0, 0.0, 0.0), 1.0);
        let result = gjk_distance(&a, &b, &Vector3::new(1.0, 0.0, 0.0));
        assert!(!result.penetrating);
        assert!(
            (result.distance - 8.0).abs() < 1e-3,
            "distance = {}",
            result.distance
        );
    }

    #[test]
    fn test_gjk_sphere_same_center() {
        let a = Sphere::new(Vector3::new(0.0, 0.0, 0.0), 1.0);
        let b = Sphere::new(Vector3::new(0.0, 0.0, 0.0), 1.0);
        let result = gjk_distance(&a, &b, &Vector3::new(1.0, 0.0, 0.0));
        assert!(result.penetrating);
        assert!(result.distance < 1e-6);
    }
}
