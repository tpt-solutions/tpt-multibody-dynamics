//! Archard wear law: volume loss proportional to contact pressure times sliding distance.

/// Parameters for the Archard wear model.
#[derive(Clone, Debug, PartialEq)]
pub struct WearParams {
    /// Dimensionless wear coefficient.
    pub wear_coefficient: f64,
    /// Hardness of the softer material (pressure units).
    pub hardness: f64,
}

impl Default for WearParams {
    fn default() -> Self {
        Self {
            wear_coefficient: 1e-6,
            hardness: 1.0e9,
        }
    }
}

/// Archard wear law: volume_loss = k * (normal_force * sliding_distance) / hardness.
pub struct ArchardWear {
    /// Wear model parameters.
    pub params: WearParams,
}

impl ArchardWear {
    /// Create a new Archard wear model.
    pub fn new(params: WearParams) -> Self {
        Self { params }
    }

    /// Compute accumulated volume loss for a given normal force and sliding distance.
    pub fn volume_loss(&self, normal_force: f64, sliding_distance: f64) -> f64 {
        if self.params.hardness <= 0.0 {
            return 0.0;
        }
        self.params.wear_coefficient * normal_force * sliding_distance / self.params.hardness
    }
}

/// Computes the wear rate given normal pressure and sliding velocity.
///
/// Returns volume loss per unit time.
pub fn compute_wear_rate(normal_pressure: f64, sliding_velocity: f64, params: &WearParams) -> f64 {
    if params.hardness <= 0.0 || sliding_velocity <= 0.0 {
        return 0.0;
    }
    params.wear_coefficient * normal_pressure * sliding_velocity / params.hardness
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_archard_wear() {
        let aw = ArchardWear::new(WearParams {
            wear_coefficient: 1e-6,
            hardness: 1.0e9,
        });
        let normal_force = 100.0;
        let sliding_distance = 0.1;
        let loss = aw.volume_loss(normal_force, sliding_distance);
        let expected = 1e-6 * normal_force * sliding_distance / 1.0e9;
        assert!((loss - expected).abs() < 1e-20);
    }

    #[test]
    fn test_wear_rate() {
        let rate = compute_wear_rate(
            1e6,
            0.01,
            &WearParams {
                wear_coefficient: 1e-6,
                hardness: 1.0e9,
            },
        );
        let expected = 1e-6 * 1e6 * 0.01 / 1.0e9;
        assert!((rate - expected).abs() < 1e-20);
    }

    #[test]
    fn test_wear_zero_hardness() {
        let rate = compute_wear_rate(1e6, 0.01, &WearParams {
            wear_coefficient: 1e-6,
            hardness: 0.0,
        });
        assert_eq!(rate, 0.0);
    }

    #[test]
    fn test_wear_zero_velocity() {
        let rate = compute_wear_rate(1e6, 0.0, &WearParams {
            wear_coefficient: 1e-6,
            hardness: 1.0e9,
        });
        assert_eq!(rate, 0.0);
    }

    #[test]
    fn test_archard_wear_zero_force() {
        let aw = ArchardWear::new(WearParams::default());
        assert_eq!(aw.volume_loss(0.0, 0.1), 0.0);
    }

    #[test]
    fn test_wear_proportional_to_force() {
        let aw = ArchardWear::new(WearParams {
            wear_coefficient: 1e-6,
            hardness: 1.0e9,
        });
        let loss1 = aw.volume_loss(100.0, 0.1);
        let loss2 = aw.volume_loss(200.0, 0.1);
        assert!((loss2 - 2.0 * loss1).abs() < 1e-20);
    }
}
