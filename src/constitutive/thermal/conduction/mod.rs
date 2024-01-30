//! Thermal conduction constitutive models.

#[cfg(test)]
mod test;

mod fourier;

use super::*;

pub use self::fourier::FourierModel;

/// Required methods for thermal conduction constitutive models.
pub trait ThermalConduction<'a>
where
    Self: ConstitutiveModel<'a>
{
    /// Calculates and returns the heat flux.
    fn calculate_heat_flux(&self, temperature_gradient: &TemperatureGradient) -> HeatFlux;
}