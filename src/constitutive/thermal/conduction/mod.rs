//! Thermal conduction constitutive models.

#[cfg(test)]
pub mod test;

mod fourier;

use super::
{
    ConstitutiveModel,
    ConstitutiveModelParameters,
    HeatFlux,
    Scalar,
    TemperatureGradient,
    Thermal
};

pub use fourier::Fourier;

/// Required methods for thermal conduction constitutive models.
pub trait ThermalConduction<'a>
where
    Self: ConstitutiveModel<'a> + Thermal<'a>
{
    /// Calculates and returns the heat flux.
    fn calculate_heat_flux(&self, temperature_gradient: &TemperatureGradient) -> HeatFlux;
}