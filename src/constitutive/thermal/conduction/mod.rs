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
    ThermalConstitutiveModel
};

pub use fourier::FourierModel;

/// Required methods for thermal conduction constitutive models.
pub trait ThermalConductionConstitutiveModel<'a>
where
    Self: ConstitutiveModel<'a> + ThermalConstitutiveModel<'a>
{
    /// Calculates and returns the heat flux.
    fn calculate_heat_flux(&self, temperature_gradient: &TemperatureGradient) -> HeatFlux;
}