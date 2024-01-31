#[cfg(test)]
mod test;

use super::*;

/// The Fourier thermal conduction constitutive model.
///
/// **Parameters**
/// - The thermal conductivity $`k`$.
///
/// **External variables**
/// - The temperature gradient $`\nabla T`$.
///
/// **Internal variables**
/// - None.
pub struct FourierModel<'a>
{
    parameters: ConstitutiveModelParameters<'a>
}

/// Inherent implementation of the Fourier thermal conduction constitutive model.
impl<'a> FourierModel<'a>
{
    fn get_thermal_conductivity(&self) -> &Scalar
    {
        &self.parameters[0]
    }
}

/// Constitutive model implementation of the Fourier thermal conduction constitutive model.
impl<'a> ConstitutiveModel<'a> for FourierModel<'a>
{
    fn new(parameters: ConstitutiveModelParameters<'a>) -> Self
    {
        Self
        {
            parameters
        }
    }
}
/// Thermal constitutive model implementation of the Fourier thermal conduction constitutive model.
impl<'a> ThermalConstitutiveModel<'a> for FourierModel<'a> {}

/// Thermal conduction constitutive model implementation of the Fourier thermal conduction constitutive model.
impl<'a> ThermalConduction<'a> for FourierModel<'a>
{
    /// Calculates and returns the heat flux.
    ///
    /// ```math
    /// \mathbf{q}(\nabla T) = -k\nabla T
    /// ```
    fn calculate_heat_flux(&self, temperature_gradient: &TemperatureGradient) -> HeatFlux
    {
        temperature_gradient * -self.get_thermal_conductivity()
    }
}