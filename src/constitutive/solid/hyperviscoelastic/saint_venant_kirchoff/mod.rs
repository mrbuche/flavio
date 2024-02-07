#[cfg(test)]
mod test;

use super::*;

/// The Saint Venant-Kirchoff hyperviscoelastic constitutive model.
///
/// **Parameters**
/// - The bulk modulus $`\kappa`$.
/// - The shear modulus $`\mu`$.
/// - The bulk viscosity $`\zeta`$.
/// - The shear viscosity $`\eta`$.
///
/// **External variables**
/// - The deformation gradient $`\mathbf{F}`$.
/// - The rate of deformation gradient $`\dot{\mathbf{F}}`$.
///
/// **Internal variables**
/// - None.
///
/// **Notes**
/// - The Green-Saint Venant strain measure is given by $`\mathbf{E}=\tfrac{1}{2}(\mathbf{C} - \mathbf{1})`$.
pub struct SaintVenantKirchoff<'a>
{
    parameters: Parameters<'a>
}

/// Constitutive model implementation of the Saint Venant-Kirchoff hyperviscoelastic constitutive model.
impl<'a> Constitutive<'a> for SaintVenantKirchoff<'a>
{
    fn new(parameters: Parameters<'a>) -> Self
    {
        Self
        {
            parameters
        }
    }
}

/// Solid constitutive model implementation of the Saint Venant-Kirchoff hyperviscoelastic constitutive model.
impl<'a> Solid<'a> for SaintVenantKirchoff<'a>
{
    fn get_bulk_modulus(&self) -> &Scalar
    {
        &self.parameters[0]
    }
    fn get_shear_modulus(&self) -> &Scalar
    {
        &self.parameters[1]
    }
}

/// Viscous constitutive model implementation of the Saint Venant-Kirchoff hyperviscoelastic constitutive model.
impl<'a> Viscous<'a> for SaintVenantKirchoff<'a>
{
    fn get_bulk_viscosity(&self) -> &Scalar
    {
        &self.parameters[2]
    }
    fn get_shear_viscosity(&self) -> &Scalar
    {
        &self.parameters[3]
    }
}

/// Viscoelastic constitutive model implementation of the Saint Venant-Kirchoff hyperviscoelastic constitutive model.
impl<'a> Viscoelastic<'a> for SaintVenantKirchoff<'a>
{
    /// Calculates and returns the second Piola-Kirchoff stress.
    ///
    /// ```math
    /// \mathbf{S}(\mathbf{F},\dot\mathbf{F}) = 2\mu\mathbf{E}' + \kappa\,\mathrm{tr}(\mathbf{E})\mathbf{1} + 2\eta\dot{\mathbf{E}}' + \zeta\,\mathrm{tr}(\dot{\mathbf{E}})\mathbf{1}
    /// ```
    fn calculate_second_piola_kirchoff_stress(&self, deformation_gradient: &DeformationGradient, deformation_gradient_dot: &DeformationGradientDot) -> SecondPiolaKirchoffStress
    {
        let (deviatoric_strain, strain_trace) = ((self.calculate_right_cauchy_green_deformation(deformation_gradient) - RightCauchyGreenDeformation::identity())*0.5).deviatoric_and_trace();
        let (deviatoric_strain_rate, strain_rate_trace) = ((self.calculate_right_cauchy_green_deformation(deformation_gradient_dot) - RightCauchyGreenDeformation::identity())*0.5).deviatoric_and_trace();
        deviatoric_strain*(2.0*self.get_shear_modulus()) + RightCauchyGreenDeformation::identity()*(self.get_bulk_modulus()*strain_trace) + deviatoric_strain_rate*(2.0*self.get_shear_viscosity()) + RightCauchyGreenDeformation::identity()*(self.get_bulk_viscosity()*strain_rate_trace)
    }
}

/// Hyperviscoelastic constitutive model implementation of the Saint Venant-Kirchoff hyperviscoelastic constitutive model.
impl<'a> Hyperviscoelastic<'a> for SaintVenantKirchoff<'a>
{
    /// Calculates and returns the Helmholtz free energy density.
    ///
    /// ```math
    /// a(\mathbf{F}) = \mu\,\mathrm{tr}(\mathbf{E}^2) + \left(\kappa - \frac{2}{3}\,\mu\right)\mathrm{tr}(\mathbf{E})^2
    /// ```
    fn calculate_helmholtz_free_energy_density(&self, deformation_gradient: &DeformationGradient) -> Scalar
    {
        let strain = (self.calculate_right_cauchy_green_deformation(deformation_gradient) - RightCauchyGreenDeformation::identity())*0.5;
        self.get_shear_modulus()*strain.squared_trace() + 0.5*(self.get_bulk_modulus() - 2.0/3.0*self.get_shear_modulus())*strain.trace().powi(2)
    }
    /// Calculates and returns the Helmholtz free energy density.
    ///
    /// ```math
    /// b(\dot{\mathbf{F}}) = \mu\,\mathrm{tr}(\dot{\mathbf{E}}^2) + \left(\kappa - \frac{2}{3}\,\mu\right)\mathrm{tr}(\dot{\mathbf{E}})^2
    /// ```
    fn calculate_viscous_dissipation(&self, deformation_gradient_dot: &DeformationGradientDot) -> Scalar
    {
        let strain_rate = (self.calculate_right_cauchy_green_deformation(deformation_gradient_dot) - RightCauchyGreenDeformation::identity())*0.5;
        self.get_shear_modulus()*strain_rate.squared_trace() + 0.5*(self.get_bulk_modulus() - 2.0/3.0*self.get_shear_modulus())*strain_rate.trace().powi(2)
    }
}
