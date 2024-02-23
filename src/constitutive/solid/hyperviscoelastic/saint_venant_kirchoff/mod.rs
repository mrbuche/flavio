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
/// - The deformation gradient rate $`\dot{\mathbf{F}}`$.
///
/// **Internal variables**
/// - None.
///
/// **Notes**
/// - The Green-Saint Venant strain measure is given by $`\mathbf{E}=\tfrac{1}{2}(\mathbf{C}-\mathbf{1})`$.
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
    fn calculate_second_piola_kirchoff_stress(&self, deformation_gradient: &DeformationGradient, deformation_gradient_rate: &DeformationGradientRate) -> SecondPiolaKirchoffStress
    {
        let (deviatoric_strain, strain_trace) = ((self.calculate_right_cauchy_green_deformation(deformation_gradient) - RightCauchyGreenDeformation::identity())*0.5).deviatoric_and_trace();
        let first_term = deformation_gradient_rate.transpose()*deformation_gradient;
        let (deviatoric_strain_rate, strain_rate_trace) = ((&first_term + first_term.transpose())*0.5).deviatoric_and_trace();
        deviatoric_strain*(2.0*self.get_shear_modulus()) + deviatoric_strain_rate*(2.0*self.get_shear_viscosity()) + RightCauchyGreenDeformation::identity()*(self.get_bulk_modulus()*strain_trace + self.get_bulk_viscosity()*strain_rate_trace)
    }
    /// Calculates and returns the rate tangent stiffness associated with the second Piola-Kirchoff stress.
    ///
    /// ```math
    /// \mathcal{W}_{IJkL}(\mathbf{F}) = \eta\,\delta_{JL}F_{kI} + \eta\,\delta_{IL}F_{kJ} + \left(\zeta - \frac{2}{3}\,\eta\right)\delta_{IJ}F_{kL}
    /// ```
    fn calculate_second_piola_kirchoff_rate_tangent_stiffness(&self, deformation_gradient: &DeformationGradient, _: &DeformationGradientRate) -> SecondPiolaKirchoffRateTangentStiffness
    {
        let identity = SecondPiolaKirchoffStress::identity();
        let scaled_deformation_gradient_transpose = deformation_gradient.transpose()*self.get_shear_viscosity();
        SecondPiolaKirchoffRateTangentStiffness::dyad_ik_jl(&scaled_deformation_gradient_transpose, &identity) + SecondPiolaKirchoffRateTangentStiffness::dyad_il_jk(&identity, &scaled_deformation_gradient_transpose) + SecondPiolaKirchoffRateTangentStiffness::dyad_ij_kl(&(identity*(self.get_bulk_viscosity() - 2.0/3.0*self.get_shear_viscosity())), deformation_gradient)
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
    /// Calculates and returns the viscous dissipation.
    ///
    /// ```math
    /// \phi(\mathbf{F},\dot{\mathbf{F}}) = \eta\,\mathrm{tr}(\dot{\mathbf{E}}^2) + \left(\zeta - \frac{2}{3}\,\eta\right)\mathrm{tr}(\dot{\mathbf{E}})^2
    /// ```
    fn calculate_viscous_dissipation(&self, deformation_gradient: &DeformationGradient, deformation_gradient_rate: &DeformationGradientRate) -> Scalar
    {
        let first_term = deformation_gradient_rate.transpose()*deformation_gradient;
        let strain_rate = (&first_term + first_term.transpose())*0.5;
        self.get_shear_viscosity()*strain_rate.squared_trace() + 0.5*(self.get_bulk_viscosity() - 2.0/3.0*self.get_shear_viscosity())*strain_rate.trace().powi(2)
    }
}
