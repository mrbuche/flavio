#[cfg(test)]
mod test;

use super::*;

/// The Saint Venant-Kirchoff hyperelastic constitutive model.
/// 
/// **Parameters**
/// - The bulk modulus $`\kappa`$.
/// - The shear modulus $`\mu`$.
///
/// **External variables**
/// - The deformation gradient $`\mathbf{F}`$.
///
/// **Internal variables**
/// - None.
/// 
/// **Notes**
/// - The Green-Saint Venant strain measure is given by $`\mathbf{E}=\tfrac{1}{2}(\mathbf{C} - \mathbf{1})`$.
pub struct SaintVenantKirchoffModel<'a>
{
    parameters: ConstitutiveModelParameters<'a>
}

/// Constitutive model implementation of the Saint Venant-Kirchoff hyperelastic constitutive model.
impl<'a> ConstitutiveModel<'a> for SaintVenantKirchoffModel<'a>
{
    /// Calculates and returns the second Piola-Kirchoff stress.
    ///
    /// ```math
    /// \boldsymbol{\sigma}(\mathbf{F}) = 2\mu\mathbf{E}' + \kappa\,\mathrm{tr}(\mathbf{E})\mathbf{1}
    /// ```
    fn calculate_second_piola_kirchoff_stress(&self, deformation_gradient: &DeformationGradient) -> SecondPiolaKirchoffStress
    {
        let (deviatoric_strain, strain_trace) = ((self.calculate_right_cauchy_green_deformation(deformation_gradient) - RightCauchyGreenDeformation::identity())*0.5).deviatoric_and_trace();
        deviatoric_strain*(2.0*self.get_shear_modulus()) + RightCauchyGreenDeformation::identity()*(self.get_bulk_modulus()*strain_trace)
    }
    /// Calculates and returns the tangent stiffness associated with the second Piola-Kirchoff stress.
    ///
    /// ```math
    /// \mathcal{G}_{IJkL} = 2\mu\,\delta_{JL}F_{kI} + \left(\kappa - \frac{2}{3}\,\mu\right)\delta_{IJ}F_{kL}
    /// ```
    fn calculate_second_piola_kirchoff_tangent_stiffness(&self, deformation_gradient: &DeformationGradient) -> SecondPiolaKirchoffTangentStiffness
    {
        let identity = SecondPiolaKirchoffStress::identity();
        SecondPiolaKirchoffTangentStiffness::dyad_ik_jl(&(deformation_gradient.transpose()*(2.0*self.get_shear_modulus())), &identity) + SecondPiolaKirchoffTangentStiffness::dyad_ij_kl(&(identity*(self.get_bulk_modulus() - 2.0/3.0*self.get_shear_modulus())), &deformation_gradient)
    }
    fn new(parameters: ConstitutiveModelParameters<'a>) -> Self
    {
        Self
        {
            parameters
        }
    }
}

/// Hyperelastic constitutive model implementation of the Saint Venant-Kirchoff hyperelastic constitutive model.
impl<'a> HyperelasticConstitutiveModel for SaintVenantKirchoffModel<'a>
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
    fn get_bulk_modulus(&self) -> &Scalar
    {
        &self.parameters[0]
    }
    fn get_shear_modulus(&self) -> &Scalar
    {
        &self.parameters[1]
    }
}
