#[cfg(test)]
mod test;

use super::*;

/// The Saint Venant-Kirchoff hyperelastoplastic constitutive model.
///
/// **Parameters**
/// - The bulk modulus $`\kappa`$.
/// - The shear modulus $`\mu`$.
/// - The reference plastic flow rate $`\dot{\gamma}^\mathrm{p}_0`$.
/// - The shear strength $`S_0`$.
///
/// **External variables**
/// - The deformation gradient $`\mathbf{F}`$.
///
/// **Internal variables**
/// - The plastic deformation gradient $`\mathbf{F}^\mathrm{p}`$.
///
/// **Notes**
/// - The Green-Saint Venant strain measure is given by $`\mathbf{E}=\tfrac{1}{2}(\mathbf{C}-\mathbf{1})`$.
pub struct SaintVenantKirchoff<'a>
{
    parameters: Parameters<'a>,
    plastic_deformation_gradient: DeformationGradientPlastic
}

/// Constitutive model implementation of the Saint Venant-Kirchoff hyperelastoplastic constitutive model.
impl<'a> Constitutive<'a> for SaintVenantKirchoff<'a>
{
    fn new(parameters: Parameters<'a>) -> Self
    {
        Self
        {
            parameters,
            plastic_deformation_gradient: DeformationGradientPlastic::identity()
        }
    }
}

/// Solid constitutive model implementation of the Saint Venant-Kirchoff hyperelastoplastic constitutive model.
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

/// Elastic constitutive model implementation of the Saint Venant-Kirchoff hyperelastic constitutive model.
impl<'a> Elastic<'a> for SaintVenantKirchoff<'a>
{
    /// Calculates and returns the second Piola-Kirchoff stress.
    ///
    /// ```math
    /// \mathbf{S}(\mathbf{F}^\mathrm{e}) = 2\mu{\mathbf{E}^\mathrm{e}}' + \kappa\,\mathrm{tr}(\mathbf{E}^\mathrm{e})\mathbf{1}
    /// ```
    fn calculate_second_piola_kirchoff_stress(&self, deformation_gradient: &DeformationGradient) -> SecondPiolaKirchoffStress
    {
        let elastic_deformation = self.compute_elastic_deformation(deformation_gradient).into();
        let (deviatoric_strain, strain_trace) = ((self.calculate_right_cauchy_green_deformation(&elastic_deformation) - RightCauchyGreenDeformation::identity())*0.5).deviatoric_and_trace();
        deviatoric_strain*(2.0*self.get_shear_modulus()) + RightCauchyGreenDeformation::identity()*(self.get_bulk_modulus()*strain_trace)
    }
    /// Calculates and returns the tangent stiffness associated with the second Piola-Kirchoff stress.
    ///
    /// ```math
    /// \mathcal{G}_{IJkL}(\mathbf{F}^\mathrm{e}) = \mu\,\delta_{JL}F^\mathrm{e}_{kI} + \mu\,\delta_{IL}F^\mathrm{e}_{kJ} + \left(\kappa - \frac{2}{3}\,\mu\right)\delta_{IJ}F^\mathrm{e}_{kL}
    /// ```
    fn calculate_second_piola_kirchoff_tangent_stiffness(&self, deformation_gradient: &DeformationGradient) -> SecondPiolaKirchoffTangentStiffness
    {
        let elastic_deformation: DeformationGradient = self.compute_elastic_deformation(deformation_gradient).into();
        let identity = SecondPiolaKirchoffStress::identity();
        let scaled_deformation_gradient_transpose = elastic_deformation.transpose()*self.get_shear_modulus();
        SecondPiolaKirchoffTangentStiffness::dyad_ik_jl(&scaled_deformation_gradient_transpose, &identity) + SecondPiolaKirchoffTangentStiffness::dyad_il_jk(&identity, &scaled_deformation_gradient_transpose) + SecondPiolaKirchoffTangentStiffness::dyad_ij_kl(&(identity*(self.get_bulk_modulus() - 2.0/3.0*self.get_shear_modulus())), &elastic_deformation)
    }
}

/// Hyperelastic constitutive model implementation of the Saint Venant-Kirchoff hyperelastic constitutive model.
impl<'a> Hyperelastic<'a> for SaintVenantKirchoff<'a>
{
    /// Calculates and returns the Helmholtz free energy density.
    ///
    /// ```math
    /// a(\mathbf{F}^\mathrm{e}) = \mu\,\mathrm{tr}({\mathbf{E}^\mathrm{e}}^2) + \frac{1}{2}\left(\kappa - \frac{2}{3}\,\mu\right)\mathrm{tr}(\mathbf{E}^\mathrm{e})^2
    /// ```
    fn calculate_helmholtz_free_energy_density(&self, deformation_gradient: &DeformationGradient) -> Scalar
    {
        let elastic_deformation = self.compute_elastic_deformation(deformation_gradient).into();
        let strain = (self.calculate_right_cauchy_green_deformation(&elastic_deformation) - RightCauchyGreenDeformation::identity())*0.5;
        self.get_shear_modulus()*strain.squared_trace() + 0.5*(self.get_bulk_modulus() - 2.0/3.0*self.get_shear_modulus())*strain.trace().powi(2)
    }
}

/// Plastic constitutive model implementation of the Saint Venant-Kirchoff hyperelastic constitutive model.
impl<'a> Plastic<'a> for SaintVenantKirchoff<'a>
{
    /// Calculates and returns the plastic flow rate.
    ///
    /// ```math
    /// \dot{\gamma}^\mathrm{p}(\mathbf{M}^\mathrm{e}) = \dot{\gamma}^\mathrm{p}_0 \sinh\left(\frac{|{\mathbf{M}^\mathrm{e}}'|}{\sqrt{2}\,S_0}\right)
    /// ```
    fn compute_plastic_flow_rate(&self, deviatoric_mandel_stress_norm: &Scalar) -> Scalar
    {
        self.get_reference_plastic_flow_rate() * (deviatoric_mandel_stress_norm / self.get_shear_strength() / 2.0_f64.sqrt()).sinh()
    }
    fn get_reference_plastic_flow_rate(&self) -> &Scalar
    {
        &self.parameters[2]
    }
    fn get_shear_strength(&self) -> &Scalar
    {
        &self.parameters[3]
    }
}

/// Elastoplastic constitutive model implementation of the Saint Venant-Kirchoff hyperelastic constitutive model.
impl<'a> Elastoplastic<'a> for SaintVenantKirchoff<'a>
{
    fn get_plastic_deformation_gradient(&self) -> &DeformationGradientPlastic
    {
        &self.plastic_deformation_gradient
    }
    fn set_plastic_deformation_gradient(&mut self, plastic_deformation_gradient: DeformationGradientPlastic)
    {
        self.plastic_deformation_gradient = plastic_deformation_gradient
    }
}

/// Hyperelastoplastic constitutive model implementation of the Saint Venant-Kirchoff hyperelastic constitutive model.
impl<'a> Hyperelastoplastic<'a> for SaintVenantKirchoff<'a>
{}
