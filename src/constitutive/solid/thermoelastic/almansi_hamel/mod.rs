#[cfg(test)]
mod test;

use super::*;

/// The Almansi-Hamel thermoelastic constitutive model.
///
/// **Parameters**
/// - The bulk modulus $`\kappa`$.
/// - The shear modulus $`\mu`$.
/// - The coefficient of thermal expansion $`\alpha`$.
/// - The reference temperature $`T_\mathrm{ref}`$.
///
/// **External variables**
/// - The deformation gradient $`\mathbf{F}`$.
/// - The temperature $`T`$.
///
/// **Internal variables**
/// - None.
///
/// **Notes**
/// - The Almansi-Hamel strain measure is given by $`\mathbf{e}=\tfrac{1}{2}(\mathbf{1} - \mathbf{B}^{-1})`$.
pub struct AlmansiHamel<'a>
{
    parameters: Parameters<'a>
}

/// Constitutive model implementation of the Almansi-Hamel thermoelastic constitutive model.
impl<'a> Constitutive<'a> for AlmansiHamel<'a>
{
    fn new(parameters: Parameters<'a>) -> Self
    {
        Self
        {
            parameters
        }
    }
}

/// Solid constitutive model implementation of the Almansi-Hamel thermoelastic constitutive model.
impl<'a> Solid<'a> for AlmansiHamel<'a> {}

/// Thermoelastic constitutive model implementation of the Almansi-Hamel thermoelastic constitutive model.
impl<'a> Thermoelastic<'a> for AlmansiHamel<'a>
{
    /// Calculates and returns the Cauchy stress.
    ///
    /// ```math
    /// \boldsymbol{\sigma}(\mathbf{F}, T) = \frac{2\mu}{J}\,\mathbf{e}' + \frac{\kappa}{J}\Big[\mathrm{tr}(\mathbf{e}) - 3\alpha(T - T_\mathrm{ref})\Big]\mathbf{1}
    /// ```
    fn calculate_cauchy_stress(&self, deformation_gradient: &DeformationGradient, temperature: &Scalar) -> CauchyStress
    {
        let identity = LeftCauchyGreenDeformation::identity();
        let (inverse_deformation_gradient, jacobian) = deformation_gradient.inverse_and_determinant();
        let strain = (&identity * 1.0 - inverse_deformation_gradient.transpose() * &inverse_deformation_gradient) * 0.5;
        let (deviatoric_strain, strain_trace) = strain.deviatoric_and_trace();
        deviatoric_strain * (2.0 * self.get_shear_modulus() / jacobian) + identity * (self.get_bulk_modulus() / jacobian * (strain_trace - 3.0*self.get_coefficient_of_thermal_expansion()*(temperature - self.get_reference_temperature())))
    }
    /// Calculates and returns the tangent stiffness associated with the Cauchy stress.
    ///
    /// ```math
    /// \mathcal{T}_{ijkL}(\mathbf{F}, T) = \frac{\mu}{J}\left[B_{jk}^{-1}F_{iL}^{-T} + B_{ik}^{-1}F_{jL}^{-T} - \frac{2}{3}\,\delta_{ij}B_{km}^{-1}F_{mL}^{-T} - 2e_{ij}'F_{kL}^{-T}\right] + \frac{\kappa}{J}\left\{\delta_{ij}B_{km}^{-1}F_{mL}^{-T} - \Big[\mathrm{tr}(\mathbf{e}) - 3\alpha(T - T_\mathrm{ref})\Big]\delta_{ij}F_{kL}^{-T}\right\}
    /// ```
    fn calculate_cauchy_tangent_stiffness(&self, deformation_gradient: &DeformationGradient, temperature: &Scalar) -> CauchyTangentStiffness
    {
        let identity = LeftCauchyGreenDeformation::identity();
        let (inverse_transpose_deformation_gradient, jacobian) = deformation_gradient.inverse_transpose_and_determinant();
        let inverse_left_cauchy_green_deformation = &inverse_transpose_deformation_gradient * inverse_transpose_deformation_gradient.transpose();
        let strain = (&identity * 1.0 - &inverse_left_cauchy_green_deformation) * 0.5;
        let (deviatoric_strain, strain_trace) = strain.deviatoric_and_trace();
        (CauchyTangentStiffness::dyad_il_jk(&inverse_transpose_deformation_gradient, &inverse_left_cauchy_green_deformation) + CauchyTangentStiffness::dyad_ik_jl(&inverse_left_cauchy_green_deformation, &inverse_transpose_deformation_gradient)) * (self.get_shear_modulus() / jacobian)+ CauchyTangentStiffness::dyad_ij_kl(&identity,&(inverse_left_cauchy_green_deformation * &inverse_transpose_deformation_gradient*((self.get_bulk_modulus() - self.get_shear_modulus() * 2.0 / 3.0) / jacobian))) - CauchyTangentStiffness::dyad_ij_kl(&(deviatoric_strain * (2.0 * self.get_shear_modulus() / jacobian) + identity * (self.get_bulk_modulus() / jacobian * (strain_trace - 3.0*self.get_coefficient_of_thermal_expansion()*(temperature - self.get_reference_temperature())))), &inverse_transpose_deformation_gradient)
    }
    fn get_bulk_modulus(&self) -> &Scalar
    {
        &self.parameters[0]
    }
    fn get_shear_modulus(&self) -> &Scalar
    {
        &self.parameters[1]
    }
    fn get_coefficient_of_thermal_expansion(&self) -> &Scalar
    {
        &self.parameters[2]
    }
    fn get_reference_temperature(&self) -> &Scalar
    {
        &self.parameters[3]
    }
}
