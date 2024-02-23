#[cfg(test)]
mod test;

use super::*;

/// The Almansi-Hamel elastic constitutive model.
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
/// - The Almansi-Hamel strain measure is given by $`\mathbf{e}=\tfrac{1}{2}(\mathbf{1}-\mathbf{B}^{-1})`$.
pub struct AlmansiHamel<'a>
{
    parameters: Parameters<'a>
}

/// Constitutive model implementation of the Almansi-Hamel elastic constitutive model.
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

/// Solid constitutive model implementation of the Almansi-Hamel elastic constitutive model.
impl<'a> Solid<'a> for AlmansiHamel<'a>
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

/// Elastic constitutive model implementation of the Almansi-Hamel elastic constitutive model.
impl<'a> Elastic<'a> for AlmansiHamel<'a>
{
    /// Calculates and returns the Cauchy stress.
    ///
    /// ```math
    /// \boldsymbol{\sigma}(\mathbf{F}) = 2\mu\mathbf{e}' + \kappa\,\mathrm{tr}(\mathbf{e})\mathbf{1}
    /// ```
    fn calculate_cauchy_stress(&self, deformation_gradient: &DeformationGradient) -> CauchyStress
    {
        let identity = LeftCauchyGreenDeformation::identity();
        let inverse_deformation_gradient = deformation_gradient.inverse();
        let strain = (&identity * 1.0 - inverse_deformation_gradient.transpose() * &inverse_deformation_gradient) * 0.5;
        let (deviatoric_strain, strain_trace) = strain.deviatoric_and_trace();
        deviatoric_strain * (2.0 * self.get_shear_modulus()) + identity * (self.get_bulk_modulus() * strain_trace)
    }
    /// Calculates and returns the tangent stiffness associated with the Cauchy stress.
    ///
    /// ```math
    /// \mathcal{T}_{ijkL}(\mathbf{F}) = \mu B_{jk}^{-1}F_{iL}^{-T} + \mu B_{ik}^{-1}F_{jL}^{-T} + \left(\kappa - \frac{2}{3}\,\mu\right)\delta_{ij}B_{km}^{-1}F_{mL}^{-T}
    /// ```
    fn calculate_cauchy_tangent_stiffness(&self, deformation_gradient: &DeformationGradient) -> CauchyTangentStiffness
    {
        let identity = LeftCauchyGreenDeformation::identity();
        let inverse_transpose_deformation_gradient = deformation_gradient.inverse_transpose();
        let inverse_left_cauchy_green_deformation = &inverse_transpose_deformation_gradient * inverse_transpose_deformation_gradient.transpose();
        (CauchyTangentStiffness::dyad_il_jk(&inverse_transpose_deformation_gradient, &inverse_left_cauchy_green_deformation) + CauchyTangentStiffness::dyad_ik_jl(&inverse_left_cauchy_green_deformation, &inverse_transpose_deformation_gradient)) * self.get_shear_modulus() + CauchyTangentStiffness::dyad_ij_kl(&identity, &(inverse_left_cauchy_green_deformation * &inverse_transpose_deformation_gradient*(self.get_bulk_modulus() - self.get_shear_modulus() * 2.0 / 3.0)))
    }
}
