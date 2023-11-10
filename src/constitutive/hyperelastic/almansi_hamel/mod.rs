#[cfg(test)]
mod test;

use super::*;

/// The Almansi-Hamel hyperelastic constitutive model.
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
/// - The Almansi-Hamel strain measure is given by $`\tfrac{1}{2}(\mathbf{1} - \mathbf{B}^{-1})`$.
pub struct AlmansiHamelModel<'a>
{
    parameters: ConstitutiveModelParameters<'a>
}

/// Constitutive model implementation of the Almansi-Hamel hyperelastic constitutive model.
impl<'a> ConstitutiveModel<'a> for AlmansiHamelModel<'a>
{
    /// Calculates and returns the Cauchy stress.
    ///
    /// ```math
    /// \boldsymbol{\sigma}(\mathbf{F}) = \frac{\mu}{J}\left(\mathbf{1} - \mathbf{B}^{-1}\right)' + \frac{\kappa}{2J}\,\mathrm{tr}\left(\mathbf{1} - \mathbf{B}^{-1}\right)\mathbf{1}
    /// ```
    fn calculate_cauchy_stress(&self, deformation_gradient: &DeformationGradient) -> CauchyStress
    {
        let (inverse_deformation_gradient, jacobian) = deformation_gradient.inverse_and_determinant();
        let almansi_hamel = LeftCauchyGreenDeformation::identity() - inverse_deformation_gradient.transpose() * &inverse_deformation_gradient;
        let (deviatoric_almansi_hamel, almansi_hamel_trace) = almansi_hamel.deviatoric_and_trace();
        let spherical_term = LeftCauchyGreenDeformation::identity() * almansi_hamel_trace * self.get_bulk_modulus() / 2.0;
        (deviatoric_almansi_hamel * self.get_shear_modulus() + spherical_term) / jacobian
    }
    /// Calculates and returns the tangent stiffness associated with the Cauchy stress.
    ///
    /// ```math
    /// \mathcal{T}_{ijkL}(\mathbf{F}) = ?
    /// ```
    fn calculate_cauchy_tangent_stiffness(&self, deformation_gradient: &DeformationGradient) -> CauchyTangentStiffness
    {
        CauchyTangentStiffness::zero()
    }
    /// Calculates and returns the Helmholtz free energy density.
    ///
    /// ```math
    /// a(\mathbf{F}) = ?
    /// ```
    fn calculate_helmholtz_free_energy_density(&self, deformation_gradient: &DeformationGradient) -> Scalar
    {
        0.0
    }
    fn new(parameters: ConstitutiveModelParameters<'a>) -> Self
    {
        Self
        {
            parameters
        }
    }
}

/// Hyperelastic constitutive model implementation of the Almansi-Hamel hyperelastic constitutive model.
impl<'a> HyperelasticConstitutiveModel for AlmansiHamelModel<'a>
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
