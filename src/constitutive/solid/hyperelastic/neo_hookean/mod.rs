#[cfg(test)]
mod test;

use super::*;

/// The Neo-Hookean hyperelastic constitutive model.
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
#[derive(Debug)]
pub struct NeoHookean<'a>
{
    parameters: Parameters<'a>
}

/// Constitutive model implementation of the Neo-Hookean hyperelastic constitutive model.
impl<'a> Constitutive<'a> for NeoHookean<'a>
{
    fn new(parameters: Parameters<'a>) -> Self
    {
        Self
        {
            parameters
        }
    }
}

/// Solid constitutive model implementation of the Neo-Hookean hyperelastic constitutive model.
impl<'a> Solid<'a> for NeoHookean<'a>
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

/// Elastic constitutive model implementation of the Neo-Hookean hyperelastic constitutive model.
impl<'a> Elastic<'a> for NeoHookean<'a>
{
    /// Calculates and returns the Cauchy stress.
    ///
    /// ```math
    /// \boldsymbol{\sigma}(\mathbf{F}) = \frac{\mu}{J}\,{\mathbf{B}^*}' + \frac{\kappa}{2}\left(J - \frac{1}{J}\right)\mathbf{1}
    /// ```
    fn calculate_cauchy_stress(&self, deformation_gradient: &DeformationGradient) -> CauchyStress
    {
        let jacobian = deformation_gradient.determinant();
        self.calculate_left_cauchy_green_deformation(deformation_gradient).deviatoric()/jacobian.powf(FIVE_THIRDS)*self.get_shear_modulus() + IDENTITY*self.get_bulk_modulus()*0.5*(jacobian - 1.0/jacobian)
    }
    /// Calculates and returns the tangent stiffness associated with the Cauchy stress.
    ///
    /// ```math
    /// \mathcal{T}_{ijkL}(\mathbf{F}) = \frac{\mu}{J^{5/3}}\left(\delta_{ik}F_{jL} + \delta_{jk}F_{iL} - \frac{2}{3}\,\delta_{ij}F_{kL} - \frac{5}{3} \, B_{ij}'F_{kL}^{-T} \right) + \frac{\kappa}{2} \left(J + \frac{1}{J}\right)\delta_{ij}F_{kL}^{-T}
    /// ```
    fn calculate_cauchy_tangent_stiffness(&self, deformation_gradient: &DeformationGradient) -> CauchyTangentStiffness
    {
        let (inverse_transpose_deformation_gradient, jacobian) = deformation_gradient.inverse_transpose_and_determinant();
        let scaled_shear_modulus = self.get_shear_modulus()/jacobian.powf(FIVE_THIRDS);
        (CauchyTangentStiffness::dyad_ik_jl(&IDENTITY, deformation_gradient) + CauchyTangentStiffness::dyad_il_jk(deformation_gradient, &IDENTITY) - CauchyTangentStiffness::dyad_ij_kl(&IDENTITY, deformation_gradient)*(TWO_THIRDS))*scaled_shear_modulus + CauchyTangentStiffness::dyad_ij_kl(&(IDENTITY*(0.5*self.get_bulk_modulus()*(jacobian + 1.0/jacobian)) - self.calculate_left_cauchy_green_deformation(deformation_gradient).deviatoric()*(scaled_shear_modulus*FIVE_THIRDS)), &inverse_transpose_deformation_gradient)
    }
}

/// Hyperelastic constitutive model implementation of the Neo-Hookean hyperelastic constitutive model.
impl<'a> Hyperelastic<'a> for NeoHookean<'a>
{
    /// Calculates and returns the Helmholtz free energy density.
    ///
    /// ```math
    /// a(\mathbf{F}) = \frac{\mu}{2}\left[\mathrm{tr}(\mathbf{B}^*) - 3\right] + \frac{\kappa}{2}\left[\frac{1}{2}\left(J^2 - 1\right) - \ln J\right]
    /// ```
    fn calculate_helmholtz_free_energy_density(&self, deformation_gradient: &DeformationGradient) -> Result<Scalar, ConstitutiveError>
    {
        let jacobian = deformation_gradient.determinant();
        if jacobian > 0.0 {
            Ok(0.5*(self.get_shear_modulus()*(self.calculate_left_cauchy_green_deformation(deformation_gradient).trace()/jacobian.powf(TWO_THIRDS) - 3.0) + self.get_bulk_modulus()*(0.5*(jacobian.powi(2) - 1.0) - jacobian.ln())))
        } else {
            Err(ConstitutiveError::InvalidJacobianElastic(jacobian, deformation_gradient.copy(), format!("{:?}", &self)))
        }
    }
}
