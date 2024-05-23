#[cfg(test)]
mod test;

use super::*;

/// The Gent hyperelastic constitutive model.
/// 
/// **Parameters**
/// - The bulk modulus $`\kappa`$.
/// - The shear modulus $`\mu`$.
/// - The extensibility $`J_m`$.
///
/// **External variables**
/// - The deformation gradient $`\mathbf{F}`$.
///
/// **Internal variables**
/// - None.
/// 
/// **Notes**
/// - The Gent model reduces to the [Neo-Hookean model](NeoHookean) when $`J_m\to\infty`$.
pub struct Gent<'a>
{
    parameters: Parameters<'a>
}

/// Inherent implementation of the Gent hyperelastic constitutive model.
impl<'a> Gent<'a>
{
    /// Returns the extensibility.
    fn get_extensibility(&self) -> &Scalar
    {
        &self.parameters[2]
    }
}

/// Constitutive model implementation of the Gent hyperelastic constitutive model.
impl<'a> Constitutive<'a> for Gent<'a>
{
    fn new(parameters: Parameters<'a>) -> Self
    {
        Self
        {
            parameters
        }
    }
}

/// Solid constitutive model implementation of the Gent hyperelastic constitutive model.
impl<'a> Solid<'a> for Gent<'a>
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

/// Elastic constitutive model implementation of the Gent hyperelastic constitutive model.
impl<'a> Elastic<'a> for Gent<'a>
{
    /// Calculates and returns the Cauchy stress.
    ///
    /// ```math
    /// \boldsymbol{\sigma}(\mathbf{F}) = \frac{J^{-1}\mu J_m {\mathbf{B}^* }'}{J_m - \mathrm{tr}(\mathbf{B}^* ) + 3} + \frac{\kappa}{2}\left(J - \frac{1}{J}\right)\mathbf{1}
    /// ```
    fn calculate_cauchy_stress(&self, deformation_gradient: &DeformationGradient) -> CauchyStress
    {
        let jacobian = deformation_gradient.determinant();
        let isochoric_left_cauchy_green_deformation = self.calculate_left_cauchy_green_deformation(deformation_gradient)/jacobian.powf(2.0/3.0);
        let (deviatoric_isochoric_left_cauchy_green_deformation, isochoric_left_cauchy_green_deformation_trace) = isochoric_left_cauchy_green_deformation.deviatoric_and_trace();
        let denominator = self.get_extensibility() - isochoric_left_cauchy_green_deformation_trace + 3.0;
        if denominator < 0.0
        {
            panic!()
        }
        (deviatoric_isochoric_left_cauchy_green_deformation*self.get_shear_modulus()*self.get_extensibility()/jacobian)/denominator + LeftCauchyGreenDeformation::identity()*self.get_bulk_modulus()*0.5*(jacobian - 1.0/jacobian)
    }
    /// Calculates and returns the tangent stiffness associated with the Cauchy stress.
    ///
    /// ```math
    /// \mathcal{T}_{ijkL}(\mathbf{F}) = \frac{J^{-5/3}\mu J_m}{J_m - \mathrm{tr}(\mathbf{B}^* ) + 3}\Bigg[ \delta_{ik}F_{jL} + \delta_{jk}F_{iL} - \frac{2}{3}\,\delta_{ij}F_{kL} + \frac{2{B_{ij}^* }' F_{kL}}{J_m - \mathrm{tr}(\mathbf{B}^* ) + 3} - \left(\frac{5}{3} + \frac{2}{3}\frac{\mathrm{tr}(\mathbf{B}^* )}{J_m - \mathrm{tr}(\mathbf{B}^* ) + 3}\right) J^{2/3} {B_{ij}^* }' F_{kL}^{-T} \Bigg] + \frac{\kappa}{2} \left(J + \frac{1}{J}\right)\delta_{ij}F_{kL}^{-T}
    /// ```
    fn calculate_cauchy_tangent_stiffness(&self, deformation_gradient: &DeformationGradient) -> CauchyTangentStiffness
    {
        let identity = CauchyStress::identity();
        let (inverse_transpose_deformation_gradient, jacobian) = deformation_gradient.inverse_transpose_and_determinant();
        let isochoric_left_cauchy_green_deformation = self.calculate_left_cauchy_green_deformation(deformation_gradient)/jacobian.powf(2.0/3.0);
        let (deviatoric_isochoric_left_cauchy_green_deformation, isochoric_left_cauchy_green_deformation_trace) = isochoric_left_cauchy_green_deformation.deviatoric_and_trace();
        let denominator = self.get_extensibility() - isochoric_left_cauchy_green_deformation_trace + 3.0;
        if denominator < 0.0
        {
            panic!()
        }
        let prefactor = self.get_shear_modulus()*self.get_extensibility()/jacobian/denominator;
        (CauchyTangentStiffness::dyad_ik_jl(&identity, deformation_gradient) + CauchyTangentStiffness::dyad_il_jk(deformation_gradient, &identity) - CauchyTangentStiffness::dyad_ij_kl(&identity, deformation_gradient)*(2.0/3.0) + CauchyTangentStiffness::dyad_ij_kl(&deviatoric_isochoric_left_cauchy_green_deformation, deformation_gradient)*(2.0/denominator))*(prefactor/jacobian.powf(2.0/3.0)) + CauchyTangentStiffness::dyad_ij_kl(&(identity*(0.5*self.get_bulk_modulus()*(jacobian + 1.0/jacobian)) - deviatoric_isochoric_left_cauchy_green_deformation*prefactor*((5.0 + 2.0*isochoric_left_cauchy_green_deformation_trace/denominator)/3.0)), &inverse_transpose_deformation_gradient)
    }
}

/// Hyperelastic constitutive model implementation of the Gent hyperelastic constitutive model.
impl<'a> Hyperelastic<'a> for Gent<'a>
{
    /// Calculates and returns the Helmholtz free energy density.
    ///
    /// ```math
    /// a(\mathbf{F}) = -\frac{\mu J_m}{2}\,\ln\left[1 - \frac{\mathrm{tr}(\mathbf{B}^* ) - 3}{J_m}\right] + \frac{\kappa}{2}\left[\frac{1}{2}\left(J^2 - 1\right) - \ln J\right]
    /// ```
    fn calculate_helmholtz_free_energy_density(&self, deformation_gradient: &DeformationGradient) -> Scalar
    {
        let jacobian = deformation_gradient.determinant();
        let factor = (self.calculate_left_cauchy_green_deformation(deformation_gradient).trace()/jacobian.powf(2.0/3.0) - 3.0)/self.get_extensibility();
        if factor > 1.0
        {
            panic!()
        }
        0.5*(-self.get_shear_modulus()*self.get_extensibility()*(1.0 - factor).ln() + self.get_bulk_modulus()*(0.5*(jacobian.powi(2) - 1.0) - jacobian.ln()))
    }
}
