#[cfg(test)]
mod test;

use super::*;

/// The Fung hyperelastic constitutive model.
///
/// **Parameters**
/// - The bulk modulus $`\kappa`$.
/// - The shear modulus $`\mu`$.
/// - The extra modulus $`\mu_m`$.
/// - The exponent $`a`$.
///
/// **External variables**
/// - The deformation gradient $`\mathbf{F}`$.
///
/// **Internal variables**
/// - None.
///
/// **Notes**
/// - The Fung model reduces to the [Neo-Hookean model](NeoHookean) when $`\mu_m\to 0`$ or $`a\to 0`$.
pub struct Fung<'a>
{
    parameters: Parameters<'a>
}

/// Inherent implementation of the Fung hyperelastic constitutive model.
impl<'a> Fung<'a>
{
    /// Returns the extra modulus.
    pub fn get_extra_modulus(&self) -> &Scalar
    {
        &self.parameters[2]
    }
    /// Returns the exponent.
    fn get_exponent(&self) -> &Scalar
    {
        &self.parameters[3]
    }
}

/// Constitutive model implementation of the Fung hyperelastic constitutive model.
impl<'a> Constitutive<'a> for Fung<'a>
{
    fn new(parameters: Parameters<'a>) -> Self
    {
        Self
        {
            parameters
        }
    }
}

/// Solid constitutive model implementation of the Fung hyperelastic constitutive model.
impl<'a> Solid<'a> for Fung<'a>
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

/// Elastic constitutive model implementation of the Fung hyperelastic constitutive model.
impl<'a> Elastic<'a> for Fung<'a>
{
    /// Calculates and returns the Cauchy stress.
    ///
    /// ```math
    /// \boldsymbol{\sigma}(\mathbf{F}) = \frac{1}{J}\left[\mu + \mu_m\left(e^{a[\mathrm{tr}(\mathbf{B}^* ) - 3]} - 1\right)\right]{\mathbf{B}^* }' + \frac{\kappa}{2}\left(J - \frac{1}{J}\right)\mathbf{1}
    /// ```
    fn calculate_cauchy_stress(&self, deformation_gradient: &DeformationGradient) -> CauchyStress
    {
        let jacobian = deformation_gradient.determinant();
        let isochoric_left_cauchy_green_deformation = self.calculate_left_cauchy_green_deformation(deformation_gradient)/jacobian.powf(2.0/3.0);
        let (deviatoric_isochoric_left_cauchy_green_deformation, isochoric_left_cauchy_green_deformation_trace) = isochoric_left_cauchy_green_deformation.deviatoric_and_trace();
        deviatoric_isochoric_left_cauchy_green_deformation*((self.get_shear_modulus() + self.get_extra_modulus()*((self.get_exponent()*(isochoric_left_cauchy_green_deformation_trace - 3.0)).exp() - 1.0))/jacobian) + LeftCauchyGreenDeformation::identity()*self.get_bulk_modulus()*0.5*(jacobian - 1.0/jacobian)
    }
    /// Calculates and returns the tangent stiffness associated with the Cauchy stress.
    ///
    /// ```math
    /// \mathcal{T}_{ijkL}(\mathbf{F}) = \frac{1}{J^{5/3}}\left[\mu + \mu_m\left(e^{a[\mathrm{tr}(\mathbf{B}^* ) - 3]} - 1\right)\right]\left(\delta_{ik}F_{jL} + \delta_{jk}F_{iL} - \frac{2}{3}\,\delta_{ij}F_{kL} - \frac{5}{3} \, B_{ij}'F_{kL}^{-T} \right) + \frac{2a\mu_m}{J^{7/3}}\,e^{a[\mathrm{tr}(\mathbf{B}^* ) - 3]}B_{ij}'B_{km}'F_{mL}^{-T} + \frac{\kappa}{2} \left(J + \frac{1}{J}\right)\delta_{ij}F_{kL}^{-T}
    /// ```
    fn calculate_cauchy_tangent_stiffness(&self, deformation_gradient: &DeformationGradient) -> CauchyTangentStiffness
    {
        let identity = CauchyStress::identity();
        let (inverse_transpose_deformation_gradient, jacobian) = deformation_gradient.inverse_transpose_and_determinant();
        let isochoric_left_cauchy_green_deformation = self.calculate_left_cauchy_green_deformation(deformation_gradient)/jacobian.powf(2.0/3.0);
        let (deviatoric_isochoric_left_cauchy_green_deformation, isochoric_left_cauchy_green_deformation_trace) = isochoric_left_cauchy_green_deformation.deviatoric_and_trace();
        let exponential = (self.get_exponent()*(isochoric_left_cauchy_green_deformation_trace - 3.0)).exp();
        let scaled_shear_modulus_0 = (self.get_shear_modulus() + self.get_extra_modulus()*(exponential - 1.0))/jacobian.powf(5.0/3.0);
        (CauchyTangentStiffness::dyad_ik_jl(&identity, deformation_gradient) + CauchyTangentStiffness::dyad_il_jk(deformation_gradient, &identity) - CauchyTangentStiffness::dyad_ij_kl(&identity, deformation_gradient)*(2.0/3.0))*scaled_shear_modulus_0 + CauchyTangentStiffness::dyad_ij_kl(&deviatoric_isochoric_left_cauchy_green_deformation, &((&deviatoric_isochoric_left_cauchy_green_deformation*&inverse_transpose_deformation_gradient)*(2.0*self.get_exponent()*self.get_extra_modulus()*exponential/jacobian))) + CauchyTangentStiffness::dyad_ij_kl(&(identity*(0.5*self.get_bulk_modulus()*(jacobian + 1.0/jacobian)) - self.calculate_left_cauchy_green_deformation(deformation_gradient).deviatoric()*(scaled_shear_modulus_0*5.0/3.0)), &inverse_transpose_deformation_gradient)
    }
}

/// Hyperelastic constitutive model implementation of the Fung hyperelastic constitutive model.
impl<'a> Hyperelastic<'a> for Fung<'a>
{
    /// Calculates and returns the Helmholtz free energy density.
    ///
    /// ```math
    /// a(\mathbf{F}) = \frac{\mu - \mu_m}{2}\left[\mathrm{tr}(\mathbf{B}^* ) - 3\right] + \frac{\mu_m}{2a}\left(e^{a[\mathrm{tr}(\mathbf{B}^* ) - 3]} - 1\right)
    /// ```
    fn calculate_helmholtz_free_energy_density(&self, deformation_gradient: &DeformationGradient) -> Scalar
    {
        let jacobian = deformation_gradient.determinant();
        let scalar_term = self.calculate_left_cauchy_green_deformation(deformation_gradient).trace()/jacobian.powf(2.0/3.0) - 3.0;
        0.5*((self.get_shear_modulus() - self.get_extra_modulus())*scalar_term + self.get_extra_modulus()/self.get_exponent()*((self.get_exponent()*scalar_term).exp() - 1.0) + self.get_bulk_modulus()*(0.5*(jacobian.powi(2) - 1.0) - jacobian.ln()))
    }
}
