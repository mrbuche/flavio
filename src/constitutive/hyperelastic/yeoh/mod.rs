#[cfg(test)]
mod test;

use super::*;

/// The Yeoh hyperelastic constitutive model.
/// 
/// **Parameters**
/// - The bulk modulus $`\kappa`$.
/// - The shear modulus $`\mu`$.
/// - The extra moduli $`\mu_n`$ for $`n=2\ldots N`$.
///
/// **External variables**
/// - The deformation gradient $`\mathbf{F}`$.
///
/// **Internal variables**
/// - None.
/// 
/// **Notes**
/// - The Yeoh model reduces to the [Neo-Hookean model](NeoHookeanModel) when $`\mu_n\to 0`$ for $`n=2\ldots N`$.
pub struct YeohModel<'a>
{
    parameters: ConstitutiveModelParameters<'a>
}

/// Base implementation of the Yeoh hyperelastic constitutive model.
impl<'a> YeohModel<'a>
{
    /// Returns an array of the moduli.
    pub fn get_moduli(&self) -> &[Scalar]
    {
        &self.parameters[1..]
    }
    /// Returns an array of the extra moduli.
    pub fn get_extra_moduli(&self) -> &[Scalar]
    {
        &self.parameters[2..]
    }
}

/// Constitutive model implementation of the Yeoh hyperelastic constitutive model.
impl<'a> ConstitutiveModel<'a> for YeohModel<'a>
{
    /// Calculates and returns the Cauchy stress.
    ///
    /// ```math
    /// \boldsymbol{\sigma}(\mathbf{F}) = \sum_{n=1}^N \frac{n\mu_n}{J}\left[\mathrm{tr}(\mathbf{B}^* ) - 3\right]^{n-1}\,{\mathbf{B}^*}' + \frac{\kappa}{2}\left(J - \frac{1}{J}\right)\mathbf{1}
    /// ```
    fn calculate_cauchy_stress(&self, deformation_gradient: &DeformationGradient) -> CauchyStress
    {
        let jacobian = deformation_gradient.determinant();
        let (deviatoric_left_cauchy_green_deformation, left_cauchy_green_deformation_trace) = self.calculate_left_cauchy_green_deformation(deformation_gradient).deviatoric_and_trace();
        let scalar_term = left_cauchy_green_deformation_trace/jacobian.powf(2.0/3.0) - 3.0;
        deviatoric_left_cauchy_green_deformation*self.get_moduli().iter().enumerate().map(|(n, modulus)| ((n as Scalar) + 1.0)*modulus*scalar_term.powi(n.try_into().unwrap())).sum::<Scalar>()/jacobian.powf(5.0/3.0) + LeftCauchyGreenDeformation::identity()*self.get_bulk_modulus()*0.5*(jacobian - 1.0/jacobian)
    }
    /// Calculates and returns the tangent stiffness associated with the Cauchy stress.
    ///
    /// ```math
    /// \mathcal{T}_{ijkL}(\mathbf{F}) = \sum_{n=1}^N \frac{n\mu_n}{J^{5/3}}\left[\mathrm{tr}(\mathbf{B}^* ) - 3\right]^{n-1}\left(\delta_{ik}F_{jL} + \delta_{jk}F_{iL} - \frac{2}{3}\,\delta_{ij}F_{kL}- \frac{5}{3} \, B_{ij}'F_{kL}^{-T} \right) + \sum_{n=2}^N \frac{2n(n-1)\mu_n}{J^{7/3}}\left[\mathrm{tr}(\mathbf{B}^* ) - 3\right]^{n-2}B_{ij}'B_{km}'F_{mL}^{-T} + \frac{\kappa}{2} \left(J + \frac{1}{J}\right)\delta_{ij}F_{kL}^{-T}
    /// ```
    fn calculate_cauchy_tangent_stiffness(&self, deformation_gradient: &DeformationGradient) -> CauchyTangentStiffness
    {
        let identity = CauchyStress::identity();
        let (inverse_transpose_deformation_gradient, jacobian) = deformation_gradient.inverse_transpose_and_determinant();
        let left_cauchy_green_deformation = self.calculate_left_cauchy_green_deformation(deformation_gradient);
        let scalar_term = left_cauchy_green_deformation.trace()/jacobian.powf(2.0/3.0) - 3.0;
        let scaled_modulus = self.get_moduli().iter().enumerate().map(|(n, modulus)| ((n as Scalar) + 1.0)*modulus*scalar_term.powi(n.try_into().unwrap())).sum::<Scalar>()/jacobian.powf(5.0/3.0);
        let deviatoric_left_cauchy_green_deformation = left_cauchy_green_deformation.deviatoric();
        let last_term = CauchyTangentStiffness::dyad_ij_kl(&deviatoric_left_cauchy_green_deformation, &((left_cauchy_green_deformation.deviatoric() * &inverse_transpose_deformation_gradient) * (2.0*self.get_extra_moduli().iter().enumerate().map(|(n, modulus)| ((n as Scalar) + 2.0)*((n as Scalar) + 1.0)*modulus*scalar_term.powi(n.try_into().unwrap())).sum::<Scalar>()/jacobian.powf(7.0/3.0))));
        (CauchyTangentStiffness::dyad_ik_jl(&identity, deformation_gradient) + CauchyTangentStiffness::dyad_il_jk(deformation_gradient, &identity) - CauchyTangentStiffness::dyad_ij_kl(&identity, deformation_gradient)*(2.0/3.0))*scaled_modulus + CauchyTangentStiffness::dyad_ij_kl(&(identity*(0.5*self.get_bulk_modulus()*(jacobian + 1.0/jacobian)) - deviatoric_left_cauchy_green_deformation*(scaled_modulus*5.0/3.0)), &inverse_transpose_deformation_gradient) + last_term
    }
    /// Calculates and returns the Helmholtz free energy density.
    ///
    /// ```math
    /// a(\mathbf{F}) = \sum_{n=1}^N \frac{\mu_n}{2}\left[\mathrm{tr}(\mathbf{B}^* ) - 3\right]^n + \frac{\kappa}{2}\left[\frac{1}{2}\left(J^2 - 1\right) - \ln J\right]
    /// ```
    fn calculate_helmholtz_free_energy_density(&self, deformation_gradient: &DeformationGradient) -> Scalar
    {
        let jacobian = deformation_gradient.determinant();
        let scalar_term = self.calculate_left_cauchy_green_deformation(deformation_gradient).trace()/jacobian.powf(2.0/3.0) - 3.0;
        0.5*(self.get_moduli().iter().enumerate().map(|(n, modulus)| modulus*scalar_term.powi(<usize as TryInto<i32>>::try_into(n).unwrap() + 1)).sum::<Scalar>() + self.get_bulk_modulus()*(0.5*(jacobian.powi(2) - 1.0) - jacobian.ln()))
    }
    fn new(parameters: ConstitutiveModelParameters<'a>) -> Self
    {
        Self
        {
            parameters
        }
    }
}

/// Hyperelastic constitutive model implementation of the Yeoh hyperelastic constitutive model.
impl<'a> HyperelasticConstitutiveModel for YeohModel<'a>
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