#[cfg(test)]
mod test;

use super::*;

pub struct YeohModel<'a>
{
    parameters: ConstitutiveModelParameters<'a>
}

impl<'a> YeohModel<'a>
{
    fn get_moduli(&self) -> &[Scalar]
    {
        &self.parameters[1..]
    }
    fn get_extra_moduli(&self) -> &[Scalar]
    {
        &self.parameters[2..]
    }
}

impl<'a> ConstitutiveModel<'a> for YeohModel<'a>
{
    fn calculate_helmholtz_free_energy_density(&self, deformation_gradient: &DeformationGradient) -> Scalar
    {
        let jacobian = deformation_gradient.determinant();
        let scalar_term = self.calculate_left_cauchy_green_deformation(deformation_gradient).trace()/jacobian.powf(2.0/3.0) - 3.0;
        0.5*(self.get_moduli().iter().enumerate().map(|(n, modulus)| modulus*scalar_term.powi(<usize as TryInto<i32>>::try_into(n).unwrap() + 1)).sum::<Scalar>() + self.get_bulk_modulus()*(0.5*(jacobian.powi(2) - 1.0) - jacobian.ln()))
    }
    fn calculate_cauchy_stress(&self, deformation_gradient: &DeformationGradient) -> CauchyStress
    {
        let jacobian = deformation_gradient.determinant();
        let left_cauchy_green_deformation = self.calculate_left_cauchy_green_deformation(deformation_gradient);
        let scalar_term = left_cauchy_green_deformation.trace()/jacobian.powf(2.0/3.0) - 3.0;
        left_cauchy_green_deformation.deviatoric()*self.get_moduli().iter().enumerate().map(|(n, modulus)| ((n as Scalar) + 1.0)*modulus*scalar_term.powi(n.try_into().unwrap())).sum::<Scalar>()/jacobian.powf(5.0/3.0) + LeftCauchyGreenDeformation::identity()*self.get_bulk_modulus()*0.5*(jacobian - 1.0/jacobian)
    }
    fn calculate_cauchy_tangent_stiffness(&self, deformation_gradient: &DeformationGradient) -> CauchyTangentStiffness
    {
        let identity = CauchyStress::identity();
        let jacobian = deformation_gradient.determinant();
        let inverse_transpose_deformation_gradient = deformation_gradient.inverse_transpose();
        let left_cauchy_green_deformation = self.calculate_left_cauchy_green_deformation(deformation_gradient);
        let scalar_term = left_cauchy_green_deformation.trace()/jacobian.powf(2.0/3.0) - 3.0;
        let scaled_modulus = self.get_moduli().iter().enumerate().map(|(n, modulus)| ((n as Scalar) + 1.0)*modulus*scalar_term.powi(n.try_into().unwrap())).sum::<Scalar>()/jacobian.powf(5.0/3.0);
        let deviatoric_left_cauchy_green_deformation = left_cauchy_green_deformation.deviatoric();
        let last_term = CauchyTangentStiffness::dyad_ij_kl(&deviatoric_left_cauchy_green_deformation, &((left_cauchy_green_deformation.deviatoric() * &inverse_transpose_deformation_gradient) * (2.0*self.get_extra_moduli().iter().enumerate().map(|(n, modulus)| ((n as Scalar) + 2.0)*((n as Scalar) + 1.0)*modulus*scalar_term.powi(n.try_into().unwrap())).sum::<Scalar>()/jacobian.powf(7.0/3.0))));
        (CauchyTangentStiffness::dyad_ik_jl(&identity, deformation_gradient) + CauchyTangentStiffness::dyad_il_jk(deformation_gradient, &identity) - CauchyTangentStiffness::dyad_ij_kl(&identity, deformation_gradient)*(2.0/3.0))*scaled_modulus + CauchyTangentStiffness::dyad_ij_kl(&(identity*(0.5*self.get_bulk_modulus()*(jacobian + 1.0/jacobian)) - deviatoric_left_cauchy_green_deformation*(scaled_modulus*5.0/3.0)), &inverse_transpose_deformation_gradient) + last_term
    }
    fn new(parameters: ConstitutiveModelParameters<'a>) -> Self
    {
        Self
        {
            parameters
        }
    }
}

impl<'a> HyperelasticConstitutiveModel<'a> for YeohModel<'a>
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