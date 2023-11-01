#[cfg(test)]
mod test;

use super::*;

pub struct MooneyRivlinModel<'a>
{
    parameters: ConstitutiveModelParameters<'a>
}

impl<'a> MooneyRivlinModel<'a>
{
    fn get_extra_modulus(&self) -> &Scalar
    {
        &self.parameters[2]
    }
}

impl<'a> ConstitutiveModel<'a> for MooneyRivlinModel<'a>
{
    fn calculate_cauchy_stress(&self, deformation_gradient: &DeformationGradient) -> CauchyStress
    {
        let jacobian = deformation_gradient.determinant();
        let isochoric_left_cauchy_green_deformation = self.calculate_left_cauchy_green_deformation(deformation_gradient)/jacobian.powf(2.0/3.0);
        ((isochoric_left_cauchy_green_deformation.deviatoric()*(self.get_shear_modulus() - self.get_extra_modulus()) - isochoric_left_cauchy_green_deformation.inverse().deviatoric()*self.get_extra_modulus()) + CauchyStress::identity()*(self.get_bulk_modulus()*0.5*(jacobian.powi(2) - 1.0)))/jacobian
    }
    fn calculate_cauchy_tangent_stiffness(&self, deformation_gradient: &DeformationGradient) -> CauchyTangentStiffness
    {
        let identity = CauchyStress::identity();
        let (inverse_transpose_deformation_gradient, jacobian) = deformation_gradient.inverse_transpose_and_determinant();
        let scaled_delta_shear_modulus = (self.get_shear_modulus() - self.get_extra_modulus())/jacobian.powf(5.0/3.0);
        let inverse_isochoric_left_cauchy_green_deformation = (self.calculate_left_cauchy_green_deformation(deformation_gradient)/jacobian.powf(2.0/3.0)).inverse();
        let deviatoric_inverse_isochoric_left_cauchy_green_deformation = inverse_isochoric_left_cauchy_green_deformation.deviatoric();
        let term_1 = CauchyTangentStiffness::dyad_ij_kl(&inverse_isochoric_left_cauchy_green_deformation, &inverse_transpose_deformation_gradient)*2.0/3.0 - CauchyTangentStiffness::dyad_ik_jl(&inverse_isochoric_left_cauchy_green_deformation, &inverse_transpose_deformation_gradient) - CauchyTangentStiffness::dyad_il_jk(&inverse_transpose_deformation_gradient, &inverse_isochoric_left_cauchy_green_deformation);
        let term_3 = CauchyTangentStiffness::dyad_ij_kl(&deviatoric_inverse_isochoric_left_cauchy_green_deformation, &inverse_transpose_deformation_gradient);
        let term_2 = CauchyTangentStiffness::dyad_ij_kl(&identity, &((deviatoric_inverse_isochoric_left_cauchy_green_deformation * 2.0/3.0) * &inverse_transpose_deformation_gradient));
        (CauchyTangentStiffness::dyad_ik_jl(&identity, deformation_gradient) + CauchyTangentStiffness::dyad_il_jk(deformation_gradient, &identity) - CauchyTangentStiffness::dyad_ij_kl(&identity, deformation_gradient)*(2.0/3.0))*scaled_delta_shear_modulus + CauchyTangentStiffness::dyad_ij_kl(&(identity*(0.5*self.get_bulk_modulus()*(jacobian + 1.0/jacobian)) - self.calculate_left_cauchy_green_deformation(deformation_gradient).deviatoric()*(scaled_delta_shear_modulus*5.0/3.0)), &inverse_transpose_deformation_gradient) - (term_1 + term_2 - term_3)*self.get_extra_modulus()/jacobian
    }
    fn calculate_helmholtz_free_energy_density(&self, deformation_gradient: &DeformationGradient) -> Scalar
    {
        let jacobian = deformation_gradient.determinant();
        let isochoric_left_cauchy_green_deformation = self.calculate_left_cauchy_green_deformation(deformation_gradient)/jacobian.powf(2.0/3.0);
        0.5*((self.get_shear_modulus() - self.get_extra_modulus())*(isochoric_left_cauchy_green_deformation.trace() - 3.0) + self.get_extra_modulus()*(isochoric_left_cauchy_green_deformation.second_invariant() - 3.0) + self.get_bulk_modulus()*(0.5*(jacobian.powi(2) - 1.0) - jacobian.ln()))
    }
    fn new(parameters: ConstitutiveModelParameters<'a>) -> Self
    {
        Self
        {
            parameters
        }
    }
}

impl<'a> HyperelasticConstitutiveModel for MooneyRivlinModel<'a>
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
