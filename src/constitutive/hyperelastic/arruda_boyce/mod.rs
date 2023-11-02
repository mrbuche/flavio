#[cfg(test)]
mod test;

use crate::math::inverse_langevin;
use super::*;

pub struct ArrudaBoyceModel<'a>
{
    parameters: ConstitutiveModelParameters<'a>
}

impl<'a> ArrudaBoyceModel<'a>
{
    fn get_number_of_links(&self) -> &Scalar
    {
        &self.parameters[2]
    }
}

impl<'a> ConstitutiveModel<'a> for ArrudaBoyceModel<'a>
{
    fn calculate_cauchy_stress(&self, deformation_gradient: &DeformationGradient) -> CauchyStress
    {
        let jacobian = deformation_gradient.determinant();
        let (deviatoric_isochoric_left_cauchy_green_deformation, isochoric_left_cauchy_green_deformation_trace) = (self.calculate_left_cauchy_green_deformation(deformation_gradient)/jacobian.powf(2.0/3.0)).deviatoric_and_trace();
        let gamma = (isochoric_left_cauchy_green_deformation_trace/3.0/self.get_number_of_links()).sqrt();
        deviatoric_isochoric_left_cauchy_green_deformation*(self.get_shear_modulus()*inverse_langevin(gamma)/gamma/3.0/jacobian) + LeftCauchyGreenDeformation::identity()*self.get_bulk_modulus()*0.5*(jacobian - 1.0/jacobian)
    }
    fn calculate_cauchy_tangent_stiffness(&self, deformation_gradient: &DeformationGradient) -> CauchyTangentStiffness
    {
        let identity = CauchyStress::identity();
        let (inverse_transpose_deformation_gradient, jacobian) = deformation_gradient.inverse_transpose_and_determinant();
        let left_cauchy_green_deformation = self.calculate_left_cauchy_green_deformation(deformation_gradient);
        let (deviatoric_isochoric_left_cauchy_green_deformation, isochoric_left_cauchy_green_deformation_trace) = (left_cauchy_green_deformation/jacobian.powf(2.0/3.0)).deviatoric_and_trace();
        let gamma = (isochoric_left_cauchy_green_deformation_trace/3.0/self.get_number_of_links()).sqrt();
        let eta = inverse_langevin(gamma);
        let scaled_shear_modulus = self.get_shear_modulus()*eta/gamma/3.0/jacobian.powf(5.0/3.0);
        (CauchyTangentStiffness::dyad_ik_jl(&identity, deformation_gradient) + CauchyTangentStiffness::dyad_il_jk(deformation_gradient, &identity) - CauchyTangentStiffness::dyad_ij_kl(&identity, deformation_gradient)*(2.0/3.0))*scaled_shear_modulus + CauchyTangentStiffness::dyad_ij_kl(&(identity*(0.5*self.get_bulk_modulus()*(jacobian + 1.0/jacobian)) - deviatoric_isochoric_left_cauchy_green_deformation*(scaled_shear_modulus*5.0/3.0)), &inverse_transpose_deformation_gradient)
        //
        // need term from d/dF (eta/gamma)
    }
    fn calculate_helmholtz_free_energy_density(&self, deformation_gradient: &DeformationGradient) -> Scalar
    {
        let jacobian = deformation_gradient.determinant();
        let isochoric_left_cauchy_green_deformation = self.calculate_left_cauchy_green_deformation(deformation_gradient)/jacobian.powf(2.0/3.0);
        let gamma = (isochoric_left_cauchy_green_deformation.trace()/3.0/self.get_number_of_links()).sqrt();
        let eta = inverse_langevin(gamma);
        let gamma_0 = (1.0/self.get_number_of_links()).sqrt();
        let eta_0 = inverse_langevin(gamma_0);
        self.get_shear_modulus() * self.get_number_of_links() * (gamma * eta - gamma_0 * eta_0 - (eta_0*eta.sinh()/(eta*eta_0.sinh())).ln()) + 0.5*self.get_bulk_modulus()*(0.5*(jacobian.powi(2) - 1.0) - jacobian.ln())

    }
    fn new(parameters: ConstitutiveModelParameters<'a>) -> Self
    {
        Self
        {
            parameters
        }
    }
}

impl<'a> HyperelasticConstitutiveModel for ArrudaBoyceModel<'a>
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
