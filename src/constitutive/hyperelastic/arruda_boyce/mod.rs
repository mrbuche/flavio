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
        todo!()
    }
    fn calculate_cauchy_tangent_stiffness(&self, deformation_gradient: &DeformationGradient) -> CauchyTangentStiffness
    {
        todo!()
    }
    fn calculate_helmholtz_free_energy_density(&self, deformation_gradient: &DeformationGradient) -> Scalar
    {
        let jacobian = deformation_gradient.determinant();
        let isochoric_left_cauchy_green_deformation = self.calculate_left_cauchy_green_deformation(deformation_gradient)/jacobian.powf(2.0/3.0);
        let gamma = (isochoric_left_cauchy_green_deformation.trace()/3.0/self.get_number_of_links()).sqrt();
        let eta = inverse_langevin(gamma);
        self.get_shear_modulus() * self.get_number_of_links() * (gamma * eta - (eta.sinh() / eta).ln()) + 0.5 * self.get_bulk_modulus() * (0.5 * (jacobian.powi(2) - 1.0) - jacobian.ln())

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
