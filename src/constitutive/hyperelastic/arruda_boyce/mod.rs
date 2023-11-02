#[cfg(test)]
mod test;

use super::*;

pub struct ArrudaBoyceModel<'a>
{
    parameters: ConstitutiveModelParameters<'a>
}

impl<'a> ArrudaBoyceModel<'a>
{
    fn get_number_density(&self) -> &Scalar
    {
        &self.parameters[2]
    }
    fn get_number_of_links(&self) -> &Scalar
    {
        &self.parameters[3]
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
        todo!()
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
