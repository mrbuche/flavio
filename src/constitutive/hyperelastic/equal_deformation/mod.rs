#[cfg(test)]
mod test;

use super::*;

pub struct CompositeHyperelasticConstitutiveModelEqualDeformation<C1, C2>
{
    hyperelastic_constitutive_model_1: C1,
    hyperelastic_constitutive_model_2: C2
}

impl<'a, C1, C2> ConstitutiveModel<'a> for CompositeHyperelasticConstitutiveModelEqualDeformation<C1, C2>
where
    C1: ConstitutiveModel<'a>,
    C2: ConstitutiveModel<'a>
{
    fn calculate_cauchy_stress(&self, deformation_gradient: &DeformationGradient) -> CauchyStress
    {
        self.get_constitutive_model_1().calculate_cauchy_stress(deformation_gradient) + self.get_constitutive_model_2().calculate_cauchy_stress(deformation_gradient)
    }
    fn calculate_cauchy_tangent_stiffness(&self, deformation_gradient: &DeformationGradient) -> CauchyTangentStiffness
    {
        self.get_constitutive_model_1().calculate_cauchy_tangent_stiffness(deformation_gradient) + self.get_constitutive_model_2().calculate_cauchy_tangent_stiffness(deformation_gradient)
    }
    fn calculate_helmholtz_free_energy_density(&self, deformation_gradient: &DeformationGradient) -> Scalar
    {
        self.get_constitutive_model_1().calculate_helmholtz_free_energy_density(deformation_gradient) + self.get_constitutive_model_2().calculate_helmholtz_free_energy_density(deformation_gradient)
    }
    fn new(_parameters: ConstitutiveModelParameters<'a>) -> Self
    {
        panic!()
    }
}

impl<C1, C2> HyperelasticConstitutiveModel for CompositeHyperelasticConstitutiveModelEqualDeformation<C1, C2>
{
    fn get_bulk_modulus(&self) -> &Scalar
    {
        panic!()
    }
    fn get_shear_modulus(&self) -> &Scalar
    {
        panic!()
    }
}

impl<C1, C2> CompositeConstitutiveModel<C1, C2> for CompositeHyperelasticConstitutiveModelEqualDeformation<C1, C2>
{

    fn construct(hyperelastic_constitutive_model_1: C1, hyperelastic_constitutive_model_2: C2) -> Self
    {
        Self
        {
            hyperelastic_constitutive_model_1,
            hyperelastic_constitutive_model_2
        }
    }
    fn get_constitutive_model_1(&self) -> &C1
    {
        &self.hyperelastic_constitutive_model_1
    }
    fn get_constitutive_model_2(&self) -> &C2
    {
        &self.hyperelastic_constitutive_model_2
    }
}
