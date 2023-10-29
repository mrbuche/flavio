#[cfg(test)]
mod test;

use super::*;

pub struct CompositeHyperelasticConstitutiveModelMultiplicativeDecomposition<C1, C2>
{
    hyperelastic_constitutive_model_1: C1,
    hyperelastic_constitutive_model_2: C2
}

impl<'a, C1, C2> CompositeHyperelasticConstitutiveModelMultiplicativeDecomposition<C1, C2>
where
    C1: ConstitutiveModel<'a>,
    C2: ConstitutiveModel<'a>
{
    fn calculate_deformation_gradient_2(&self, deformation_gradient: &DeformationGradient) -> DeformationGradient2
    {
        let step_size = 0.999;
        let tolerance = 1e-8;
        let mut deformation_gradient_2 = DeformationGradient2::identity();
        let mut residual = self.calculate_residual(deformation_gradient, &deformation_gradient_2);
        let mut residual_norm = residual.norm();
        while residual_norm > tolerance
        {
            residual = self.calculate_residual(deformation_gradient, &deformation_gradient_2);
            residual_norm = residual.norm();
            deformation_gradient_2 -= self.calculate_residual_tangent(deformation_gradient, &deformation_gradient_2).inverse().contract_third_fourth_indices_with_first_second_indices_of(&residual) * step_size;
        }
        deformation_gradient_2
    }
}

impl<'a, C1, C2> ConstitutiveModel<'a> for CompositeHyperelasticConstitutiveModelMultiplicativeDecomposition<C1, C2>
where
    C1: ConstitutiveModel<'a>,
    C2: ConstitutiveModel<'a>
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
        let deformation_gradient_2 = self.calculate_deformation_gradient_2(deformation_gradient);
        self.calculate_objective(deformation_gradient, &deformation_gradient_2)
    }
    fn new(_parameters: ConstitutiveModelParameters<'a>) -> Self
    {
        panic!()
    }
}

impl<C1, C2> HyperelasticConstitutiveModel for CompositeHyperelasticConstitutiveModelMultiplicativeDecomposition<C1, C2>
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

impl<C1, C2> CompositeConstitutiveModel<C1, C2> for CompositeHyperelasticConstitutiveModelMultiplicativeDecomposition<C1, C2>
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

impl<'a, C1, C2> CompositeConstitutiveModelMultiplicativeDecomposition<'a, C1, C2> for CompositeHyperelasticConstitutiveModelMultiplicativeDecomposition<C1, C2>
where
    C1: ConstitutiveModel<'a>,
    C2: ConstitutiveModel<'a>
{
    fn calculate_objective(&self, deformation_gradient: &DeformationGradient, deformation_gradient_2: &DeformationGradient2) -> Scalar
    {
        todo!()
    }
    fn calculate_residual(&self, deformation_gradient: &DeformationGradient, deformation_gradient_2: &DeformationGradient2) -> FirstPiolaKirchoffStress2
    {
        todo!()
    }
    fn calculate_residual_tangent(&self, deformation_gradient: &DeformationGradient, deformation_gradient_2: &DeformationGradient2) -> FirstPiolaKirchoffTangentStiffness2
    {
        todo!()
    }
}