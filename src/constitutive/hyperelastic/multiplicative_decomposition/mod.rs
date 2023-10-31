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
        let mut deformation_gradient_2 = DeformationGradient2::identity();
        let mut residual = self.calculate_residual(deformation_gradient, &deformation_gradient_2);
        let mut residual_norm = residual.norm();
        while residual_norm > 1e-3
        {
            residual = self.calculate_residual(deformation_gradient, &deformation_gradient_2);
            residual_norm = residual.norm();
            deformation_gradient_2 -= self.calculate_residual_tangent(deformation_gradient, &deformation_gradient_2).inverse().contract_third_fourth_indices_with_first_second_indices_of(&residual) * 0.999;
            if deformation_gradient_2.determinant() < crate::REL_TOL
            {
                panic!()
            }
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
        let (deformation_gradient_2_inverse, deformation_gradient_2_determinant) = self.calculate_deformation_gradient_2(deformation_gradient).inverse_and_determinant();
        self.get_constitutive_model_1().calculate_cauchy_stress(&(deformation_gradient * deformation_gradient_2_inverse).into())/deformation_gradient_2_determinant
    }
    fn calculate_cauchy_tangent_stiffness(&self, deformation_gradient: &DeformationGradient) -> CauchyTangentStiffness
    {
        let (deformation_gradient_2_inverse_transpose, deformation_gradient_2_determinant) = self.calculate_deformation_gradient_2(deformation_gradient).inverse_transpose_and_determinant();
        <CauchyTangentStiffness as Into<CauchyTangentStiffness1>>::into(self.get_constitutive_model_1().calculate_cauchy_tangent_stiffness(&(deformation_gradient * deformation_gradient_2_inverse_transpose.transpose()).into())) * (deformation_gradient_2_inverse_transpose / deformation_gradient_2_determinant)
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
        self.get_constitutive_model_1().calculate_helmholtz_free_energy_density(&(deformation_gradient * deformation_gradient_2.inverse()).into()) + self.get_constitutive_model_2().calculate_helmholtz_free_energy_density(&deformation_gradient_2.convert())
    }
    fn calculate_residual(&self, deformation_gradient: &DeformationGradient, deformation_gradient_2: &DeformationGradient2) -> FirstPiolaKirchoffStress2
    {
        let deformation_gradient_1 = deformation_gradient * deformation_gradient_2.inverse();
        <FirstPiolaKirchoffStress as Into<FirstPiolaKirchoffStress2>>::into(self.get_constitutive_model_2().calculate_first_piola_kirchoff_stress(&deformation_gradient_2.convert())) - deformation_gradient_1.transpose() * <FirstPiolaKirchoffStress as Into<FirstPiolaKirchoffStress1>>::into(self.get_constitutive_model_1().calculate_first_piola_kirchoff_stress(&deformation_gradient_1.into())) * deformation_gradient_2.inverse_transpose()
    }
    fn calculate_residual_tangent(&self, deformation_gradient: &DeformationGradient, deformation_gradient_2: &DeformationGradient2) -> FirstPiolaKirchoffTangentStiffness2
    {
        let deformation_gradient_1 = deformation_gradient * deformation_gradient_2.inverse();
        let deformation_gradient_1_convert = deformation_gradient_1.convert();
        let deformation_gradient_2_inverse_transpose = deformation_gradient_2.inverse_transpose();
        let first_piola_kirchoff_stress_1: FirstPiolaKirchoffStress1 = self.get_constitutive_model_1().calculate_first_piola_kirchoff_stress(&deformation_gradient_1_convert).into();
        <FirstPiolaKirchoffTangentStiffness as Into<FirstPiolaKirchoffTangentStiffness2>>::into(self.get_constitutive_model_2().calculate_first_piola_kirchoff_tangent_stiffness(&deformation_gradient_2.convert())) + <FirstPiolaKirchoffTangentStiffness as Into<FirstPiolaKirchoffTangentStiffness1>>::into(self.get_constitutive_model_1().calculate_first_piola_kirchoff_tangent_stiffness(&deformation_gradient_1_convert)).contract_all_indices_with_first_indices_of(&deformation_gradient_1, &deformation_gradient_2_inverse_transpose, &deformation_gradient_1, &deformation_gradient_2_inverse_transpose) + FirstPiolaKirchoffTangentStiffness2::dyad_il_kj(&deformation_gradient_2_inverse_transpose, &(deformation_gradient_1.transpose() * &first_piola_kirchoff_stress_1 * &deformation_gradient_2_inverse_transpose)) + FirstPiolaKirchoffTangentStiffness2::dyad_il_kj(&(deformation_gradient_1.transpose() * &first_piola_kirchoff_stress_1 * &deformation_gradient_2_inverse_transpose), &deformation_gradient_2_inverse_transpose)
    }
}