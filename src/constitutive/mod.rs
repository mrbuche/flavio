#[cfg(test)]
pub mod test;

mod hyperelastic;

pub use hyperelastic::
{
    HyperelasticConstitutiveModel,
    gent::GentModel,
    mooney_rivlin::MooneyRivlinModel,
    neo_hookean::NeoHookeanModel,
    yeoh::YeohModel
};

use crate::
{
    math::
    {
        ContractSecondIndexWithFirstIndexOf,
        Convert,
        TensorRank2Trait,
        TensorRank4Trait
    },
    mechanics::
    {
        CauchyStress,
        CauchyTangentStiffness,
        DeformationGradient,
        DeformationGradient1,
        DeformationGradient2,
        FirstPiolaKirchoffStress,
        FirstPiolaKirchoffStress2,
        FirstPiolaKirchoffTangentStiffness,
        FirstPiolaKirchoffTangentStiffness2,
        LeftCauchyGreenDeformation,
        MandelStress,
        Scalar
    }
};

pub type ConstitutiveModelParameters<'a> = &'a [Scalar];

pub trait ConstitutiveModel<'a>
{
    fn calculate_cauchy_stress(&self, deformation_gradient: &DeformationGradient) -> CauchyStress;
    fn calculate_cauchy_tangent_stiffness(&self, deformation_gradient: &DeformationGradient) -> CauchyTangentStiffness;
    fn calculate_left_cauchy_green_deformation(&self, deformation_gradient: &DeformationGradient) -> LeftCauchyGreenDeformation
    {
        deformation_gradient.iter()
        .map(|deformation_gradient_i|
            deformation_gradient.iter()
            .map(|deformation_gradient_j|
                deformation_gradient_i * deformation_gradient_j
            ).collect()
        ).collect()
    }
    fn calculate_first_piola_kirchoff_stress(&self, deformation_gradient: &DeformationGradient) -> FirstPiolaKirchoffStress
    {
        self.calculate_cauchy_stress(deformation_gradient)*deformation_gradient.inverse_transpose()*deformation_gradient.determinant()
    }
    fn calculate_first_piola_kirchoff_tangent_stiffness(&self, deformation_gradient: &DeformationGradient) -> FirstPiolaKirchoffTangentStiffness
    {
        let deformation_gradient_inverse_transpose = deformation_gradient.inverse_transpose();
        let first_piola_kirchoff_stress = self.calculate_first_piola_kirchoff_stress(deformation_gradient);
        self.calculate_cauchy_tangent_stiffness(deformation_gradient).contract_second_index_with_first_index_of(&deformation_gradient_inverse_transpose)*deformation_gradient.determinant() + FirstPiolaKirchoffTangentStiffness::dyad_ij_kl(&first_piola_kirchoff_stress, &deformation_gradient_inverse_transpose) - FirstPiolaKirchoffTangentStiffness::dyad_il_kj(&first_piola_kirchoff_stress, &deformation_gradient_inverse_transpose)
    }
    fn calculate_helmholtz_free_energy_density(&self, deformation_gradient: &DeformationGradient) -> Scalar;
    fn new(parameters: ConstitutiveModelParameters<'a>) -> Self;
}

pub trait CompositeConstitutiveModel<C1, C2>
{
    fn construct(constitutive_model_1: C1, constitutive_model_2: C2) -> Self;
    fn get_constitutive_model_1(&self) -> &C1;
    fn get_constitutive_model_2(&self) -> &C2;
}

pub trait CompositeConstitutiveModelMultiplicativeDecomposition<'a, C1, C2>:
    CompositeConstitutiveModel<C1, C2>
where
    C1: ConstitutiveModel<'a>,
    C2: ConstitutiveModel<'a>
{
    fn calculate_mandel_stress(&self, deformation_gradient_1: &DeformationGradient1) -> MandelStress
    {
        deformation_gradient_1.transpose()*self.get_constitutive_model_1().calculate_first_piola_kirchoff_stress(&deformation_gradient_1.convert()).convert()
    }
    fn calculate_objective(&self, deformation_gradient: &DeformationGradient, deformation_gradient_2: &DeformationGradient2) -> Scalar;
    fn calculate_residual(&self, deformation_gradient: &DeformationGradient, deformation_gradient_2: &DeformationGradient2) -> FirstPiolaKirchoffStress2;
    fn calculate_residual_tangent(&self, deformation_gradient: &DeformationGradient, deformation_gradient_2: &DeformationGradient2) -> FirstPiolaKirchoffTangentStiffness2;
}
