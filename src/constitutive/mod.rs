#[cfg(test)]
pub mod test;

/// Hyperelastic constitutive models.
pub mod hyperelastic;

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
    /// Calculates and returns the Cauchy stress.
    ///
    /// ```math
    /// \boldsymbol{\sigma} = \frac{1}{J}\frac{\partial a}{\partial\mathbf{F}}\cdot\mathbf{F}^T
    /// ```
    fn calculate_cauchy_stress(&self, deformation_gradient: &DeformationGradient) -> CauchyStress;
    /// Calculates and returns the tangent stiffness associated with the Cauchy stress.
    ///
    /// ```math
    /// \boldsymbol{\mathcal{T}} = \frac{\partial\boldsymbol{\sigma}}{\partial\mathbf{F}}
    /// ```
    fn calculate_cauchy_tangent_stiffness(&self, deformation_gradient: &DeformationGradient) -> CauchyTangentStiffness;
    /// Calculates and returns the left Cauchy-Green deformation.
    ///
    /// ```math
    /// \mathbf{B} = \mathbf{F}\cdot\mathbf{F}^T
    /// ```
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
    /// Calculates and returns the tangent stiffness associated with the first Piola-Kirchoff stress.
    ///
    /// ```math
    /// \mathbf{P} =
    /// \frac{\partial a}{\partial\mathbf{F}} =
    /// J\boldsymbol{\sigma}\cdot\mathbf{F}^{-T}
    /// ```
    fn calculate_first_piola_kirchoff_stress(&self, deformation_gradient: &DeformationGradient) -> FirstPiolaKirchoffStress
    {
        self.calculate_cauchy_stress(deformation_gradient)*deformation_gradient.inverse_transpose()*deformation_gradient.determinant()
    }
    /// Calculates and returns the tangent stiffness associated with the first Piola-Kirchoff stress.
    ///
    /// ```math
    /// \mathcal{C}_{iJkL} =
    /// \frac{\partial P_{iJ}}{\partial F_{kL}} = 
    /// J \mathcal{T}_{iskL} F_{sJ}^{-T} + P_{iJ} F_{kL}^{-T} - P_{iL} F_{kJ}^{-T}
    /// ```
    fn calculate_first_piola_kirchoff_tangent_stiffness(&self, deformation_gradient: &DeformationGradient) -> FirstPiolaKirchoffTangentStiffness
    {
        let deformation_gradient_inverse_transpose = deformation_gradient.inverse_transpose();
        let first_piola_kirchoff_stress = self.calculate_first_piola_kirchoff_stress(deformation_gradient);
        self.calculate_cauchy_tangent_stiffness(deformation_gradient).contract_second_index_with_first_index_of(&deformation_gradient_inverse_transpose)*deformation_gradient.determinant() + FirstPiolaKirchoffTangentStiffness::dyad_ij_kl(&first_piola_kirchoff_stress, &deformation_gradient_inverse_transpose) - FirstPiolaKirchoffTangentStiffness::dyad_il_kj(&first_piola_kirchoff_stress, &deformation_gradient_inverse_transpose)
    }
    /// Calculates and returns the Helmholtz free energy density.
    ///
    /// ```math
    /// a = a(\mathbf{F})
    /// ```
    fn calculate_helmholtz_free_energy_density(&self, deformation_gradient: &DeformationGradient) -> Scalar;
    /// Creates and returns a new constitutive model.
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
