#[cfg(test)]
pub mod test;

pub mod elastic;
pub mod hyperelastic;

use crate::
{
    math::
    {
        ContractFirstSecondIndicesWithSecondIndicesOf,
        ContractSecondIndexWithFirstIndexOf,
        TensorRank2Trait,
        TensorRank4Trait
    },
    mechanics::
    {
        CauchyStress,
        CauchyTangentStiffness,
        DeformationGradient,
        FirstPiolaKirchoffStress,
        FirstPiolaKirchoffTangentStiffness,
        LeftCauchyGreenDeformation,
        RightCauchyGreenDeformation,
        Scalar,
        SecondPiolaKirchoffStress,
        SecondPiolaKirchoffTangentStiffness
    }
};

/// Array of constitutive model parameters.
pub type ConstitutiveModelParameters<'a> = &'a [Scalar];

/// Required methods for constitutive models.
pub trait ConstitutiveModel<'a>
{
    /// Calculates and returns the Cauchy stress.
    ///
    /// ```math
    /// \boldsymbol{\sigma} = \frac{1}{J}\frac{\partial a}{\partial\mathbf{F}}\cdot\mathbf{F}^T
    /// ```
    fn calculate_cauchy_stress(&self, deformation_gradient: &DeformationGradient) -> CauchyStress
    {
        deformation_gradient*self.calculate_second_piola_kirchoff_stress(deformation_gradient)*deformation_gradient.transpose()/deformation_gradient.determinant()
    }
    /// Calculates and returns the tangent stiffness associated with the Cauchy stress.
    ///
    /// ```math
    /// \mathcal{T}_{ijkL} = \frac{\partial\sigma_{ij}}{\partial F_{kL}} = J^{-1} \mathcal{G}_{MNkL} F_{iM} F_{jN} - \sigma_{ij} F_{kL}^{-T} + \left(\delta_{jk}\sigma_{is} + \delta_{ik}\sigma_{js}\right)F_{sL}^{-T}
    /// ```
    fn calculate_cauchy_tangent_stiffness(&self, deformation_gradient: &DeformationGradient) -> CauchyTangentStiffness
    {
        let deformation_gradient_inverse_transpose = deformation_gradient.inverse_transpose();
        let cauchy_stress = self.calculate_cauchy_stress(deformation_gradient);
        let identity = CauchyStress::identity();
        let some_stress = &cauchy_stress * &deformation_gradient_inverse_transpose;
        self.calculate_second_piola_kirchoff_tangent_stiffness(deformation_gradient).contract_first_second_indices_with_second_indices_of(deformation_gradient, deformation_gradient)/deformation_gradient.determinant() - CauchyTangentStiffness::dyad_ij_kl(&cauchy_stress, &deformation_gradient_inverse_transpose) + CauchyTangentStiffness::dyad_il_kj(&some_stress, &identity) + CauchyTangentStiffness::dyad_ik_jl(&identity, &some_stress)
    }
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
    /// Calculates and returns the first Piola-Kirchoff stress.
    ///
    /// ```math
    /// \mathbf{P} = \frac{\partial a}{\partial\mathbf{F}} = J\boldsymbol{\sigma}\cdot\mathbf{F}^{-T}
    /// ```
    fn calculate_first_piola_kirchoff_stress(&self, deformation_gradient: &DeformationGradient) -> FirstPiolaKirchoffStress
    {
        self.calculate_cauchy_stress(deformation_gradient)*deformation_gradient.inverse_transpose()*deformation_gradient.determinant()
    }
    /// Calculates and returns the tangent stiffness associated with the first Piola-Kirchoff stress.
    ///
    /// ```math
    /// \mathcal{C}_{iJkL} = \frac{\partial P_{iJ}}{\partial F_{kL}} = J \mathcal{T}_{iskL} F_{sJ}^{-T} + P_{iJ} F_{kL}^{-T} - P_{iL} F_{kJ}^{-T}
    /// ```
    fn calculate_first_piola_kirchoff_tangent_stiffness(&self, deformation_gradient: &DeformationGradient) -> FirstPiolaKirchoffTangentStiffness
    {
        let deformation_gradient_inverse_transpose = deformation_gradient.inverse_transpose();
        let first_piola_kirchoff_stress = self.calculate_first_piola_kirchoff_stress(deformation_gradient);
        self.calculate_cauchy_tangent_stiffness(deformation_gradient).contract_second_index_with_first_index_of(&deformation_gradient_inverse_transpose)*deformation_gradient.determinant() + FirstPiolaKirchoffTangentStiffness::dyad_ij_kl(&first_piola_kirchoff_stress, &deformation_gradient_inverse_transpose) - FirstPiolaKirchoffTangentStiffness::dyad_il_kj(&first_piola_kirchoff_stress, &deformation_gradient_inverse_transpose)
    }
    /// Calculates and returns the right Cauchy-Green deformation.
    ///
    /// ```math
    /// \mathbf{C} = \mathbf{F}^T\cdot\mathbf{F}
    /// ```
    fn calculate_right_cauchy_green_deformation(&self, deformation_gradient: &DeformationGradient) -> RightCauchyGreenDeformation
    {
        let deformation_gradient_transpose = deformation_gradient.transpose();
        deformation_gradient_transpose.iter()
        .map(|deformation_gradient_transpose_i|
            deformation_gradient_transpose.iter()
            .map(|deformation_gradient_transpose_j|
                deformation_gradient_transpose_i * deformation_gradient_transpose_j
            ).collect()
        ).collect()
    }
    /// Calculates and returns the second Piola-Kirchoff stress.
    ///
    /// ```math
    /// \mathbf{S} = \frac{\partial a}{\partial\mathbf{E}} = \mathbf{F}^{-1}\cdot\mathbf{P}
    /// ```
    fn calculate_second_piola_kirchoff_stress(&self, deformation_gradient: &DeformationGradient) -> SecondPiolaKirchoffStress
    {
        deformation_gradient.inverse()*self.calculate_cauchy_stress(deformation_gradient)*deformation_gradient.inverse_transpose()*deformation_gradient.determinant()
    }
    /// Calculates and returns the tangent stiffness associated with the second Piola-Kirchoff stress.
    ///
    /// ```math
    /// \mathcal{G}_{IJkL} = \frac{\partial S_{IJ}}{\partial F_{kL}} = \mathcal{C}_{mJkL}F_{mI}^{-T} - S_{LJ}F_{kI}^{-T} = J \mathcal{T}_{mnkL} F_{mI}^{-T} F_{nJ}^{-T} + S_{IJ} F_{kL}^{-T} - S_{IL} F_{kJ}^{-T} -S_{LJ} F_{kI}^{-T}
    /// ```
    fn calculate_second_piola_kirchoff_tangent_stiffness(&self, deformation_gradient: &DeformationGradient) -> SecondPiolaKirchoffTangentStiffness
    {
        let deformation_gradient_inverse_transpose = deformation_gradient.inverse_transpose();
        let deformation_gradient_inverse = deformation_gradient_inverse_transpose.transpose();
        let second_piola_kirchoff_stress = self.calculate_second_piola_kirchoff_stress(deformation_gradient);
        self.calculate_cauchy_tangent_stiffness(deformation_gradient).contract_first_second_indices_with_second_indices_of(&deformation_gradient_inverse, &deformation_gradient_inverse)*deformation_gradient.determinant() + SecondPiolaKirchoffTangentStiffness::dyad_ij_kl(&second_piola_kirchoff_stress, &deformation_gradient_inverse_transpose) - SecondPiolaKirchoffTangentStiffness::dyad_il_kj(&second_piola_kirchoff_stress, &deformation_gradient_inverse_transpose) - SecondPiolaKirchoffTangentStiffness::dyad_ik_jl(&deformation_gradient_inverse, &second_piola_kirchoff_stress)
    }
    /// Constructs and returns a new constitutive model.
    fn new(parameters: ConstitutiveModelParameters<'a>) -> Self;
}

/// Required methods for composite constitutive models.
pub trait CompositeConstitutiveModel<C1, C2>
{
    /// Constructs and returns a new composite constitutive model.
    fn construct(constitutive_model_1: C1, constitutive_model_2: C2) -> Self;
    /// Returns a reference to the first constitutive model.
    fn get_constitutive_model_1(&self) -> &C1;
    /// Returns a reference to the second constitutive model.
    fn get_constitutive_model_2(&self) -> &C2;
}
