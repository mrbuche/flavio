#[cfg(test)]
pub mod test;

pub mod elastic;
pub mod hyperelastic;
pub mod thermoelastic;
pub mod thermohyperelastic;

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
        SecondPiolaKirchoffTangentStiffness,
        Temperature
    }
};

/// Array of constitutive model parameters.
pub type ConstitutiveModelParameters<'a> = &'a [Scalar];

// why give ConstitutiveModel any methods if not going to ever call them until the type of model (i.e. hyperelastic) is known?
//     -> doesnt seem necessary or useful now that external (and later, internal) state variables vary among models
//     -> new() is probably the only thing you need besides the overal Trait just for typing fields in FEM
//     -> can keep stuff like calculate_left_cauchy_green_deformation()
//     -> can probably avoid this external_state_variables: &Y stuff
//         -> which will let you get back to default impls
//     -> see in FEM, never call these methods (you really cant) until you know that its hyperelastic

/// Required methods for constitutive models.
pub trait ConstitutiveModel<'a>
{
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
