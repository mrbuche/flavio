//! Solid constitutive models.

pub mod elastic;
pub mod elastic_hyperviscous;
pub mod hyperelastic;
pub mod hyperviscoelastic;
pub mod thermoelastic;
pub mod thermohyperelastic;
pub mod viscoelastic;

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
        CauchyRateTangentStiffness,
        DeformationGradient,
        DeformationGradientRate,
        FirstPiolaKirchoffStress,
        FirstPiolaKirchoffTangentStiffness,
        FirstPiolaKirchoffRateTangentStiffness,
        LeftCauchyGreenDeformation,
        RightCauchyGreenDeformation,
        Scalar,
        SecondPiolaKirchoffStress,
        SecondPiolaKirchoffTangentStiffness,
        SecondPiolaKirchoffRateTangentStiffness
    }
};
use super::
{
    Constitutive,
    Parameters
};

/// Required methods for solid constitutive models.
pub trait Solid<'a>
where
    Self: Constitutive<'a>
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
    /// Returns the bulk modulus.
    fn get_bulk_modulus(&self) -> &Scalar;
    /// Returns the shear modulus.
    fn get_shear_modulus(&self) -> &Scalar;
}
