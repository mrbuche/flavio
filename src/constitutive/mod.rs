#[cfg(test)]
mod test;

mod hyperelastic;

use crate::mechanics::
{
    DeformationGradient,
    LeftCauchyGreenDeformation,
    Scalar
};

pub type ConstitutiveModelParameters<'a> = &'a [Scalar];

pub trait ConstitutiveModel
{
    fn calculate_left_cauchy_green_deformation(&self, deformation_gradient: &DeformationGradient) -> LeftCauchyGreenDeformation;
}
