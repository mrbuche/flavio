#[cfg(test)]
pub mod test;

use crate::math::
{
    TensorRank0,
    TensorRank2,
    TensorRank4
};

pub type CauchyStress = TensorRank2<3, 1, 1>;
pub type CauchyTangentStiffness = TensorRank4<3, 1, 1, 1, 0>;
pub type CauchyTangentStiffness1 = TensorRank4<3, 1, 1, 1, 2>;
pub type DeformationGradient = TensorRank2<3, 1, 0>;
pub type DeformationGradient1 = TensorRank2<3, 1, 2>;
pub type DeformationGradient2 = TensorRank2<3, 2, 0>;
pub type FirstPiolaKirchoffStress = TensorRank2<3, 1, 0>;
pub type FirstPiolaKirchoffStress1 = TensorRank2<3, 1, 2>;
pub type FirstPiolaKirchoffStress2 = TensorRank2<3, 2, 0>;
pub type FirstPiolaKirchoffTangentStiffness = TensorRank4<3, 1, 0, 1, 0>;
pub type FirstPiolaKirchoffTangentStiffness1 = TensorRank4<3, 1, 2, 1, 2>;
pub type FirstPiolaKirchoffTangentStiffness2 = TensorRank4<3, 2, 0, 2, 0>;
pub type LeftCauchyGreenDeformation = TensorRank2<3, 1, 1>;
pub type MandelStress = TensorRank2<3, 2, 2>;
pub type RotationCurrentConfiguration = TensorRank2<3, 1, 1>;
pub type RotationIntermediateConfiguration = TensorRank2<3, 2, 2>;
pub type RotationReferenceConfiguration = TensorRank2<3, 0, 0>;
pub type Scalar = TensorRank0;
