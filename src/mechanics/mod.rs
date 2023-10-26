use crate::math::
{
    TensorRank0,
    TensorRank2
};

pub type CauchyStress = TensorRank2<3, 1, 1>;
pub type DeformationGradient = TensorRank2<3, 1, 0>;
pub type DeformationGradient1 = TensorRank2<3, 1, 2>;
pub type DeformationGradient2 = TensorRank2<3, 2, 0>;
pub type FirstPiolaKirchoffStress = TensorRank2<3, 1, 0>;
pub type Scalar = TensorRank0;
