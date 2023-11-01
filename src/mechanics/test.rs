use crate::math::TensorRank2Trait;
use super::
{
    DeformationGradient,
    DeformationGradient2,
    RotationCurrentConfiguration,
    RotationIntermediateConfiguration,
    RotationReferenceConfiguration,
    Scalar
};

pub fn get_deformation_gradient() -> DeformationGradient
{
    DeformationGradient::new([
        [0.63595746, 0.69157849, 0.71520784],
        [0.80589604, 0.83687323, 0.19312595],
        [0.05387420, 0.86551549, 0.41880244]
    ])
}

pub fn get_deformation_gradient_2() -> DeformationGradient2
{
    DeformationGradient2::new([
        [0.97902607, 0.09463390, 0.48164011],
        [0.19335510, 0.89379463, 0.49423509],
        [0.02648444, 0.28479496, 0.34719272]
    ])
}

pub fn get_deformation_gradient_rotated() -> DeformationGradient
{
    get_rotation_current_configuration() * get_deformation_gradient() * get_rotation_reference_configuration().transpose()
}

pub fn get_deformation_gradient_2_rotated() -> DeformationGradient2
{
    get_rotation_intermediate_configuration() * get_deformation_gradient_2() * get_rotation_reference_configuration().transpose()
}

pub fn get_rotation_current_configuration() -> RotationCurrentConfiguration
{
    let sqrt_2 = (2.0 as Scalar).sqrt();
    let sqrt_3 = (3.0 as Scalar).sqrt();
    let sqrt_6 = (6.0 as Scalar).sqrt();
    RotationCurrentConfiguration::new([
        [0.25*sqrt_2, -0.25*sqrt_6, 0.5*sqrt_2],
        [0.125*sqrt_2 + 0.75, -0.125*sqrt_6 + 0.25*sqrt_3, -0.25*sqrt_2],
        [-0.125*sqrt_6 + 0.25*sqrt_3, 0.25 + 0.375*sqrt_2, 0.25*sqrt_6]
    ])
}

pub fn get_rotation_intermediate_configuration() -> RotationIntermediateConfiguration
{
    RotationIntermediateConfiguration::new([
        [ 0.36,  0.48, -0.80],
        [-0.80,  0.60,  0.00],
        [ 0.48,  0.64,  0.60]
    ])
}

pub fn get_rotation_reference_configuration() -> RotationReferenceConfiguration
{
    let sqrt_2 = (2.0 as Scalar).sqrt();
    let sqrt_3 = (3.0 as Scalar).sqrt();
    let sqrt_6 = (6.0 as Scalar).sqrt();
    RotationReferenceConfiguration::new([
        [0.25*sqrt_6, -0.25*sqrt_6, 0.5],
        [0.125*sqrt_6 + 0.25*sqrt_2, -0.125*sqrt_6 + 0.25*sqrt_2, -0.75],
        [-0.125*sqrt_2 + 0.25*sqrt_6, 0.125*sqrt_2 + 0.25*sqrt_6, 0.25*sqrt_3]
    ])
}