use crate::math::
{
    TensorRank1Trait,
    TensorRank2Trait
};
use super::
{
    CurrentCoordinate,
    DeformationGradient,
    RotationCurrentConfiguration,
    RotationReferenceConfiguration,
    ReferenceCoordinate,
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

pub fn get_deformation_gradient_rotated() -> DeformationGradient
{
    get_rotation_current_configuration() * get_deformation_gradient() * get_rotation_reference_configuration().transpose()
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

pub fn get_translation_current_configuration() -> CurrentCoordinate
{
    CurrentCoordinate::new([1.1, 2.2, 3.3])
}

pub fn get_translation_reference_configuration() -> ReferenceCoordinate
{
    ReferenceCoordinate::new([4.4, 5.5, 6.6])
}

#[test]
fn size()
{
    assert_eq!(
        std::mem::size_of::<Scalar>(),
        std::mem::size_of::<f64>()
    );
    assert_eq!(
        std::mem::size_of::<DeformationGradient>(),
        std::mem::size_of::<[[f64; 3]; 3]>()
    );
}