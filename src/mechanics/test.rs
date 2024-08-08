use crate::
{
    math::
    {
        TensorRank1Trait,
        TensorRank2Trait
    },
    test::assert_eq_within_tols
};
use super::
{
    CurrentCoordinate,
    CurrentVelocity,
    DeformationGradient,
    DeformationGradientRate,
    FrameSpin,
    RotationCurrentConfiguration,
    RotationRateCurrentConfiguration,
    RotationReferenceConfiguration,
    ReferenceCoordinate,
    TemperatureGradient,
    Scalar,
    IDENTITY,
    IDENTITY_00,
    IDENTITY_10
};

pub fn get_deformation_gradient() -> DeformationGradient
{
    DeformationGradient::new([
        [0.63595746, 0.69157849, 0.71520784],
        [0.80589604, 0.83687323, 0.19312595],
        [0.05387420, 0.86551549, 0.41880244]
    ])
}

pub fn get_deformation_gradient_rate() -> DeformationGradientRate
{
    DeformationGradientRate::new([
        [0.17414455, 0.97269465, 0.87856299],
        [0.96651379, 0.84670298, 0.66739030],
        [0.16205052, 0.85112927, 0.38711266]
    ])
}

pub fn get_deformation_gradient_rotated() -> DeformationGradient
{
    get_rotation_current_configuration() * get_deformation_gradient() * get_rotation_reference_configuration().transpose()
}

pub fn get_deformation_gradient_rotated_undeformed() -> DeformationGradient
{
    get_rotation_current_configuration() * IDENTITY_10 * get_rotation_reference_configuration().transpose()
}

pub fn get_deformation_gradient_rate_rotated() -> DeformationGradientRate
{
    (get_rotation_current_configuration() * get_deformation_gradient_rate()
        + get_rotation_rate_current_configuration() * get_deformation_gradient()
    ) * get_rotation_reference_configuration().transpose()
}

pub fn get_deformation_gradient_rate_rotated_undeformed() -> DeformationGradientRate
{
    get_rotation_rate_current_configuration() * get_deformation_gradient() * get_rotation_reference_configuration().transpose()
}

pub fn get_frame_spin() -> FrameSpin
{
    FrameSpin::new([
        [ 0.0,  0.1, -0.2],
        [-0.1,  0.0,  0.8],
        [ 0.2, -0.8,  0.0]
    ])
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

pub fn get_rotation_rate_current_configuration() -> RotationRateCurrentConfiguration
{
    get_frame_spin() * get_rotation_current_configuration()
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

pub fn get_translation_rate_current_configuration() -> CurrentVelocity
{
    CurrentVelocity::new([0.4, -0.7, 1.2])
}

pub fn get_translation_reference_configuration() -> ReferenceCoordinate
{
    ReferenceCoordinate::new([4.4, 5.5, 6.6])
}

pub fn get_temperature() -> Scalar
{
    123.0
}

pub fn get_temperature_gradient() -> TemperatureGradient
{
    TemperatureGradient::new([12.3, -5.0, 8.8])
}

#[test]
fn frame_spin_tensor()
{
    get_frame_spin().iter()
    .zip(get_frame_spin().transpose().iter())
    .for_each(|(frame_spin_i, frame_spin_transpose_i)|
        frame_spin_i.iter()
        .zip(frame_spin_transpose_i.iter())
        .for_each(|(frame_spin_ij, frame_spin_transpose_ij)|
            assert_eq!(frame_spin_ij, &(-frame_spin_transpose_ij))
        )
    )
}

#[test]
fn rotation_current_configuration()
{
    (get_rotation_current_configuration() *
        get_rotation_current_configuration().transpose()
    ).iter().zip(IDENTITY.iter())
    .for_each(|(res_i, identity_i)|
        res_i.iter().zip(identity_i.iter())
        .for_each(|(res_ij, identity_ij)|
            assert_eq_within_tols(res_ij, identity_ij)
        )
    )
}

#[test]
fn rotation_reference_configuration()
{
    (get_rotation_reference_configuration() *
        get_rotation_reference_configuration().transpose()
    ).iter().zip(IDENTITY_00.iter())
    .for_each(|(res_i, identity_i)|
        res_i.iter().zip(identity_i.iter())
        .for_each(|(res_ij, identity_ij)|
            assert_eq_within_tols(res_ij, identity_ij)
        )
    )
}

#[test]
fn rotation_rate_current_configuration()
{
    get_rotation_rate_current_configuration().iter()
    .zip((get_frame_spin() * get_rotation_current_configuration()).iter())
    .for_each(|(rotation_rate_i, res_rotation_rate_i)|
        rotation_rate_i.iter()
        .zip(res_rotation_rate_i.iter())
        .for_each(|(rotation_rate_ij, res_rotation_rate_ij)|
            assert_eq_within_tols(rotation_rate_ij, res_rotation_rate_ij)
        )
    )
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
