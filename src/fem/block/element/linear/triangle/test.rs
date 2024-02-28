use crate::fem::block::element::
{
    test::test_finite_element,
    linear::
    {
        test::test_linear_finite_element
    }
};
use super::*;

fn get_reference_coordinates() -> ReferenceNodalCoordinates<N>
{
    ReferenceNodalCoordinates::new([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0]
    ])
}

use crate::
{
    constitutive::solid::hyperelastic::
    {
        NeoHookean,
        test::NEOHOOKEANPARAMETERS
    },
    mechanics::
    {
        RotationCurrentConfiguration
    }
};

fn get_deformation_gradient_surface() -> DeformationGradient
{
    // needs to be simplified and precise most likely to get tests to pass
    let rotation = RotationCurrentConfiguration::new([
        [-0.96152505,  0.15428784,  0.22729901],
        [-0.26402731, -0.29043967, -0.91974691],
        [-0.07588912, -0.94437284,  0.32000122]
    ]);
    let deformation_gradient = DeformationGradient::new([
        [0.61926467, 0.20573410, 0.0],
        [0.32283599, 0.98534316, 0.0],
        [0.00000000, 0.00000000, 1.0]
    ]);
    rotation * deformation_gradient
}

#[test]
fn temporary_1()
{
    let element = Triangle::<NeoHookean>::new(
        NEOHOOKEANPARAMETERS,
        get_reference_coordinates()
    );
    element.calculate_deformation_gradient(
        &(get_deformation_gradient_surface()*get_reference_coordinates())
    ).iter().zip(get_deformation_gradient_surface().iter()).for_each(|(f_i, ff_i)|
        f_i.iter().zip(ff_i.iter()).for_each(|(f_ij, ff_ij)|
            assert!((f_ij/ff_ij - 1.0).abs() < 1e-8)
        )
    );
}

#[test]
fn test_finite_element()
{
    todo!()
}

#[test]
fn test_linear_finite_element()
{
    todo!()
}