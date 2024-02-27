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

use crate::constitutive::solid::hyperelastic::
{
    NeoHookean,
    test::NEOHOOKEANPARAMETERS
};

#[test]
fn temporary_1()
{
    let element = Triangle::<NeoHookean>::new(
        NEOHOOKEANPARAMETERS,
        get_reference_coordinates()
    );
    element.calculate_deformation_gradient(
        &get_reference_coordinates().convert()
    ).iter().for_each(|f_i|
        f_i.iter().for_each(|f_ij|
            println!("{}", f_ij)
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