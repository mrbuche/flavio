use crate::fem::block::element::linear::localization::test::
{
    setup_for_test_linear_surface_element_with_constitutive_model,
    test_linear_localization_element
};
use super::*;

fn get_reference_coordinates() -> ReferenceNodalCoordinates<N>
{
    ReferenceNodalCoordinates::new([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0]
    ])
}

test_linear_localization_element!(Wedge);

#[test]
fn do_block()
{
    todo!()
}

// should use random ref/curr coords for element/block tests

// maybe go back to CurrentNodalCoordinates and stuff and call Coordinates<I> NodalCoordinates<I>