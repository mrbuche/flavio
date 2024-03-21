use crate::fem::block::
{
    element::composite::surface::test::
    {
        setup_for_test_composite_surface_element_with_constitutive_model,
        test_composite_surface_element
    },
    test::test_finite_element_block
};
use super::*;

const D: usize = 16;
const E: usize = 18;

fn get_connectivity() -> Connectivity<E, N>
{
    todo!()
}

fn get_coordinates_block() -> NodalCoordinates<D>
{
    todo!()
}

fn get_reference_coordinates() -> ReferenceNodalCoordinates<N>
{
    ReferenceNodalCoordinates::new([
        // [0.0, 0.0, 0.0],
        // [1.0, 0.0, 0.0],
        // [0.0, 1.0, 0.0],
        // [0.0, 0.5, 0.0],
        // [0.5, 0.5, 0.0],
        // [0.0, 0.5, 0.0]
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.0, 0.5, 0.0],
        [0.5, 0.0, 0.0]
    ])
}

#[test]
fn NEED_TO_MAKE_NODE_NUMBERING_AND_SHAPE_FUNCTIONS_CONSISTENT_TO_GET_RIGHT_ANSWERS()
{
    todo!()
}

fn get_reference_coordinates_block() -> ReferenceNodalCoordinates<D>
{
    todo!()
}

fn get_velocities_block() -> NodalVelocities<D>
{
    todo!()
}

test_composite_surface_element!(Triangle);
test_finite_element_block!(Triangle);