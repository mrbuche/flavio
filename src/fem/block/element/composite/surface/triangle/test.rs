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

const D: usize = 25;
const E: usize = 8;

fn get_connectivity() -> Connectivity<E, N>
{
    [
        [ 0,  2, 16,  3, 17, 21],
        [ 0, 16, 13, 21, 18, 15],
        [ 2,  1,  6,  4,  7, 22],
        [ 2,  6, 16, 22, 19, 17],
        [10,  9, 13, 12, 14, 23],
        [10, 13, 16, 23, 18, 20],
        [16,  6,  5, 19,  8, 24],
        [16,  5, 10, 24, 11, 20]
    ]
}

fn get_coordinates_block() -> NodalCoordinates<D>
{
    todo!()
}

fn get_reference_coordinates() -> ReferenceNodalCoordinates<N>
{
    ReferenceNodalCoordinates::new([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.5, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.0, 0.5, 0.0]
    ])
}

fn get_reference_coordinates_block() -> ReferenceNodalCoordinates<D>
{
    ReferenceNodalCoordinates::new([
        [ 0.50,  0.50,  0.00],
        [-0.50,  0.50,  0.00],
        [ 0.00,  0.50,  0.00],
        [ 0.25,  0.50,  0.00],
        [-0.25,  0.50,  0.00],
        [-0.50, -0.50,  0.00],
        [-0.50,  0.00,  0.00],
        [-0.50,  0.25,  0.00],
        [-0.50, -0.25,  0.00],
        [ 0.50, -0.50,  0.00],
        [ 0.00, -0.50,  0.00],
        [-0.25, -0.50,  0.00],
        [ 0.25, -0.50,  0.00],
        [ 0.50,  0.00,  0.00],
        [ 0.50, -0.25,  0.00],
        [ 0.50,  0.25,  0.00],
        [ 0.00,  0.00,  0.00],
        [ 0.00,  0.25,  0.00],
        [ 0.25,  0.00,  0.00],
        [-0.25,  0.00,  0.00],
        [ 0.00, -0.25,  0.00],
        [ 0.25,  0.25,  0.00],
        [-0.25,  0.25,  0.00],
        [ 0.25, -0.25,  0.00],
        [-0.25, -0.25,  0.00]
    ])
}

fn get_velocities_block() -> NodalVelocities<D>
{
    todo!()
}

test_composite_surface_element!(Triangle);
test_finite_element_block!(Triangle);