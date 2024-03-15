use crate::fem::block::element::linear::surface::test::
{
    setup_for_test_linear_surface_element_with_constitutive_model,
    test_linear_surface_element
};
use super::*;

const D: usize = 16;
const E: usize = 18;

fn get_connectivity() -> Connectivity<E, N>
{
    [
        [ 6,  4,  8],
        [ 8, 15,  6],
        [15,  8,  9],
        [ 9, 14, 15],
        [14,  9,  7],
        [ 7, 10, 14],
        [ 5,  6, 15],
        [15, 13,  5],
        [13, 15, 14],
        [14, 12, 13],
        [12, 14, 10],
        [10, 11, 12],
        [ 1,  5, 13],
        [13,  3,  1],
        [ 3, 13, 12],
        [12,  2,  3],
        [ 2, 12, 11],
        [11,  0,  2]
    ]
}

fn get_reference_coordinates() -> ReferenceNodalCoordinates<N>
{
    ReferenceNodalCoordinates::new([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0]
    ])
}

fn get_reference_coordinates_block() -> ReferenceNodalCoordinates<D>
{
    ReferenceNodalCoordinates::new([
        [     0.5,      0.5, 0.0],
        [    -0.5,      0.5, 0.0],
        [ 1.0/6.0,      0.5, 0.0],
        [-1.0/6.0,      0.5, 0.0],
        [    -0.5,     -0.5, 0.0],
        [    -0.5,  1.0/6.0, 0.0],
        [    -0.5, -1.0/6.0, 0.0],
        [     0.5,     -0.5, 0.0],
        [-1.0/6.0,     -0.5, 0.0],
        [ 1.0/6.0,     -0.5, 0.0],
        [     0.5, -1.0/6.0, 0.0],
        [     0.5,  1.0/6.0, 0.0],
        [ 1.0/6.0,  1.0/6.0, 0.0],
        [-1.0/6.0,  1.0/6.0, 0.0],
        [ 1.0/6.0, -1.0/6.0, 0.0],
        [-1.0/6.0, -1.0/6.0, 0.0]
    ])
}

test_linear_surface_element!(Triangle);

#[test]
fn do_block()
{
    todo!()
}
