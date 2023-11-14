use crate::fem::block::
{
    test::test_finite_element_block,
    element::
    {
        test::test_finite_element,
        linear::
        {
            test::test_linear_finite_element
        }
    }
};
use super::*;

const D: usize = 14;
const E: usize = 24;

fn get_connectivity() -> Connectivity<E, N>
{
    [
        [13, 12,  8,  1],
        [10,  3,  0,  8],
        [11, 10,  8,  3],
        [12, 11,  8,  2],
        [11,  2,  3,  8],
        [12,  2,  8,  1],
        [13, 10,  5,  0],
        [13, 11, 10,  8],
        [10,  6,  9,  5],
        [12,  7,  4,  9],
        [12, 11,  7,  9],
        [11,  7,  9,  6],
        [13,  1,  8,  0],
        [13,  9,  4,  5],
        [13, 12,  1,  4],
        [11, 10,  6,  9],
        [11, 10,  3,  6],
        [12, 11,  2,  7],
        [13, 11,  9, 10],
        [13, 12,  4,  9],
        [13, 10,  0,  8],
        [13, 10,  9,  5],
        [13, 12, 11,  8],
        [13, 12,  9, 11]
    ]
}

fn get_current_coordinates_block() -> CurrentNodalCoordinates<D>
{
    CurrentNodalCoordinates::new([
       [ 0.48419081, -0.52698494,  0.42026988],
       [ 0.43559430,  0.52696224,  0.54477963],
       [-0.56594965,  0.57076191,  0.51683869],
       [-0.56061746, -0.42795457,  0.55275658],
       [ 0.41878700,  0.53190268, -0.44744274],
       [ 0.47232357, -0.57252738, -0.42946606],
       [-0.45168197, -0.5102938 , -0.57959825],
       [-0.41776733,  0.41581785, -0.45911886],
       [ 0.05946988,  0.03773822,  0.44149305],
       [-0.08478334, -0.09009810, -0.46105872],
       [-0.04039882, -0.58201398,  0.09346960],
       [-0.57820738,  0.08325131,  0.03614415],
       [-0.04145077,  0.56406301,  0.09988905],
       [ 0.52149656, -0.08553510, -0.03187069]
    ])
}

fn get_reference_coordinates() -> ReferenceNodalCoordinates<N>
{
    ReferenceNodalCoordinates::new([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0]
    ])
}

fn get_reference_coordinates_block() -> ReferenceNodalCoordinates<D>
{
    ReferenceNodalCoordinates::new([
        [ 0.5, -0.5,  0.5],
        [ 0.5,  0.5,  0.5],
        [-0.5,  0.5,  0.5],
        [-0.5, -0.5,  0.5],
        [ 0.5,  0.5, -0.5],
        [ 0.5, -0.5, -0.5],
        [-0.5, -0.5, -0.5],
        [-0.5,  0.5, -0.5],
        [ 0.0,  0.0,  0.5],
        [ 0.0,  0.0, -0.5],
        [ 0.0, -0.5,  0.0],
        [-0.5,  0.0,  0.0],
        [ 0.0,  0.5,  0.0],
        [ 0.5,  0.0,  0.0]
    ])
}

test_finite_element!(LinearTetrahedron);
test_finite_element_block!(LinearTetrahedron);
test_linear_finite_element!(LinearTetrahedron);
