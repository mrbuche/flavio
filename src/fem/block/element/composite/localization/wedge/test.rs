use crate::fem::block::
{
    element::composite::localization::test::
    {
        setup_for_test_composite_element_with_constitutive_model,
        setup_for_test_composite_surface_element_with_constitutive_model,
        test_composite_localization_element
    },
    test::test_finite_element_block
};
use super::*;

const D: usize = 30;
const E: usize = 4;

fn get_connectivity() -> Connectivity<E, N>
{
    [
        [1, 5, 4, 16, 20, 19, 6,  7,  9, 21, 22, 24],
        [4, 3, 0, 19, 18, 15, 8, 10, 11, 23, 25, 26],
        [5, 1, 2, 20, 16, 17, 6, 12, 13, 21, 27, 28],
        [3, 4, 5, 18, 19, 20, 8,  7, 14, 23, 22, 29]
    ]
}

fn get_coordinates_block() -> NodalCoordinates<D>
{
    NodalCoordinates::new([
        [ 0.50255857,  0.24347417, -0.01615211],
        [-0.49942154,  0.26094478,  0.00416315],
        [-0.50759508, -0.25724750, -0.00599802],
        [ 0.49289393, -0.23468213,  0.00699716],
        [ 0.01579151,  0.25376321, -0.01834600],
        [-0.00883798, -0.23920464, -0.01280460],
        [-0.26988463, -0.00257442, -0.00873978],
        [-0.00771332,  0.01914662,  0.00953987],
        [ 0.26484067, -0.01948331,  0.01835137],
        [-0.25308186,  0.25771941, -0.01380911],
        [ 0.50834657, -0.01583413,  0.01780817],
        [ 0.23008039,  0.25421365,  0.00361367],
        [-0.51641813,  0.00163187,  0.00276274],
        [-0.23039358, -0.24000316, -0.01582590],
        [ 0.26513732, -0.23561404, -0.00743661],
        [ 0.50842270,  0.26881443,  0.01983390],
        [-0.48069686,  0.25571890, -0.01976750],
        [-0.51436448, -0.24251337, -0.01184794],
        [ 0.51641935, -0.25062896,  0.01178048],
        [ 0.00617594,  0.26697957, -0.00101055],
        [ 0.01663657, -0.26921294, -0.00171786],
        [-0.23589937,  0.00107257,  0.01023294],
        [-0.00685220,  0.01063963,  0.00358080],
        [ 0.23680031,  0.01264953,  0.01477911],
        [-0.24138476,  0.25920237, -0.01562145],
        [ 0.49751986,  0.01590411, -0.00847248],
        [ 0.23183119,  0.26866829, -0.00683692],
        [-0.49078100, -0.01555533,  0.00133412],
        [-0.26674660, -0.24334801, -0.01123147],
        [ 0.25780943, -0.23445127,  0.01898350]
    ])
}

fn get_reference_coordinates() -> ReferenceNodalCoordinates<N>
{
    ReferenceNodalCoordinates::new([
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.0, 0.5, 0.0],
        [0.5, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.0, 0.5, 0.0],
        [0.5, 0.0, 0.0]
    ])
}

fn get_reference_coordinates_block() -> ReferenceNodalCoordinates<D>
{
    ReferenceNodalCoordinates::new([
        [ 0.50,  0.25, 0.00],
        [-0.50,  0.25, 0.00],
        [-0.50, -0.25, 0.00],
        [ 0.50, -0.25, 0.00],
        [ 0.00,  0.25, 0.00],
        [ 0.00, -0.25, 0.00],
        [-0.25,  0.00, 0.00],
        [ 0.00,  0.00, 0.00],
        [ 0.25,  0.00, 0.00],
        [-0.25,  0.25, 0.00],
        [ 0.50,  0.00, 0.00],
        [ 0.25,  0.25, 0.00],
        [-0.50,  0.00, 0.00],
        [-0.25, -0.25, 0.00],
        [ 0.25, -0.25, 0.00],
        [ 0.50,  0.25, 0.00],
        [-0.50,  0.25, 0.00],
        [-0.50, -0.25, 0.00],
        [ 0.50, -0.25, 0.00],
        [ 0.00,  0.25, 0.00],
        [ 0.00, -0.25, 0.00],
        [-0.25,  0.00, 0.00],
        [ 0.00,  0.00, 0.00],
        [ 0.25,  0.00, 0.00],
        [-0.25,  0.25, 0.00],
        [ 0.50,  0.00, 0.00],
        [ 0.25,  0.25, 0.00],
        [-0.50,  0.00, 0.00],
        [-0.25, -0.25, 0.00],
        [ 0.25, -0.25, 0.00]
    ])
}

fn get_velocities_block() -> NodalVelocities<D>
{
    NodalVelocities::new([
        [0.01431412, 0.01804358, 0.07496300],
        [0.09962419, 0.04549100, 0.02328813],
        [0.07403234, 0.03900588, 0.02604932],
        [0.08251378, 0.00446619, 0.06083218],
        [0.07505936, 0.06427368, 0.00360294],
        [0.03396493, 0.05562344, 0.08084064],
        [0.03920721, 0.09974366, 0.00835274],
        [0.06322438, 0.04215104, 0.07457214],
        [0.00214199, 0.01956594, 0.09590513],
        [0.00158532, 0.05293008, 0.04495069],
        [0.02439420, 0.05474652, 0.06894332],
        [0.04635582, 0.02833468, 0.05637430],
        [0.08323779, 0.09134513, 0.00270436],
        [0.04526144, 0.09675598, 0.07204048],
        [0.09081726, 0.05236692, 0.06802666],
        [0.08887790, 0.04621419, 0.01521392],
        [0.02594481, 0.03407545, 0.09746018],
        [0.02060314, 0.03239568, 0.09716370],
        [0.00317384, 0.01425607, 0.09150469],
        [0.03201858, 0.09884590, 0.04427266],
        [0.09162680, 0.02950793, 0.06983756],
        [0.03955233, 0.00967575, 0.00726536],
        [0.00362971, 0.04974553, 0.03535663],
        [0.04301489, 0.05819536, 0.07788356],
        [0.02063973, 0.03599813, 0.00392894],
        [0.04542459, 0.03250225, 0.00345714],
        [0.01380998, 0.08523858, 0.00614054],
        [0.00438224, 0.02015492, 0.03752064],
        [0.09450938, 0.07154234, 0.02198726],
        [0.09789419, 0.04941774, 0.02401724]
    ])
}

test_composite_localization_element!(Wedge);
test_finite_element_block!(Wedge);