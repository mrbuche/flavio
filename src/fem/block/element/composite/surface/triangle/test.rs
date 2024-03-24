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
    NodalCoordinates::new([
        [ 0.51258167,  0.47870565, -0.01444634],
        [-0.47661871,  0.51241734,  0.02185399],
        [-0.03854309,  0.49586239, -0.02292243],
        [ 0.26983095,  0.51164101,  0.01936741],
        [-0.20585152,  0.53887608,  0.04927917],
        [-0.50576668, -0.50332372, -0.00952088],
        [-0.46094531, -0.01265329,  0.04323734],
        [-0.45062844,  0.24744830, -0.03467063],
        [-0.53795563, -0.27080582, -0.04915048],
        [ 0.50746161, -0.45226548, -0.04821875],
        [ 0.04403169, -0.51287674,  0.00866001],
        [-0.26720914, -0.53368917,  0.02890714],
        [ 0.26024402, -0.51554662, -0.00472360],
        [ 0.52483278,  0.00068077, -0.04668617],
        [ 0.49658227, -0.28263177,  0.04286684],
        [ 0.45599904,  0.25718471, -0.02772888],
        [-0.03687598,  0.00097344,  0.01425001],
        [-0.04034087,  0.26015899, -0.01316944],
        [ 0.25718593,  0.03223932,  0.03989718],
        [-0.28504236, -0.03604217, -0.01963145],
        [ 0.01291068, -0.25169747, -0.04760444],
        [ 0.24822420,  0.22998421, -0.01714032],
        [-0.25044620,  0.26047683,  0.01504145],
        [ 0.28576016, -0.20234698, -0.04457882],
        [-0.24189247, -0.22695437, -0.04827701]
    ])
}

fn get_reference_coordinates() -> ReferenceNodalCoordinates<N>
{
    ReferenceNodalCoordinates::new([
        // [0.0, 0.0, 0.0],
        // [1.0, 0.0, 0.0],
        // [0.0, 1.0, 0.0],
        // [0.5, 0.0, 0.0],
        // [0.5, 0.5, 0.0],
        // [0.0, 0.5, 0.0]
        // [1.0, 0.0, 0.0],
        // [0.0, 1.0, 0.0],
        // [0.0, 0.0, 0.0],
        // [0.5, 0.5, 0.0],
        // [0.0, 0.5, 0.0],
        // [0.5, 0.0, 0.0]
        [0.97, 0.02, 0.01],
        [-0.05, 1.03, 0.05],
        [0.0, 0.0,-0.001],
        [0.53, 0.5, 0.0],
        [-0.05, 0.5, 0.05],
        [0.48, 0.0, -0.05]
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
    NodalVelocities::new([
        [0.03898679, 0.01794692, 0.08647724],
        [0.01143928, 0.01148802, 0.09322342],
        [0.02331058, 0.02625677, 0.0765758 ],
        [0.01930480, 0.03435541, 0.07127602],
        [0.00690091, 0.02234559, 0.04352873],
        [0.09760132, 0.02284640, 0.00156334],
        [0.07934113, 0.07350421, 0.00786410],
        [0.00501462, 0.05750917, 0.03766337],
        [0.07037483, 0.04200882, 0.04983480],
        [0.04192370, 0.02365160, 0.08386066],
        [0.03517049, 0.07349123, 0.07568995],
        [0.04616277, 0.01503475, 0.04150365],
        [0.06207289, 0.01097912, 0.01624590],
        [0.05584760, 0.09164884, 0.01524798],
        [0.05475846, 0.02749090, 0.01587344],
        [0.03600224, 0.07138283, 0.04143816],
        [0.09025223, 0.09655957, 0.09073203],
        [0.08663481, 0.05015400, 0.00432735],
        [0.08613135, 0.02408223, 0.04805299],
        [0.05386770, 0.07591057, 0.02818401],
        [0.09307117, 0.06591787, 0.06368264],
        [0.04242519, 0.06270905, 0.04721817],
        [0.03438826, 0.09398549, 0.05947200],
        [0.02522213, 0.04261638, 0.02932986],
        [0.07692136, 0.07541025, 0.07702821]
    ])
}

test_composite_surface_element!(Triangle);
test_finite_element_block!(Triangle);