use super::*;
use crate::fem::block::{
    element::{
        composite::surface::test::{
            setup_for_test_composite_element_with_constitutive_model,
            setup_for_test_composite_surface_element_with_constitutive_model,
            test_composite_surface_element,
        },
        linear::surface::test::{
            setup_for_test_finite_element_block_with_elastic_constitutive_model,
            setup_for_test_finite_element_with_elastic_constitutive_model,
        },
    },
    test::test_finite_element_block,
};

const D: usize = 25;
const E: usize = 8;

fn get_connectivity() -> Connectivity<E, N> {
    [
        [0, 2, 16, 3, 17, 21],
        [0, 16, 13, 21, 18, 15],
        [2, 1, 6, 4, 7, 22],
        [2, 6, 16, 22, 19, 17],
        [10, 9, 13, 12, 14, 23],
        [10, 13, 16, 23, 18, 20],
        [16, 6, 5, 19, 8, 24],
        [16, 5, 10, 24, 11, 20],
    ]
}

fn get_coordinates_block() -> NodalCoordinates<D> {
    NodalCoordinates::new([
        [0.49663681, 0.48487177, 0.00470544],
        [-0.51103016, 0.49201974, 0.00373833],
        [0.00955038, 0.50231237, 0.00266678],
        [0.26843046, 0.51634742, -0.00528073],
        [-0.25992891, 0.49576645, 0.01744170],
        [-0.51783719, -0.48436913, 0.00253641],
        [-0.49364474, -0.00190810, 0.00929366],
        [-0.51135401, 0.24340846, 0.01328358],
        [-0.48397380, -0.26270101, 0.01084105],
        [0.51356010, -0.48228234, 0.01219576],
        [0.00858244, -0.51827679, -0.01784131],
        [-0.26890474, -0.51830564, 0.00891536],
        [0.26901624, -0.51941609, 0.00380471],
        [0.50702926, -0.00537950, -0.01900236],
        [0.49003802, -0.23435359, 0.01796925],
        [0.50573890, 0.23902918, 0.01252039],
        [0.00245407, -0.01773899, -0.01080003],
        [-0.01484150, 0.23000923, 0.00186402],
        [0.25434975, 0.01811168, 0.00228300],
        [-0.24457155, -0.01450823, 0.01306671],
        [-0.01233367, -0.24302578, -0.01923283],
        [0.24975239, 0.24588480, 0.01715625],
        [-0.24681722, 0.26888657, -0.01611686],
        [0.23716171, -0.24092213, -0.01764789],
        [-0.23431565, -0.26743560, -0.01928865],
    ])
}

fn get_reference_coordinates() -> ReferenceNodalCoordinates<N> {
    ReferenceNodalCoordinates::new([
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.0, 0.5, 0.0],
        [0.5, 0.0, 0.0],
    ])
}

fn get_reference_coordinates_block() -> ReferenceNodalCoordinates<D> {
    ReferenceNodalCoordinates::new([
        [0.50, 0.50, 0.00],
        [-0.50, 0.50, 0.00],
        [0.00, 0.50, 0.00],
        [0.25, 0.50, 0.00],
        [-0.25, 0.50, 0.00],
        [-0.50, -0.50, 0.00],
        [-0.50, 0.00, 0.00],
        [-0.50, 0.25, 0.00],
        [-0.50, -0.25, 0.00],
        [0.50, -0.50, 0.00],
        [0.00, -0.50, 0.00],
        [-0.25, -0.50, 0.00],
        [0.25, -0.50, 0.00],
        [0.50, 0.00, 0.00],
        [0.50, -0.25, 0.00],
        [0.50, 0.25, 0.00],
        [0.00, 0.00, 0.00],
        [0.00, 0.25, 0.00],
        [0.25, 0.00, 0.00],
        [-0.25, 0.00, 0.00],
        [0.00, -0.25, 0.00],
        [0.25, 0.25, 0.00],
        [-0.25, 0.25, 0.00],
        [0.25, -0.25, 0.00],
        [-0.25, -0.25, 0.00],
    ])
}

fn get_velocities_block() -> NodalVelocities<D> {
    NodalVelocities::new([
        [0.03898679, 0.01794692, 0.08647724],
        [0.01143928, 0.01148802, 0.09322342],
        [0.02331058, 0.02625677, 0.07657580],
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
        [0.07692136, 0.07541025, 0.07702821],
    ])
}

test_composite_surface_element!(Triangle);
test_finite_element_block!(Triangle);
