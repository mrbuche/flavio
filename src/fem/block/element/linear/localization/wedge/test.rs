use super::*;
use crate::fem::block::{
    element::linear::{
        localization::test::{
            setup_for_test_linear_surface_element_with_constitutive_model,
            test_linear_localization_element,
        },
        surface::test::{
            setup_for_test_finite_element_block_with_elastic_constitutive_model,
            setup_for_test_finite_element_with_elastic_constitutive_model,
            setup_for_test_linear_element_with_constitutive_model,
        },
    },
    test::test_finite_element_block,
};

const D: usize = 32;
const E: usize = 18;

fn get_connectivity() -> Connectivity<E, N> {
    [
        [6, 4, 8, 22, 20, 24],
        [8, 15, 6, 24, 31, 22],
        [15, 8, 9, 31, 24, 25],
        [9, 14, 15, 25, 30, 31],
        [14, 9, 7, 30, 25, 23],
        [7, 10, 14, 23, 26, 30],
        [5, 6, 15, 21, 22, 31],
        [15, 13, 5, 31, 29, 21],
        [13, 15, 14, 29, 31, 30],
        [14, 12, 13, 30, 28, 29],
        [12, 14, 10, 28, 30, 26],
        [10, 11, 12, 26, 27, 28],
        [1, 5, 13, 17, 21, 29],
        [13, 3, 1, 29, 19, 17],
        [3, 13, 12, 19, 29, 28],
        [12, 2, 3, 28, 18, 19],
        [2, 12, 11, 18, 28, 27],
        [11, 0, 2, 27, 16, 18],
    ]
}

fn get_coordinates_block() -> NodalCoordinates<D> {
    NodalCoordinates::new([
        [0.52493242, 0.46007179, -0.01313593],
        [-0.49249732, 0.45861313, 0.00692950],
        [0.13681751, 0.49739797, -0.03913867],
        [-0.16906010, 0.45570952, 0.04155203],
        [-0.50784484, -0.51557480, 0.02400643],
        [-0.52697789, 0.15806322, -0.02101903],
        [-0.50460505, -0.20384408, 0.04339776],
        [0.50284742, -0.52930851, 0.01780512],
        [-0.12117202, -0.47170986, -0.01554877],
        [0.17068755, -0.48362566, 0.04805858],
        [0.49315351, -0.21084858, 0.00677581],
        [0.51913636, 0.18622272, -0.02365121],
        [0.13880631, 0.18701065, 0.02217795],
        [-0.17167971, 0.12074383, 0.01228280],
        [0.18742507, -0.18877085, -0.03995387],
        [-0.20089730, -0.17737422, 0.03563169],
        [0.53224136, 0.47503067, 0.96814942],
        [-0.52879305, 0.48946303, 1.01865459],
        [0.17381915, 0.48720033, 0.96261128],
        [-0.13847523, 0.46435300, 0.99918085],
        [-0.50177781, -0.45903988, 0.95547808],
        [-0.52394573, 0.11984925, 1.03092470],
        [-0.47797934, -0.16358887, 0.99753093],
        [0.46629365, -0.53002434, 0.96840737],
        [-0.20840967, -0.47733345, 0.97379936],
        [0.12606825, -0.47465160, 0.97819479],
        [0.49362273, -0.21543374, 0.96389391],
        [0.47883486, 0.17160061, 1.03083090],
        [0.15057969, 0.17550796, 1.00824486],
        [-0.13764710, 0.17932179, 1.00134590],
        [0.18300281, -0.17029714, 0.97287949],
        [-0.14469107, -0.12224519, 0.99747609],
    ])
}

fn get_reference_coordinates() -> ReferenceNodalCoordinates<N> {
    ReferenceNodalCoordinates::new([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    ])
}

fn get_reference_coordinates_block() -> ReferenceNodalCoordinates<D> {
    ReferenceNodalCoordinates::new([
        [0.5, 0.5, 0.0],
        [-0.5, 0.5, 0.0],
        [1.0 / 6.0, 0.5, 0.0],
        [-1.0 / 6.0, 0.5, 0.0],
        [-0.5, -0.5, 0.0],
        [-0.5, 1.0 / 6.0, 0.0],
        [-0.5, -1.0 / 6.0, 0.0],
        [0.5, -0.5, 0.0],
        [-1.0 / 6.0, -0.5, 0.0],
        [1.0 / 6.0, -0.5, 0.0],
        [0.5, -1.0 / 6.0, 0.0],
        [0.5, 1.0 / 6.0, 0.0],
        [1.0 / 6.0, 1.0 / 6.0, 0.0],
        [-1.0 / 6.0, 1.0 / 6.0, 0.0],
        [1.0 / 6.0, -1.0 / 6.0, 0.0],
        [-1.0 / 6.0, -1.0 / 6.0, 0.0],
        [0.5, 0.5, 0.0],
        [-0.5, 0.5, 0.0],
        [1.0 / 6.0, 0.5, 0.0],
        [-1.0 / 6.0, 0.5, 0.0],
        [-0.5, -0.5, 0.0],
        [-0.5, 1.0 / 6.0, 0.0],
        [-0.5, -1.0 / 6.0, 0.0],
        [0.5, -0.5, 0.0],
        [-1.0 / 6.0, -0.5, 0.0],
        [1.0 / 6.0, -0.5, 0.0],
        [0.5, -1.0 / 6.0, 0.0],
        [0.5, 1.0 / 6.0, 0.0],
        [1.0 / 6.0, 1.0 / 6.0, 0.0],
        [-1.0 / 6.0, 1.0 / 6.0, 0.0],
        [1.0 / 6.0, -1.0 / 6.0, 0.0],
        [-1.0 / 6.0, -1.0 / 6.0, 0.0],
    ])
}

fn get_velocities_block() -> NodalVelocities<D> {
    NodalCoordinates::new([
        [-0.06433848, -0.01926187, -0.03743756],
        [-0.0870451, -0.04597518, -0.03565326],
        [0.04400981, -0.09000301, 0.04678533],
        [0.00173760, -0.08503324, -0.07858199],
        [-0.05268299, -0.04592905, 0.00647059],
        [-0.03494861, -0.00688540, -0.07195594],
        [0.03256779, -0.05282913, -0.08925343],
        [-0.07227557, 0.01213563, 0.00693205],
        [0.08459964, 0.05677483, 0.01987572],
        [-0.02275489, -0.07373460, -0.08081551],
        [-0.03142081, -0.08566530, 0.06210498],
        [-0.06069344, 0.03923413, 0.06299218],
        [0.06499705, -0.07443175, -0.00382237],
        [-0.06030412, -0.08955708, -0.04993707],
        [0.02621134, -0.06516815, -0.01265165],
        [-0.03223677, 0.01629867, 0.08890064],
        [-0.00062922, 0.09861261, -0.02108821],
        [0.00093658, 0.01677123, 0.07564730],
        [-0.03366722, -0.00382467, -0.01066520],
        [0.03326620, -0.05894294, -0.09974025],
        [-0.08501372, -0.00149077, 0.05720255],
        [-0.04737467, 0.05557303, 0.03014414],
        [0.04400757, 0.02143862, 0.08174897],
        [-0.09654050, -0.01515334, -0.02669490],
        [0.04375730, -0.09169403, 0.07701589],
        [0.08591825, -0.00565911, 0.02977132],
        [0.05940797, 0.04279079, -0.01290276],
        [-0.05899700, 0.01084053, -0.07248295],
        [0.00876680, -0.06769134, 0.00247064],
        [-0.03368043, -0.07116854, 0.09420690],
        [-0.01459746, -0.06751504, -0.07558666],
        [-0.08412322, -0.09261289, -0.04835088],
    ])
}

fn get_dirichlet_places<'a>() -> [&'a [usize]; 10] {
    panic!()
}

fn get_dirichlet_values(_x: Scalar) -> [Scalar; 10] {
    panic!()
}

test_linear_localization_element!(Wedge);
test_finite_element_block!(Wedge);
