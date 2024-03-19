use crate::fem::block::
{
    element::composite::test::test_composite_element,
    test::test_finite_element_block
};
use super::*;

const D: usize = 63;
const E: usize = 24;

fn get_connectivity() -> Connectivity<E, N>
{
    [
        [13, 12,  8,  1, 50, 51, 52, 29, 43, 31],
        [10,  3,  0,  8, 47, 22, 49, 53, 32, 30],
        [11, 10,  8,  3, 54, 53, 55, 40, 47, 32],
        [12, 11,  8,  2, 56, 55, 51, 44, 39, 33],
        [11,  2,  3,  8, 39, 15, 40, 55, 33, 32],
        [12,  2,  8,  1, 44, 33, 51, 43, 21, 31],
        [13, 10,  5,  0, 57, 48, 26, 27, 49, 14],
        [13, 11, 10,  8, 58, 54, 57, 52, 55, 53],
        [10,  6,  9,  5, 46, 36, 59, 48, 17, 34],
        [12,  7,  4,  9, 45, 20, 42, 60, 37, 35],
        [12, 11,  7,  9, 56, 38, 45, 60, 61, 37],
        [11,  7,  9,  6, 38, 37, 61, 41, 19, 36],
        [13,  1,  8,  0, 29, 31, 52, 27, 25, 30],
        [13,  9,  4,  5, 62, 35, 28, 26, 34, 23],
        [13, 12,  1,  4, 50, 43, 29, 28, 42, 18],
        [11, 10,  6,  9, 54, 46, 41, 61, 59, 36],
        [11, 10,  3,  6, 54, 47, 40, 41, 46, 16],
        [12, 11,  2,  7, 56, 39, 44, 45, 38, 24],
        [13, 11,  9, 10, 58, 61, 62, 57, 54, 59],
        [13, 12,  4,  9, 50, 42, 28, 62, 60, 35],
        [13, 10,  0,  8, 57, 49, 27, 52, 53, 30],
        [13, 10,  9,  5, 57, 59, 62, 26, 48, 34],
        [13, 12, 11,  8, 50, 56, 58, 52, 51, 55],
        [13, 12,  9, 11, 50, 60, 62, 58, 56, 61]
    ]
}

fn get_reference_coordinates_block() -> ReferenceNodalCoordinates<D>
{
    ReferenceNodalCoordinates::new([
        [ 0.50, -0.50,  0.50],
        [ 0.50,  0.50,  0.50],
        [-0.50,  0.50,  0.50],
        [-0.50, -0.50,  0.50],
        [ 0.50,  0.50, -0.50],
        [ 0.50, -0.50, -0.50],
        [-0.50, -0.50, -0.50],
        [-0.50,  0.50, -0.50],
        [ 0.00,  0.00,  0.50],
        [ 0.00,  0.00, -0.50],
        [ 0.00, -0.50,  0.00],
        [-0.50,  0.00,  0.00],
        [ 0.00,  0.50,  0.00],
        [ 0.50,  0.00,  0.00],
        [ 0.25, -0.25,  0.50],
        [ 0.25,  0.25,  0.50],
        [-0.25, -0.25,  0.50],
        [-0.25,  0.25,  0.50],
        [ 0.50,  0.00,  0.50],
        [-0.50,  0.00,  0.50],
        [ 0.00, -0.50,  0.50],
        [ 0.00,  0.50,  0.50],
        [ 0.25, -0.25, -0.50],
        [ 0.25,  0.25, -0.50],
        [-0.25, -0.25, -0.50],
        [-0.25,  0.25, -0.50],
        [ 0.50,  0.00, -0.50],
        [ 0.00, -0.50, -0.50],
        [-0.50,  0.00, -0.50],
        [ 0.00,  0.50, -0.50],
        [-0.50,  0.25, -0.25],
        [-0.50,  0.25,  0.25],
        [-0.50, -0.25,  0.25],
        [-0.50, -0.25, -0.25],
        [-0.50,  0.50,  0.00],
        [-0.50, -0.50,  0.00],
        [ 0.25,  0.50, -0.25],
        [ 0.25,  0.50,  0.25],
        [-0.25,  0.50,  0.25],
        [-0.25,  0.50, -0.25],
        [ 0.50,  0.50,  0.00],
        [-0.25, -0.50, -0.25],
        [-0.25, -0.50,  0.25],
        [ 0.25, -0.50, -0.25],
        [ 0.25, -0.50,  0.25],
        [ 0.50, -0.50,  0.00],
        [ 0.50, -0.25, -0.25],
        [ 0.50, -0.25,  0.25],
        [ 0.50,  0.25, -0.25],
        [ 0.50,  0.25,  0.25],
        [ 0.25,  0.25,  0.00],
        [ 0.00,  0.25,  0.25],
        [ 0.25,  0.25,  0.25],
        [ 0.00, -0.25,  0.25],
        [-0.25, -0.25,  0.00],
        [-0.25, -0.25,  0.25],
        [-0.25,  0.25,  0.00],
        [ 0.25, -0.25,  0.00],
        [ 0.00,  0.00,  0.00],
        [ 0.00, -0.25, -0.25],
        [ 0.00,  0.25, -0.25],
        [-0.25, -0.25, -0.25],
        [ 0.25,  0.25, -0.25]
    ])
}

test_composite_element!(Tetrahedron);
test_finite_element_block!(Tetrahedron);





use crate::constitutive::solid::hyperelastic::
{
    NeoHookean,
    test::NEOHOOKEANPARAMETERS
};
use crate::
{
    EPSILON,
    test::assert_eq_within_tols
};

fn get_reference_coordinates() -> ReferenceNodalCoordinates<N>
{
    ReferenceNodalCoordinates::new([
        [ 0.00,  0.00,  0.00],
        [1.0,  0.00,  0.00],
        [ 0.00, 1.0,  0.00],
        [ 0.00,  0.00, 1.0],
        [ 0.50,  0.00,  0.00],
        [ 0.50,  0.50,  0.00],
        [ 0.00,  0.50,  0.00],
        [ 0.00,  0.00,  0.50],
        [ 0.50,  0.00,  0.50],
        [ 0.00,  0.50,  0.50]
    ])
}

fn get_element<'a>() -> Tetrahedron::<NeoHookean<'a>>
{
    Tetrahedron::new(
        NEOHOOKEANPARAMETERS, get_reference_coordinates()
    )
}

fn get_coordinatess() -> NodalCoordinates<N>
{
    get_deformation_gradient() * get_reference_coordinates()
}

fn get_reference_coordinates_crazy() -> ReferenceNodalCoordinates<N>
{
    ReferenceNodalCoordinates::new([
        [ 0.87900980,  0.73292419,  0.67130546],
        [ 0.89557721,  0.40375233,  0.94588738],
        [ 0.98513155,  0.501115036,  0.502239128],
        [ 0.42737582,  0.87092470,  0.20908760],
        [ 0.78803900,  0.83651295,  0.69155253],
        [ 0.23296807,  0.41294924,  0.26352866],
        [ 0.26965179,  0.63112838,  0.008563159],
        [ 0.19258050,  0.508921056,  0.87556116],
        [ 0.45519136,  0.90772602,  0.26139220],
        [ 0.47893746,  0.64353023,  0.004217477]
    ])
}

fn get_element_crazy<'a>() -> Tetrahedron::<NeoHookean<'a>>
{
    Tetrahedron::new(
        NEOHOOKEANPARAMETERS, get_reference_coordinates_crazy()
    )
}

fn get_coordinates_crazy() -> NodalCoordinates<N>
{
    get_deformation_gradient() * get_reference_coordinates_crazy()
}

#[test]
fn partition_of_unity<'a>()
{
    Tetrahedron::<NeoHookean<'a>>::calculate_standard_gradient_operators()
    .iter()
    .for_each(|standard_gradient_operator|{
        let mut sum = [ 0.00_f64; 3];
        standard_gradient_operator.iter()
        .for_each(|row|
            row.iter()
            .zip(sum.iter_mut())
            .for_each(|(entry, sum_i)|
                *sum_i += entry
            )
        );
        sum.iter()
        .for_each(|sum_i|
            assert_eq_within_tols(sum_i, &0.00)
        )
    })
}

#[test]
fn undeformed_crazy()
{
    get_element_crazy().calculate_deformation_gradients(
        &get_reference_coordinates_crazy().convert()
    ).iter()
    .for_each(|deformation_gradient_g|
        deformation_gradient_g.iter()
        .enumerate()
        .for_each(|(i, deformation_gradient_g_i)|
            deformation_gradient_g_i.iter()
            .enumerate()
            .for_each(|(j, deformation_gradient_g_ij)|
                assert_eq_within_tols(
                    deformation_gradient_g_ij, &((i == j) as u8 as Scalar)
                )
            )
        )
    )
}

#[test]
fn deformed()
{
    get_element().calculate_deformation_gradients(
        &get_coordinatess()
    ).iter()
    .for_each(|deformation_gradient_g|
        deformation_gradient_g.iter()
        .zip(get_deformation_gradient().iter())
        .for_each(|(deformation_gradient_g_i, get_deformation_gradient_i)|
            deformation_gradient_g_i.iter()
            .zip(get_deformation_gradient_i.iter())
            .for_each(|(deformation_gradient_g_ij, get_deformation_gradient_ij)|
                assert_eq_within_tols(
                    deformation_gradient_g_ij, get_deformation_gradient_ij
                )
            )
        )
    )
}

fn get_finite_difference_of_helmholtz_free_energy() -> NodalForces<N>
{
    let element = get_element();
    let mut finite_difference = 0.00;
    (0..N).map(|node|
        (0..3).map(|i|{
            let mut nodal_coordinates = get_coordinatess();
            nodal_coordinates[node][i] += 0.50 * EPSILON;
            finite_difference = element.calculate_helmholtz_free_energy(
                &nodal_coordinates
            );
            nodal_coordinates[node][i] -= EPSILON;
            finite_difference -= element.calculate_helmholtz_free_energy(
                &nodal_coordinates
            );
            finite_difference/EPSILON
        }).collect()
    ).collect()
}

#[test]
fn finite_difference()
{
    get_element().calculate_nodal_forces(
        &get_coordinatess()
    ).iter()
    .zip(get_finite_difference_of_helmholtz_free_energy().iter())
    .for_each(|(nodal_force, fd_nodal_force)|
        nodal_force.iter()
        .zip(fd_nodal_force.iter())
        .for_each(|(nodal_force_i, fd_nodal_force_i)|
            assert!(
                (nodal_force_i/fd_nodal_force_i - 1.0).abs() < EPSILON
            )
        )
    )
}

// can probably test sum of shape functions at integration points is unity as well

#[test]
fn todo()
{
    todo!()
}