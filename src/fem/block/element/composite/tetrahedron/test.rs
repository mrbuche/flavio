use crate::fem::block::
{
    element::composite::test::test_composite_element,
    // test::test_finite_element_block
};
use super::*;

// const D: usize = 14;
// const E: usize = 24;

test_composite_element!(Tetrahedron);
// test_finite_element_block!(Tetrahedron);





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
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.5, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.0, 0.5, 0.0],
        [0.0, 0.0, 0.5],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5]
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
        [0.87900980, 0.73292419, 0.67130546],
        [0.89557721, 0.40375233, 0.94588738],
        [0.98513155, 0.51115036, 0.52239128],
        [0.42737582, 0.87092470, 0.20908760],
        [0.78803900, 0.83651295, 0.69155253],
        [0.23296807, 0.41294924, 0.26352866],
        [0.26965179, 0.63112838, 0.08563159],
        [0.19258050, 0.58921056, 0.87556116],
        [0.45519136, 0.90772602, 0.26139220],
        [0.47893746, 0.64353023, 0.04217477]
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
        let mut sum = [0.0_f64; 3];
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
            assert_eq_within_tols(sum_i, &0.0)
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
    let mut finite_difference = 0.0;
    (0..N).map(|node|
        (0..3).map(|i|{
            let mut nodal_coordinates = get_coordinatess();
            nodal_coordinates[node][i] += 0.5 * EPSILON;
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