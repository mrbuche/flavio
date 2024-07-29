use crate::
{
    EPSILON,
    constitutive::cohesive::
    {
        SmithFerrante,
        test::SMITHFERRANTEPARAMETERS
    }
};

use super::*;

fn get_reference_coordinates() -> ReferenceNodalCoordinates<N>
{
    ReferenceNodalCoordinates::new([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0]
    ])
}

fn get_coordinates() -> NodalCoordinates<N>
{
    NodalCoordinates::new([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 2.0],
        [0.0, 1.0, 3.0]
    ])
}

#[test]
fn zero()
{
    let element = Wedge::<SmithFerrante>::new(
        SMITHFERRANTEPARAMETERS,
        get_reference_coordinates()
    );
    element.calculate_nodal_forces(
        &get_reference_coordinates().into()
    ).iter()
    .for_each(|nodal_force|
        nodal_force.iter()
        .for_each(|nodal_force_i|
            assert_eq!(nodal_force_i, &0.0)
        )
    )
}

#[test]
fn finite_difference()
{
    let element = Wedge::<SmithFerrante>::new(
        SMITHFERRANTEPARAMETERS,
        get_reference_coordinates()
    );
    let mut finite_difference = 0.0;
    element.calculate_nodal_stiffnesses(
        &get_coordinates()
    ).iter()
    .enumerate()
    .for_each(|(a, nodal_stiffness_a)|
        nodal_stiffness_a.iter()
        .enumerate()
        .for_each(|(b, nodal_stiffness_ab)|
            nodal_stiffness_ab.iter()
            .enumerate()
            .for_each(|(i, nodal_stiffness_ab_i)|
                nodal_stiffness_ab_i.iter()
                .enumerate()
                .for_each(|(j, nodal_stiffness_ab_ij)|{
                    let mut nodal_coordinates = get_coordinates();
                    nodal_coordinates[b][j] += 0.5 * EPSILON;
                    finite_difference = element.calculate_nodal_forces(
                        &nodal_coordinates
                    )[a][i];
                    nodal_coordinates[b][j] -= EPSILON;
                    finite_difference -= element.calculate_nodal_forces(
                        &nodal_coordinates
                    )[a][i];
                    finite_difference /= EPSILON;
                    assert!(
                        (nodal_stiffness_ab_ij/finite_difference - 1.0).abs() < EPSILON
                    )
                })
            )
        )
    )
}
