use crate::constitutive::cohesive::SmithFerrante;

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
        &[1.0, 1.0],
        get_reference_coordinates()
    );
    element.calculate_nodal_forces(
        &get_reference_coordinates().convert()
    ).iter()
    .for_each(|nodal_force|
        nodal_force.iter()
        .for_each(|nodal_force_i|
            assert_eq!(nodal_force_i, &0.0)
        )
    )
}

#[test]
fn nonzero()
{
    let element = Wedge::<SmithFerrante>::new(
        &[1.0, 1.0, 1.0],
        get_reference_coordinates()
    );
    element.calculate_nodal_forces(
        &get_coordinates()
    ).iter()
    .for_each(|nodal_force|
        println!("{:?}", (nodal_force[0], nodal_force[1], nodal_force[2]))
    )
}