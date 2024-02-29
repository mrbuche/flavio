use crate::fem::block::element::
{
    test::test_finite_element,
    linear::
    {
        test::test_linear_finite_element
    }
};
use super::*;

fn get_reference_coordinates() -> ReferenceNodalCoordinates<N>
{
    ReferenceNodalCoordinates::new([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0]
    ])
}

// consider putting a deformation gradient calculation test in linear/test.rs
// also a similar one for deformation gradient rate
// is for displacements according to affine deformation of course

test_finite_element!(Triangle);
test_linear_finite_element!(Triangle);

use crate::constitutive::solid::elastic::AlmansiHamel;

#[test]
fn size()
{
    assert_eq!(
        std::mem::size_of::<Tetrahedron::<AlmansiHamel>>(),
        std::mem::size_of::<AlmansiHamel>()
        + std::mem::size_of::<GradientVectors<N>>()
        + std::mem::size_of::<ReferenceNormal>()
    )
}






use crate::
{
    constitutive::solid::hyperelastic::
    {
        NeoHookean,
        test::NEOHOOKEANPARAMETERS
    },
    mechanics::
    {
        RotationCurrentConfiguration
    },
    EPSILON
};

fn get_deformation_gradient_surface() -> DeformationGradient
{
    // needs to be simplified and precise most likely to get tests to pass
    let rotation = RotationCurrentConfiguration::new([
        [-0.96152505,  0.15428784,  0.22729901],
        [-0.26402731, -0.29043967, -0.91974691],
        [-0.07588912, -0.94437284,  0.32000122]
    ]);
    let deformation_gradient = DeformationGradient::new([
        [0.61926467, 0.20573410, 0.0],
        [0.32283599, 0.98534316, 0.0],
        [0.00000000, 0.00000000, 1.0]
    ]);
    rotation * deformation_gradient
}

fn get_coordinates() -> NodalCoordinates<N>
{
    get_deformation_gradient_surface()*get_reference_coordinates()
}

fn get_element<'a>() -> Triangle<NeoHookean<'a>>
{
    Triangle::<NeoHookean>::new(
        NEOHOOKEANPARAMETERS,
        get_reference_coordinates()
    )
}

fn get_finite_difference_of_helmholtz_free_energy(is_deformed: bool) -> NodalForces<N>
{
    let element = get_element();
    let mut finite_difference = 0.0;
    (0..N).map(|node|
        (0..3).map(|i|{
            let mut nodal_coordinates = 
            if is_deformed
            {
                get_coordinates()
            }
            else
            {
                get_reference_coordinates().convert()
            };
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
fn temporary_1()
{
    get_element().calculate_deformation_gradient(
        &get_coordinates()
    ).iter().zip(get_deformation_gradient_surface().iter()).for_each(|(f_i, ff_i)|
        f_i.iter().zip(ff_i.iter()).for_each(|(f_ij, ff_ij)|
            assert!((f_ij/ff_ij - 1.0).abs() < 1e-8)
        )
    );
}

#[test]
fn temporary_2()
{
    get_element().calculate_nodal_forces(
        &get_coordinates()
    ).iter()
    .zip(get_finite_difference_of_helmholtz_free_energy(true).iter())
    .for_each(|(f_a, fd_a)|
        f_a.iter()
        .zip(fd_a.iter())
        .for_each(|(f_a_i, fd_a_i)|
            assert!(
                (f_a_i/fd_a_i - 1.0).abs() < EPSILON
            )
        )
    )
}

fn get_reference_coordinates_crazy() -> ReferenceNodalCoordinates<N>
{
    ReferenceNodalCoordinates::new([
        [0.60728341, 0.65218218, -1.16792581],
        [0.67789658, -2.89610824, 0.30427157],
        [0.71685949, 0.71390134, 0.17156583]
    ])
}

fn get_element_crazy<'a>() -> Triangle<NeoHookean<'a>>
{
    Triangle::<NeoHookean>::new(
        NEOHOOKEANPARAMETERS,
        get_reference_coordinates_crazy()
    )
}

fn get_coordinates_crazy() -> NodalCoordinates<N>
{
    NodalCoordinates::new([
        [1.56613414, 0.73732276, 0.81193896],
        [0.91018069, 0.82804788, 0.95319798],
        [-2.58414216, 0.76766451, 0.26952952]
    ])
}

fn get_finite_difference_of_helmholtz_free_energy_crazy(is_deformed: bool) -> NodalForces<N>
{
    let element = get_element_crazy();
    let mut finite_difference = 0.0;
    (0..N).map(|node|
        (0..3).map(|i|{
            let mut nodal_coordinates = 
            if is_deformed
            {
                get_coordinates_crazy()
            }
            else
            {
                get_reference_coordinates_crazy().convert()
            };
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
fn temporary_3()
{
    get_element_crazy().calculate_nodal_forces(
        &get_coordinates_crazy()
    ).iter()
    .zip(get_finite_difference_of_helmholtz_free_energy_crazy(true).iter())
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
fn get_finite_difference_of_nodal_forces_crazy(is_deformed: bool) -> NodalStiffnesses<N>
{
    let element = get_element_crazy();
    let mut finite_difference = 0.0;
    (0..N).map(|node_a|
        (0..N).map(|node_b|
            (0..3).map(|i|
                (0..3).map(|j|{
                    let mut nodal_coordinates = 
                    if is_deformed
                    {
                        get_coordinates_crazy()
                    }
                    else
                    {
                        get_reference_coordinates_crazy().convert()
                    };
                    nodal_coordinates[node_a][i] += 0.5 * EPSILON;
                    finite_difference = element.calculate_nodal_forces(
                        &nodal_coordinates
                    )[node_b][j];
                    nodal_coordinates[node_a][i] -= EPSILON;
                    finite_difference -= element.calculate_nodal_forces(
                        &nodal_coordinates
                    )[node_b][j];
                    finite_difference/EPSILON
                }).collect()
            ).collect()
        ).collect()
    ).collect()
}

#[test]
fn temporary_4()
{
    get_element_crazy().calculate_nodal_stiffnesses(
        &get_coordinates_crazy()
    ).iter()
    .zip(get_finite_difference_of_nodal_forces_crazy(true).iter())
    .for_each(|(nodal_stiffness_a, fd_nodal_stiffness_a)|
        nodal_stiffness_a.iter()
        .zip(fd_nodal_stiffness_a.iter())
        .for_each(|(nodal_stiffness_ab, fd_nodal_stiffness_ab)|
            nodal_stiffness_ab.iter()
            .zip(fd_nodal_stiffness_ab.iter())
            .for_each(|(nodal_stiffness_ab_i, fd_nodal_stiffness_ab_i)|
                nodal_stiffness_ab_i.iter()
                .zip(fd_nodal_stiffness_ab_i.iter())
                .for_each(|(nodal_stiffness_ab_ij, fd_nodal_stiffness_ab_ij)|
                    // assert!(
                    //     (nodal_stiffness_ab_ij/fd_nodal_stiffness_ab_ij - 1.0).abs() < EPSILON
                    // )
                    println!("{:?}", (nodal_stiffness_ab_ij, fd_nodal_stiffness_ab_ij))
                )
            )
        )
    )
}