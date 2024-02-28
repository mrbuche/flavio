use crate::fem::block::element::
{
    test::test_finite_element,
    linear::
    {
        test::test_linear_finite_element
    }
};
use super::*;

#[test]
fn test_finite_element()
{
    todo!()
}

#[test]
fn test_linear_finite_element()
{
    todo!()
}

fn get_reference_coordinates() -> ReferenceNodalCoordinates<N>
{
    ReferenceNodalCoordinates::new([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0]
    ])
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

fn get_element<'a>() -> Triangle<'a, NeoHookean<'a>>
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

fn get_reference_coordinates_crazy() -> ReferenceNodalCoordinates<N>
{
    ReferenceNodalCoordinates::new([
        [0.60728341, 0.65218218, 0.16792581],
        [0.67789658, 0.89610824, 0.30427157],
        [0.71685949, 0.71390134, 0.17156583]
    ])
}

fn get_element_crazy<'a>() -> Triangle<'a, NeoHookean<'a>>
{
    Triangle::<NeoHookean>::new(
        NEOHOOKEANPARAMETERS,
        get_reference_coordinates_crazy()
    )
}

fn get_coordinates_crazy() -> NodalCoordinates<N>
{
    NodalCoordinates::new([
        [0.56613414, 0.73732276, 0.81193896],
        [0.91018069, 0.82804788, 0.95319798],
        [0.58414216, 0.76766451, 0.26952952]
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

#[test]
fn temporary_3()
{
    get_element_crazy().calculate_nodal_forces(
        &get_coordinates_crazy()
    ).iter()
    .zip(get_finite_difference_of_helmholtz_free_energy_crazy(true).iter())
    .for_each(|(f_a, fd_a)|
        f_a.iter()
        .zip(fd_a.iter())
        .for_each(|(f_a_i, fd_a_i)|
            println!("{:?}", (f_a_i, fd_a_i))
        )
    )
}