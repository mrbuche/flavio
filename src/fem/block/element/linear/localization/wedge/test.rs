use crate::fem::block::element::linear::localization::test::
{
    setup_for_test_linear_surface_element_with_constitutive_model,
    test_linear_localization_element
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

test_linear_localization_element!(Wedge);

// Do not forget the actual tests!

// should use random ref/curr coords for element/block tests
// and should apply surface/localiz. block tests of course

// maybe go back to CurrentNodalCoordinates and stuff and call Coordinates<I> NodalCoordinates<I>





use crate::
{
    constitutive::solid::hyperelastic::
    {
        NeoHookean,
        test::NEOHOOKEANPARAMETERS
    }
};

fn get_element<'a>() -> Wedge<NeoHookean<'a>>
{
    Wedge::<NeoHookean>::new(
        NEOHOOKEANPARAMETERS,
        get_reference_coordinates()
    )
}

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

fn get_coordinatesss() -> NodalCoordinates<N>
{
    get_deformation_gradient_surface()*get_reference_coordinates()
}

#[test]
fn temporary_1()
{
    get_element().calculate_deformation_gradient(
        &get_coordinatesss()
    ).iter().zip(get_deformation_gradient_surface().iter()).for_each(|(f_i, ff_i)|
        f_i.iter().zip(ff_i.iter()).for_each(|(f_ij, ff_ij)|
            assert!((f_ij/ff_ij - 1.0).abs() < 1e-8)
        )
    );
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
                get_coordinatesss()
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
fn temporary_2()
{
    get_element().calculate_nodal_forces(
        &get_coordinatesss()
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

fn get_coordinates_jump() -> NodalCoordinates<N>
{
    NodalCoordinates::new([
        [0.0, 0.0, -1.0],
        [1.0, 0.0, -1.0],
        [0.0, 1.0, -1.0],
        [0.0, 0.0,  1.0],
        [1.0, 0.0,  1.0],
        [0.0, 1.0,  1.0]
    ])
}

fn get_finite_difference_of_helmholtz_free_energy_jump(is_deformed: bool) -> NodalForces<N>
{
    let element = get_element();
    let mut finite_difference = 0.0;
    (0..N).map(|node|
        (0..3).map(|i|{
            let mut nodal_coordinates = 
            if is_deformed
            {
                get_coordinates_jump()
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
fn temporary_3()
{
    get_element().calculate_nodal_forces(
        &get_coordinates_jump()
    ).iter()
    .zip(get_finite_difference_of_helmholtz_free_energy_jump(true).iter())
    .for_each(|(f_a, fd_a)|
        f_a.iter()
        .zip(fd_a.iter())
        .for_each(|(f_a_i, fd_a_i)|
            assert!(
                (f_a_i/fd_a_i - 1.0).abs() < EPSILON || f_a_i == &0.0
            )
        )
    )
}

fn get_reference_coordinates_crazy() -> ReferenceNodalCoordinates<N>
{
    ReferenceNodalCoordinates::new([
        [0.67914347, 0.68755701, 0.62163983],
        [0.78461572, 0.42216947, 0.00126725],
        [0.32961799, 0.2031503 , 0.26778553],
        [0.67914347, 0.68755701, 0.62163983],
        [0.78461572, 0.42216947, 0.00126725],
        [0.32961799, 0.2031503 , 0.26778553]
    ])
}

fn get_element_crazy<'a>() -> Wedge<NeoHookean<'a>>
{
    Wedge::<NeoHookean>::new(
        NEOHOOKEANPARAMETERS,
        get_reference_coordinates_crazy()
    )
}

fn get_coordinates_crazy() -> NodalCoordinates<N>
{
    NodalCoordinates::new([
        [0.36328258, 0.64146886, 0.96629975],
        [0.77892449, 0.74881768, 0.98804891],
        [0.69872272, 0.75507056, 0.2489533 ],
        // [0.36328258, 0.64146886, 0.96629975],
        // [0.77892449, 0.74881768, 0.98804891],
        // [0.69872272, 0.75507056, 0.2489533 ],
        [0.82831486, 0.84934763, 0.44556648],
        [0.40287625, 0.36786279, 0.17704429],
        [0.37918941, 0.34669909, 0.16809232]
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
fn temporary_4()
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