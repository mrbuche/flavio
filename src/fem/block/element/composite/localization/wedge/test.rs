use super::*;

fn get_reference_coordinates() -> ReferenceNodalCoordinates<N>
{
    ReferenceNodalCoordinates::new([
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.0, 0.5, 0.0],
        [0.5, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.0, 0.5, 0.0],
        [0.5, 0.0, 0.0]
    ])
}

#[test]
fn test_composite_localization_element()
{
    todo!()
}

#[test]
fn test_finite_element_block()
{
    todo!()
}

fn get_element<'a>() -> Wedge<crate::constitutive::solid::hyperelastic::NeoHookean<'a>>
{
    Wedge::new(
        crate::constitutive::solid::hyperelastic::test::NEOHOOKEANPARAMETERS,
        get_reference_coordinates()
    )
}

#[test]
fn temporary_1()
{
    get_element().calculate_deformation_gradients(
        &get_reference_coordinates().convert()
    ).iter()
    .for_each(|deformation_gradient|
        deformation_gradient.iter()
        .enumerate()
        .for_each(|(i, deformation_gradient_i)|
            deformation_gradient_i.iter()
            .enumerate()
            .for_each(|(j, deformation_gradient_ij)|
                crate::test::assert_eq_within_tols(
                    deformation_gradient_ij, &(((i == j) as u8) as Scalar)
                )
            )
        )
    )
}

fn get_deformation_gradient_planar() -> crate::mechanics::DeformationGradient
{
    crate::mechanics::DeformationGradient::new([
        [0.9, 0.5, 0.0],
        [0.7, 0.9, 0.0],
        [0.0, 0.0, 1.0]
    ])
}

fn get_coordinates_from_deformation_gradient_planar() -> NodalCoordinates<N>
{
    NodalCoordinates::new([
        [0.90, 0.70, 0.0],
        [0.50, 0.90, 0.0],
        [0.00, 0.00, 0.0],
        [0.90, 0.70, 0.0],
        [0.50, 0.90, 0.0],
        [0.00, 0.00, 0.0],
        [0.70, 0.80, 0.0],
        [0.25, 0.45, 0.0],
        [0.45, 0.35, 0.0],
        [0.70, 0.80, 0.0],
        [0.25, 0.45, 0.0],
        [0.45, 0.35, 0.0]
    ])
}

#[test]
fn temporary_2()
{
    get_element().calculate_deformation_gradients(
        &get_coordinates_from_deformation_gradient_planar()
    ).iter()
    .for_each(|calculated_deformation_gradient|
        calculated_deformation_gradient.iter()
        .zip(get_deformation_gradient_planar().iter())
        .for_each(|(calculated_deformation_gradient_i, deformation_gradient_i)|
            calculated_deformation_gradient_i.iter()
            .zip(deformation_gradient_i.iter())
            .for_each(|(calculated_deformation_gradient_ij, deformation_gradient_ij)|
                crate::test::assert_eq_within_tols(
                    calculated_deformation_gradient_ij, deformation_gradient_ij
                )
            )
        )
    )
}

fn get_deformation_gradient_both() -> crate::mechanics::DeformationGradient
{
    crate::mechanics::DeformationGradient::new([
        [0.9, 0.5, 0.2],
        [0.7, 0.9, 0.4],
        [0.0, 0.0, 1.6]
    ])
}

fn get_coordinates_from_deformation_gradient_both() -> NodalCoordinates<N>
{
    NodalCoordinates::new([
        [ 0.80,  0.50, -0.3],
        [ 0.40,  0.70, -0.3],
        [-0.10, -0.20, -0.3],
        [ 1.00,  0.90,  0.3],
        [ 0.60,  1.10,  0.3],
        [ 0.10,  0.20,  0.3],
        [ 0.60,  0.60, -0.3],
        [ 0.15,  0.25, -0.3],
        [ 0.35,  0.15, -0.3],
        [ 0.80,  1.00,  0.3],
        [ 0.35,  0.65,  0.3],
        [ 0.55,  0.55,  0.3]
    ])
}

#[test]
fn temporary_3()
{
    get_element().calculate_deformation_gradients(
        &get_coordinates_from_deformation_gradient_both()
    ).iter()
    .for_each(|calculated_deformation_gradient|
        calculated_deformation_gradient.iter()
        .zip(get_deformation_gradient_both().iter())
        .for_each(|(calculated_deformation_gradient_i, deformation_gradient_i)|
            calculated_deformation_gradient_i.iter()
            .zip(deformation_gradient_i.iter())
            .for_each(|(calculated_deformation_gradient_ij, deformation_gradient_ij)|
                crate::test::assert_eq_within_tols(
                    calculated_deformation_gradient_ij, deformation_gradient_ij
                )
            )
        )
    )
}

fn get_finite_difference_of_helmholtz_free_energy_planar() -> NodalForces<N>
{
    let element = get_element();
    let mut finite_difference = 0.0;
    (0..N).map(|node|
        (0..3).map(|i|{
            let mut nodal_coordinates = get_coordinates_from_deformation_gradient_planar();
            nodal_coordinates[node][i] += 0.5 * crate::EPSILON;
            finite_difference = element.calculate_helmholtz_free_energy(
                &nodal_coordinates
            );
            nodal_coordinates[node][i] -= crate::EPSILON;
            finite_difference -= element.calculate_helmholtz_free_energy(
                &nodal_coordinates
            );
            finite_difference/crate::EPSILON
        }).collect()
    ).collect()
}

#[test]
fn temporary_4()
{
    get_element().calculate_nodal_forces(
        &get_coordinates_from_deformation_gradient_planar()
    ).iter()
    .zip(get_finite_difference_of_helmholtz_free_energy_planar().iter())
    .for_each(|(nodal_force, fd_nodal_force)|
        nodal_force.iter()
        .zip(fd_nodal_force.iter())
        .for_each(|(nodal_force_i, fd_nodal_force_i)|
            assert!(
                (nodal_force_i/fd_nodal_force_i - 1.0).abs() < crate::EPSILON
            )
        )
    )
}

fn get_reference_coordinates_distorted() -> ReferenceNodalCoordinates<N>
{
    ReferenceNodalCoordinates::new([
        [ 0.99926759, -0.00492477,  0.00473545],
        [ 0.00395293,  1.01517263,  0.00579895],
        [ 0.0084215 ,  0.01449531, -0.01263632],
        [ 0.98528723, -0.00700808,  0.01811121],
        [ 0.01031553,  0.98143822, -0.0055957 ],
        [ 0.00599955, -0.01477466, -0.00564895],
        [ 0.51456154,  0.5120805 , -0.01684216],
        [ 0.00283855,  0.50947624, -0.00354348],
        [ 0.50324758, -0.00521494, -0.01022159],
        [ 0.50758104,  0.50456292, -0.01654888],
        [ 0.00106299,  0.51187361,  0.01741694],
        [ 0.49211042, -0.01282023,  0.01630969]
    ])
}

fn get_element_distorted<'a>() -> Wedge<crate::constitutive::solid::hyperelastic::NeoHookean<'a>>
{
    Wedge::new(
        crate::constitutive::solid::hyperelastic::test::NEOHOOKEANPARAMETERS,
        get_reference_coordinates_distorted()
    )
}

fn get_coordinates_distorted() -> NodalCoordinates<N>
{
    NodalCoordinates::new([
        [ 0.78126924,  0.48248601, -0.31196753],
        [ 0.39495532,  0.68137112, -0.29016521],
        [-0.11102563, -0.19592438, -0.3130127 ],
        [ 1.01338446,  0.88479018,  0.29831012],
        [ 0.61690543,  1.11171078,  0.30851072],
        [ 0.1077119 ,  0.20115359,  0.30302711],
        [ 0.60260084,  0.58991677, -0.28643558],
        [ 0.15742686,  0.2610388 , -0.30296679],
        [ 0.35750349,  0.14907168, -0.31712706],
        [ 0.78279103,  1.00726579,  0.31962874],
        [ 0.34129542,  0.66540902,  0.3152397 ],
        [ 0.53464796,  0.56976664,  0.30936124]
    ])
}

fn get_finite_difference_of_helmholtz_free_energy_distorted() -> NodalForces<N>
{
    let element = get_element_distorted();
    let mut finite_difference = 0.0;
    (0..N).map(|node|
        (0..3).map(|i|{
            let mut nodal_coordinates = get_coordinates_distorted();
            nodal_coordinates[node][i] += 0.5 * crate::EPSILON;
            finite_difference = element.calculate_helmholtz_free_energy(
                &nodal_coordinates
            );
            nodal_coordinates[node][i] -= crate::EPSILON;
            finite_difference -= element.calculate_helmholtz_free_energy(
                &nodal_coordinates
            );
            finite_difference / crate::EPSILON
        }).collect()
    ).collect()
}

#[test]
fn temporary_5()
{
    get_element_distorted().calculate_nodal_forces(
        &get_coordinates_distorted()
    ).iter()
    .zip(get_finite_difference_of_helmholtz_free_energy_distorted().iter())
    .for_each(|(nodal_force, fd_nodal_force)|
        nodal_force.iter()
        .zip(fd_nodal_force.iter())
        .for_each(|(nodal_force_i, fd_nodal_force_i)|
            assert!(
                (nodal_force_i/fd_nodal_force_i - 1.0).abs() < crate::EPSILON
            )
        )
    )
}

fn get_finite_difference_of_nodal_forces_distorted() -> NodalStiffnesses<N>
{
    let element = get_element_distorted();
    let mut finite_difference = 0.0;
    (0..N).map(|a|
        (0..N).map(|b|
            (0..3).map(|i|
                (0..3).map(|j|{
                    let mut nodal_coordinates = get_coordinates_distorted();
                    nodal_coordinates[b][j] += 0.5 * crate::EPSILON;
                    finite_difference = element.calculate_nodal_forces(
                        &nodal_coordinates
                    )[a][i];
                    nodal_coordinates[b][j] -= crate::EPSILON;
                    finite_difference -= element.calculate_nodal_forces(
                        &nodal_coordinates
                    )[a][i];
                    finite_difference / crate::EPSILON
                }).collect()
            ).collect()
        ).collect()
    ).collect()
}

#[test]
fn temporary_6()
{
    get_element_distorted().calculate_nodal_stiffnesses(
        &get_coordinates_distorted()
    ).iter()
    .zip(get_finite_difference_of_nodal_forces_distorted().iter())
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
                    assert!(
                        (nodal_stiffness_ab_ij/fd_nodal_stiffness_ab_ij - 1.0).abs() < crate::EPSILON ||
                        (nodal_stiffness_ab_ij - fd_nodal_stiffness_ab_ij).abs() < crate::EPSILON / 10.0
                    )
                )
            )
        )
    )
}