macro_rules! test_finite_element
{
    ($element: ident) =>
    {
        mod finite_element
        {
            use crate::
            {
                EPSILON,
                fem::block::element::test::
                {
                    test_finite_element_with_elastic_constitutive_model,
                    test_finite_element_with_hyperelastic_constitutive_model,
                    test_finite_element_with_hyperviscoelastic_constitutive_model
                },
                math::
                {
                    Convert,
                    TensorRank2
                },
                mechanics::test::
                {
                    get_deformation_gradient,
                    get_rotation_current_configuration,
                    get_rotation_reference_configuration,
                    get_translation_current_configuration,
                    get_translation_reference_configuration
                },
                test::assert_eq_within_tols
            };
            use super::*;
            mod elastic
            {
                use super::*;
                mod almansi_hamel
                {
                    use crate::
                    {
                        constitutive::solid::elastic::
                        {
                            AlmansiHamel,
                            test::ALMANSIHAMELPARAMETERS
                        }
                    };
                    use super::*;
                    test_finite_element_with_elastic_constitutive_model!($element, AlmansiHamel, ALMANSIHAMELPARAMETERS);
                }
            }
            mod hyperelastic
            {
                use crate::
                {
                    constitutive::solid::hyperelastic::
                    {
                        ArrudaBoyce,
                        Fung,
                        Gent,
                        MooneyRivlin,
                        NeoHookean,
                        SaintVenantKirchoff,
                        Yeoh,
                        test::
                        {
                            ARRUDABOYCEPARAMETERS,
                            FUNGPARAMETERS,
                            GENTPARAMETERS,
                            MOONEYRIVLINPARAMETERS,
                            NEOHOOKEANPARAMETERS,
                            SAINTVENANTKIRCHOFFPARAMETERS,
                            YEOHPARAMETERS
                        }
                    }
                };
                use super::*;
                mod arruda_boyce
                {
                    use super::*;
                    test_finite_element_with_hyperelastic_constitutive_model!($element, ArrudaBoyce, ARRUDABOYCEPARAMETERS);
                }
                mod fung
                {
                    use super::*;
                    test_finite_element_with_hyperelastic_constitutive_model!($element, Fung, FUNGPARAMETERS);
                }
                mod gent
                {
                    use super::*;
                    test_finite_element_with_hyperelastic_constitutive_model!($element, Gent, GENTPARAMETERS);
                }
                mod mooney_rivlin
                {
                    use super::*;
                    test_finite_element_with_hyperelastic_constitutive_model!($element, MooneyRivlin, MOONEYRIVLINPARAMETERS);
                }
                mod neo_hookean
                {
                    use super::*;
                    test_finite_element_with_hyperelastic_constitutive_model!($element, NeoHookean, NEOHOOKEANPARAMETERS);
                }
                mod saint_venant_kirchoff
                {
                    use super::*;
                    test_finite_element_with_hyperelastic_constitutive_model!($element, SaintVenantKirchoff, SAINTVENANTKIRCHOFFPARAMETERS);
                }
                mod yeoh
                {
                    use super::*;
                    test_finite_element_with_hyperelastic_constitutive_model!($element, Yeoh, YEOHPARAMETERS);
                }
            }
            mod hyperviscoelastic
            {
                use crate::
                {
                    constitutive::solid::hyperviscoelastic::
                    {
                        SaintVenantKirchoff,
                        test::SAINTVENANTKIRCHOFFPARAMETERS
                    }
                };
                use super::*;
                mod saint_venant_kirchoff
                {
                    use super::*;
                    test_finite_element_with_hyperviscoelastic_constitutive_model!($element, SaintVenantKirchoff, SAINTVENANTKIRCHOFFPARAMETERS);
                }
            }
        }
    }
}
pub(crate) use test_finite_element;

macro_rules! test_finite_element_with_solid_constitutive_model
{
    ($element: ident, $constitutive_model: ident, $constitutive_model_parameters: ident) =>
    {
        fn get_coordinates() -> NodalCoordinates<N>
        {
            get_reference_coordinates().iter()
            .map(|reference_coordinate|
                get_deformation_gradient() * reference_coordinate
            ).collect()
        }
        fn get_coordinates_transformed() -> NodalCoordinates<N>
        {
            get_coordinates().iter()
            .map(|current_coordinate|
                get_rotation_current_configuration() * current_coordinate
                + get_translation_current_configuration()
            ).collect()
        }
        fn get_element<'a>() -> $element<$constitutive_model<'a>>
        {
            $element::new(
                $constitutive_model_parameters,
                get_reference_coordinates()
            )
        }
        fn get_element_transformed<'a>() -> $element<$constitutive_model<'a>>
        {
            $element::<$constitutive_model>::new
            (
                $constitutive_model_parameters,
                get_reference_coordinates_transformed()
            )
        }
        fn get_reference_coordinates_transformed() -> ReferenceNodalCoordinates<N>
        {
            get_reference_coordinates().iter()
            .map(|reference_coordinate|
                get_rotation_reference_configuration() * reference_coordinate
                + get_translation_reference_configuration()
            ).collect()
        }
        fn get_velocities() -> NodalVelocities<N>
        {
            todo!()
        }
        #[test]
        fn integration_weights_sum_to_one()
        {
            assert_eq!(get_element().get_integration_weights().iter().sum::<Scalar>(), 1.0)
        }
        mod nodal_forces
        {
            use super::*;
            mod deformed
            {
                use super::*;
                #[test]
                fn finite_difference()
                {
                    get_nodal_stiffnesses(true, false).iter()
                    .zip(get_finite_difference_of_nodal_forces(true).iter())
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
                                        (nodal_stiffness_ab_ij/fd_nodal_stiffness_ab_ij - 1.0).abs() < EPSILON
                                    )
                                )
                            )
                        )
                    )
                }
                #[test]
                fn objectivity()
                {
                    get_nodal_forces(true, false).iter()
                    .zip(get_nodal_forces(true, true).iter())
                    .for_each(|(nodal_force, res_nodal_force)|
                        nodal_force.iter()
                        .zip(res_nodal_force.iter())
                        .for_each(|(nodal_force_i, res_nodal_force_i)|
                            assert_eq_within_tols(
                                nodal_force_i, res_nodal_force_i
                            )
                        )
                    )
                }
            }
            mod undeformed
            {
                use super::*;
                #[test]
                fn finite_difference()
                {
                    get_nodal_stiffnesses(false, false).iter()
                    .zip(get_finite_difference_of_nodal_forces(false).iter())
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
                                        (nodal_stiffness_ab_ij/fd_nodal_stiffness_ab_ij - 1.0).abs() < EPSILON ||
                                        (nodal_stiffness_ab_ij.abs() < EPSILON && fd_nodal_stiffness_ab_ij.abs() < EPSILON)
                                    )
                                )
                            )
                        )
                    )
                }
                #[test]
                fn objectivity()
                {
                    get_nodal_forces(false, false).iter()
                    .for_each(|nodal_force|
                        nodal_force.iter()
                        .for_each(|nodal_force_i|
                            assert_eq_within_tols(
                                nodal_force_i, &0.0
                            )
                        )
                    )
                }
                #[test]
                fn zero()
                {
                    get_nodal_forces(false, false).iter()
                    .for_each(|nodal_force|
                        nodal_force.iter()
                        .for_each(|nodal_force_i|
                            assert_eq_within_tols(
                                nodal_force_i, &0.0
                            )
                        )
                    )
                }
            }
        }
        mod nodal_stiffnesses
        {
            use super::*;
            mod deformed
            {
                use super::*;
                #[test]
                fn objectivity()
                {
                    get_nodal_stiffnesses(true, false).iter()
                    .zip(get_nodal_stiffnesses(true, true).iter())
                    .for_each(|(nodal_stiffness_a, res_nodal_stiffness_a)|
                        nodal_stiffness_a.iter()
                        .zip(res_nodal_stiffness_a.iter())
                        .for_each(|(nodal_stiffness_ab, res_nodal_stiffness_ab)|
                            nodal_stiffness_ab.iter()
                            .zip(res_nodal_stiffness_ab.iter())
                            .for_each(|(nodal_stiffness_ab_i, res_nodal_stiffness_ab_i)|
                                nodal_stiffness_ab_i.iter()
                                .zip(res_nodal_stiffness_ab_i.iter())
                                .for_each(|(nodal_stiffness_ab_ij, res_nodal_stiffness_ab_ij)|
                                    assert_eq_within_tols(
                                        nodal_stiffness_ab_ij, res_nodal_stiffness_ab_ij
                                    )
                                )
                            )
                        )
                    )
                }
                #[test]
                fn symmetry()
                {
                    let nodal_stiffness = get_nodal_stiffnesses(true, false);
                    nodal_stiffness.iter()
                    .enumerate()
                    .for_each(|(a, nodal_stiffness_a)|
                        nodal_stiffness_a.iter()
                        .enumerate()
                        .for_each(|(b, nodal_stiffness_ab)|
                            nodal_stiffness_ab.iter()
                            .enumerate()
                            .zip(nodal_stiffness_ab.transpose().iter())
                            .for_each(|((i, nodal_stiffness_ab_i), nodal_stiffness_ab_j)|
                                nodal_stiffness_ab_i.iter()
                                .enumerate()
                                .zip(nodal_stiffness_ab_j.iter())
                                .for_each(|((j, nodal_stiffness_ab_ij), nodal_stiffness_ab_ji)|
                                    if a == b
                                    {
                                        assert_eq_within_tols(
                                            nodal_stiffness_ab_ij,
                                            &nodal_stiffness_ab_ji
                                        )
                                    }
                                    else if i == j
                                    {
                                        assert_eq_within_tols(
                                            nodal_stiffness_ab_ij, &nodal_stiffness[b][a][i][j]
                                        )
                                    }
                                )
                            )
                        )
                    )
                }
            }
            mod undeformed
            {
                use super::*;
                #[test]
                fn objectivity()
                {
                    get_nodal_stiffnesses(false, false).iter()
                    .zip(get_nodal_stiffnesses(false, true).iter())
                    .for_each(|(nodal_stiffness_a, res_nodal_stiffness_a)|
                        nodal_stiffness_a.iter()
                        .zip(res_nodal_stiffness_a.iter())
                        .for_each(|(nodal_stiffness_ab, res_nodal_stiffness_ab)|
                            nodal_stiffness_ab.iter()
                            .zip(res_nodal_stiffness_ab.iter())
                            .for_each(|(nodal_stiffness_ab_i, res_nodal_stiffness_ab_i)|
                                nodal_stiffness_ab_i.iter()
                                .zip(res_nodal_stiffness_ab_i.iter())
                                .for_each(|(nodal_stiffness_ab_ij, res_nodal_stiffness_ab_ij)|
                                    assert_eq_within_tols(
                                        nodal_stiffness_ab_ij, res_nodal_stiffness_ab_ij
                                    )
                                )
                            )
                        )
                    )
                }
                #[test]
                fn symmetry()
                {
                    let nodal_stiffness = get_nodal_stiffnesses(false, false);
                    nodal_stiffness.iter()
                    .enumerate()
                    .for_each(|(a, nodal_stiffness_a)|
                        nodal_stiffness_a.iter()
                        .enumerate()
                        .for_each(|(b, nodal_stiffness_ab)|
                            nodal_stiffness_ab.iter()
                            .enumerate()
                            .zip(nodal_stiffness_ab.transpose().iter())
                            .for_each(|((i, nodal_stiffness_ab_i), nodal_stiffness_ab_j)|
                                nodal_stiffness_ab_i.iter()
                                .enumerate()
                                .zip(nodal_stiffness_ab_j.iter())
                                .for_each(|((j, nodal_stiffness_ab_ij), nodal_stiffness_ab_ji)|
                                    if a == b
                                    {
                                        assert_eq_within_tols(
                                            nodal_stiffness_ab_ij,
                                            &nodal_stiffness_ab_ji
                                        )
                                    }
                                    else if i == j
                                    {
                                        assert_eq_within_tols(
                                            nodal_stiffness_ab_ij, &nodal_stiffness[b][a][i][j]
                                        )
                                    }
                                )
                            )
                        )
                    )
                }
            }
        }
    }
}
pub(crate) use test_finite_element_with_solid_constitutive_model;

macro_rules! test_finite_element_with_elastic_constitutive_model
{
    ($element: ident, $constitutive_model: ident, $constitutive_model_parameters: ident) =>
    {
        fn get_nodal_forces(is_deformed: bool, is_rotated: bool) -> NodalForces<N>
        {
            if is_rotated
            {
                if is_deformed
                {
                    get_rotation_current_configuration().transpose() *
                    get_element_transformed().calculate_nodal_forces(
                        &get_coordinates_transformed()
                    )
                }
                else
                {
                    panic!()
                }
            }
            else
            {
                if is_deformed
                {
                    get_element().calculate_nodal_forces(
                        &get_coordinates()
                    )
                }
                else
                {
                    get_element().calculate_nodal_forces(
                        &get_reference_coordinates().convert()
                    )
                }
            }
        }
        fn get_nodal_stiffnesses(is_deformed: bool, is_rotated: bool) -> NodalStiffnesses<N>
        {
            if is_rotated
            {
                if is_deformed
                {
                    get_rotation_current_configuration().transpose() *
                    get_element_transformed().calculate_nodal_stiffnesses(
                        &get_coordinates_transformed()
                    ) * get_rotation_current_configuration()
                }
                else
                {
                    let converted: TensorRank2<3, 1, 1> = get_rotation_reference_configuration().convert();
                    converted.transpose() *
                    get_element_transformed().calculate_nodal_stiffnesses(
                        &get_reference_coordinates_transformed().convert()
                    ) * converted
                }
            }
            else
            {
                if is_deformed
                {
                    get_element().calculate_nodal_stiffnesses(
                        &get_coordinates()
                    )
                }
                else
                {
                    get_element().calculate_nodal_stiffnesses(
                        &get_reference_coordinates().convert()
                    )
                }
            }
        }
        fn get_finite_difference_of_nodal_forces(is_deformed: bool) -> NodalStiffnesses<N>
        {
            let element = get_element();
            let mut finite_difference = 0.0;
            (0..N).map(|node_a|
                (0..N).map(|node_b|
                    (0..3).map(|i|
                        (0..3).map(|j|{
                            let mut nodal_coordinates = 
                            if is_deformed
                            {
                                get_coordinates()
                            }
                            else
                            {
                                get_reference_coordinates().convert()
                            };
                            nodal_coordinates[node_a][i] += 0.5 * EPSILON;
                            finite_difference = element.calculate_nodal_forces(&nodal_coordinates)[node_b][j];
                            nodal_coordinates[node_a][i] -= EPSILON;
                            finite_difference -= element.calculate_nodal_forces(&nodal_coordinates)[node_b][j];
                            finite_difference/EPSILON
                        }).collect()
                    ).collect()
                ).collect()
            ).collect()
        }
        crate::fem::block::element::test::test_finite_element_with_solid_constitutive_model!(
            $element, $constitutive_model, $constitutive_model_parameters
        );
    }
}
pub(crate) use test_finite_element_with_elastic_constitutive_model;

macro_rules! test_finite_element_with_hyperelastic_constitutive_model
{
    ($element: ident, $constitutive_model: ident, $constitutive_model_parameters: ident) =>
    {
        crate::fem::block::element::test::test_finite_element_with_elastic_constitutive_model!(
            $element, $constitutive_model, $constitutive_model_parameters
        );
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
                    finite_difference = element.calculate_helmholtz_free_energy(&nodal_coordinates);
                    nodal_coordinates[node][i] -= EPSILON;
                    finite_difference -= element.calculate_helmholtz_free_energy(&nodal_coordinates);
                    finite_difference/EPSILON
                }).collect()
            ).collect()
        }
        mod helmholtz_free_energy
        {
            use super::*;
            mod deformed
            {
                use super::*;
                #[test]
                fn finite_difference()
                {
                    get_element().calculate_nodal_forces(
                        &get_coordinates()
                    ).iter()
                    .zip(get_finite_difference_of_helmholtz_free_energy(
                        true
                    ).iter())
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
                #[test]
                fn minimized()
                {
                    let element = get_element();
                    let nodal_forces = element.calculate_nodal_forces(
                        &get_coordinates()
                    );
                    let minimum = element.calculate_helmholtz_free_energy(
                        &get_coordinates()
                    ) - nodal_forces.dot(
                        &get_coordinates()
                    );
                    let mut perturbed_coordinates = get_coordinates();
                    (0..N).for_each(|node|
                        (0..3).for_each(|i|{
                            perturbed_coordinates = get_coordinates();
                            perturbed_coordinates[node][i] += 0.5 * EPSILON;
                            assert!(
                                element.calculate_helmholtz_free_energy(
                                    &perturbed_coordinates
                                ) - nodal_forces.dot(
                                    &perturbed_coordinates
                                ) > minimum
                            );
                            perturbed_coordinates[node][i] -= EPSILON;
                            assert!(
                                element.calculate_helmholtz_free_energy(
                                    &perturbed_coordinates
                                ) - nodal_forces.dot(
                                    &perturbed_coordinates
                                ) > minimum
                            );
                        })
                    )
                }
                #[test]
                fn objectivity()
                {
                    assert_eq_within_tols(
                        &get_element()
                        .calculate_helmholtz_free_energy(
                            &get_coordinates()
                        ),
                        &get_element_transformed()
                        .calculate_helmholtz_free_energy(
                            &get_coordinates_transformed()
                        )
                    )
                }
                #[test]
                fn positive()
                {
                    assert!(
                        get_element()
                        .calculate_helmholtz_free_energy(
                            &get_coordinates()
                        ) > 0.0
                    )
                }
            }
            mod undeformed
            {
                use super::*;
                #[test]
                fn finite_difference()
                {
                    get_finite_difference_of_helmholtz_free_energy(
                        false
                    ).iter()
                    .for_each(|fd_nodal_force|
                        fd_nodal_force.iter()
                        .for_each(|fd_nodal_force_i|
                            assert!(
                                fd_nodal_force_i.abs() < EPSILON
                            )
                        )
                    )
                }
                #[test]
                fn minimized()
                {
                    let element = get_element();
                    let minimum = element.calculate_helmholtz_free_energy(
                        &get_reference_coordinates().convert()
                    );
                    let mut perturbed_coordinates = get_reference_coordinates();
                    (0..N).for_each(|node|
                        (0..3).for_each(|i|{
                            perturbed_coordinates = get_reference_coordinates();
                            perturbed_coordinates[node][i] += 0.5 * EPSILON;
                            assert!(
                                element.calculate_helmholtz_free_energy(
                                    &perturbed_coordinates.convert()
                                ) > minimum
                            );
                            perturbed_coordinates[node][i] -= EPSILON;
                            assert!(
                                element.calculate_helmholtz_free_energy(
                                    &perturbed_coordinates.convert()
                                ) > minimum
                            );
                        })
                    )
                }
                #[test]
                fn objectivity()
                {
                    assert_eq_within_tols(
                        &get_element()
                        .calculate_helmholtz_free_energy(
                            &get_reference_coordinates().convert()
                        ),
                        &get_element_transformed()
                        .calculate_helmholtz_free_energy(
                            &get_reference_coordinates_transformed().convert()
                        )
                    )
                }
                #[test]
                fn zero()
                {
                    assert_eq!(
                        get_element()
                        .calculate_helmholtz_free_energy(
                            &get_reference_coordinates().convert()
                        ), 0.0
                    )
                }
            }
        }
    }
}
pub(crate) use test_finite_element_with_hyperelastic_constitutive_model;

macro_rules! test_finite_element_with_viscoelastic_constitutive_model
{
    ($element: ident, $constitutive_model: ident, $constitutive_model_parameters: ident) =>
    {
        fn get_nodal_forces(is_deformed: bool, is_rotated: bool) -> NodalForces<N>
        {
            if is_deformed
            {
                get_element().calculate_nodal_forces(
                    &get_coordinates(), &get_velocities()
                )
            }
            else
            {
                get_element().calculate_nodal_forces(
                    &get_reference_coordinates().convert(), &NodalVelocities::zero()
                )
            }
        }
        fn get_nodal_stiffnesses(is_deformed: bool, is_rotated: bool) -> NodalStiffnesses<N>
        {
            if is_deformed
            {
                get_element().calculate_nodal_stiffnesses(
                    &get_coordinates(), &get_velocities()
                )
            }
            else
            {
                get_element().calculate_nodal_stiffnesses(
                    &get_reference_coordinates().convert(), &NodalVelocities::zero()
                )
            }
        }
        fn get_finite_difference_of_nodal_forces(is_deformed: bool) -> NodalStiffnesses<N>
        {
            todo!()
        }
        crate::fem::block::element::test::test_finite_element_with_solid_constitutive_model!(
            $element, $constitutive_model, $constitutive_model_parameters
        );
    }
}
pub(crate) use test_finite_element_with_viscoelastic_constitutive_model;

macro_rules! test_finite_element_with_hyperviscoelastic_constitutive_model
{
    ($element: ident, $constitutive_model: ident, $constitutive_model_parameters: ident) =>
    {
        crate::fem::block::element::test::test_finite_element_with_viscoelastic_constitutive_model!(
            $element, $constitutive_model, $constitutive_model_parameters
        );
        #[test]
        fn hyperviscoelastic()
        {
            todo!()
        }
    }
}
pub(crate) use test_finite_element_with_hyperviscoelastic_constitutive_model;