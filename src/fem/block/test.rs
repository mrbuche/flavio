macro_rules! test_finite_element_block
{
    ($element: ident) =>
    {
        mod finite_element_block
        {
            use crate::
            {
                EPSILON,
                constitutive::
                {
                    hyperelastic::
                    {
                        ArrudaBoyceModel,
                        GentModel,
                        MooneyRivlinModel,
                        NeoHookeanModel,
                        YeohModel,
                    },
                    test::
                    {
                        ARRUDABOYCEPARAMETERS,
                        GENTPARAMETERS,
                        MOONEYRIVLINPARAMETERS,
                        NEOHOOKEANPARAMETERS,
                        YEOHPARAMETERS
                    }
                },
                fem::block::test::test_finite_element_block_with_constitutive_model,
                math::
                {
                    Convert,
                    TensorRank2
                },
                mechanics::test::
                {
                    get_rotation_current_configuration,
                    get_rotation_reference_configuration,
                    get_translation_current_configuration,
                    get_translation_reference_configuration
                },
                test::assert_eq_within_tols
            };
            use super::*;
            // pub mod almansi_hamel
            // {
            //     use super::*;
            //     test_finite_element_block_with_constitutive_model!($element, AlmansiHamelModel, ALMANSIHAMELPARAMETERS);
            // }
            pub mod arruda_boyce
            {
                use super::*;
                test_finite_element_block_with_constitutive_model!($element, ArrudaBoyceModel, ARRUDABOYCEPARAMETERS);
            }
            pub mod gent
            {
                use super::*;
                test_finite_element_block_with_constitutive_model!($element, GentModel, GENTPARAMETERS);
            }
            pub mod mooney_rivlin
            {
                use super::*;
                test_finite_element_block_with_constitutive_model!($element, MooneyRivlinModel, MOONEYRIVLINPARAMETERS);
            }
            pub mod neo_hookean
            {
                use super::*;
                test_finite_element_block_with_constitutive_model!($element, NeoHookeanModel, NEOHOOKEANPARAMETERS);
            }
            pub mod yeoh
            {
                use super::*;
                test_finite_element_block_with_constitutive_model!($element, YeohModel, YEOHPARAMETERS);
            }
        }
    }
}
pub(crate) use test_finite_element_block;
macro_rules! test_finite_element_block_with_constitutive_model
{
    ($element: ident, $constitutive_model: ident, $constitutive_model_parameters: ident) =>
    {
        fn get_block<'a>() -> FiniteElementBlock<'a, $constitutive_model<'a>, D, E, $element<'a, $constitutive_model<'a>>, G, N>
        {
            FiniteElementBlock::new(
                $constitutive_model_parameters,
                get_connectivity(),
                get_reference_coordinates_block()
            )
        }
        fn get_block_transformed<'a>() -> FiniteElementBlock<'a, $constitutive_model<'a>, D, E, $element<'a, $constitutive_model<'a>>, G, N>
        {
            FiniteElementBlock::new(
                $constitutive_model_parameters,
                get_connectivity(),
                get_reference_coordinates_transformed_block()
            )
        }
        fn get_current_coordinates_transformed_block() -> CurrentNodalCoordinates<D>
        {
            get_current_coordinates_block().iter()
            .map(|current_coordinate|
                (get_rotation_current_configuration() * current_coordinate)
                + get_translation_current_configuration()
            ).collect()
        }
        fn get_fd_helmholtz_free_energy(is_deformed: bool) -> NodalForces<D>
        {
            let mut block = get_block();
            let mut finite_difference = 0.0;
            (0..D).map(|node|
                (0..3).map(|i|{
                    let mut nodal_coordinates = 
                    if is_deformed
                    {
                        get_current_coordinates_block()
                    }
                    else
                    {
                        get_reference_coordinates_block().convert()
                    };
                    nodal_coordinates[node][i] += 0.5 * EPSILON;
                    block.set_current_nodal_coordinates(nodal_coordinates);
                    finite_difference = block.calculate_helmholtz_free_energy();
                    let mut nodal_coordinates = 
                    if is_deformed
                    {
                        get_current_coordinates_block()
                    }
                    else
                    {
                        get_reference_coordinates_block().convert()
                    };
                    nodal_coordinates[node][i] -= 0.5 * EPSILON;
                    block.set_current_nodal_coordinates(nodal_coordinates);
                    finite_difference -= block.calculate_helmholtz_free_energy();
                    finite_difference/EPSILON
                }).collect()
            ).collect()
        }
        fn get_fd_nodal_forces(is_deformed: bool) -> NodalStiffnesses<D>
        {
            let mut block = get_block();
            let mut finite_difference = 0.0;
            (0..D).map(|node_a|
                (0..D).map(|node_b|
                    (0..3).map(|i|
                        (0..3).map(|j|{
                            let mut nodal_coordinates = 
                            if is_deformed
                            {
                                get_current_coordinates_block()
                            }
                            else
                            {
                                get_reference_coordinates_block().convert()
                            };
                            nodal_coordinates[node_a][i] += 0.5 * EPSILON;
                            block.set_current_nodal_coordinates(nodal_coordinates);
                            finite_difference = block.calculate_nodal_forces()[node_b][j];
                            let mut nodal_coordinates = 
                            if is_deformed
                            {
                                get_current_coordinates_block()
                            }
                            else
                            {
                                get_reference_coordinates_block().convert()
                            };
                            nodal_coordinates[node_a][i] -= 0.5 * EPSILON;
                            block.set_current_nodal_coordinates(nodal_coordinates);
                            finite_difference -= block.calculate_nodal_forces()[node_b][j];
                            finite_difference/EPSILON
                        }).collect()
                    ).collect()
                ).collect()
            ).collect()
        }
        fn get_reference_coordinates_transformed_block() -> ReferenceNodalCoordinates<D>
        {
            get_reference_coordinates_block().iter()
            .map(|reference_coordinate|
                (get_rotation_reference_configuration() * reference_coordinate)
                + get_translation_reference_configuration()
            ).collect()
        }
        #[test]
        fn temporary()
        {
            let mut block = get_block();
            block.set_current_nodal_coordinates(
                get_current_coordinates_block()
            );
            block.solve_using_gradient_descent();
            block.calculate_nodal_forces().iter()
            .for_each(|nodal_force|
                nodal_force.iter()
                .for_each(|nodal_force_i|
                    assert_eq_within_tols(
                        nodal_force_i, &0.0
                    )
                )
            );
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
                    let mut block = get_block();
                    block.set_current_nodal_coordinates(
                        get_current_coordinates_block()
                    );
                    block.calculate_nodal_forces().iter()
                    .zip(get_fd_helmholtz_free_energy(true).iter())
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
                    let mut block = get_block();
                    block.set_current_nodal_coordinates(
                        get_current_coordinates_block()
                    );
                    let nodal_forces = block.calculate_nodal_forces();
                    let minimum = block.calculate_helmholtz_free_energy()
                        - nodal_forces.dot(
                        &get_current_coordinates_block()
                    );
                    let mut perturbed_coordinates = get_current_coordinates_block();
                    (0..D).for_each(|node|
                        (0..3).for_each(|i|{
                            perturbed_coordinates = get_current_coordinates_block();
                            perturbed_coordinates[node][i] += 0.5 * EPSILON;
                            block.set_current_nodal_coordinates(
                                &perturbed_coordinates * 1.0
                            );
                            assert!(
                                block.calculate_helmholtz_free_energy()
                                - nodal_forces.dot(
                                    &perturbed_coordinates
                                ) > minimum
                            );
                            perturbed_coordinates[node][i] -= EPSILON;
                            block.set_current_nodal_coordinates(
                                &perturbed_coordinates * 1.0
                            );
                            assert!(
                                block.calculate_helmholtz_free_energy()
                                - nodal_forces.dot(
                                    &perturbed_coordinates
                                ) > minimum
                            );
                        })
                    )
                }
                #[test]
                fn objectivity()
                {
                    let mut block_1 = get_block();
                    let mut block_2 = get_block_transformed();
                    block_1.set_current_nodal_coordinates(
                        get_current_coordinates_block()
                    );
                    block_2.set_current_nodal_coordinates(
                        get_current_coordinates_transformed_block()
                    );
                    assert_eq_within_tols(
                        &block_1.calculate_helmholtz_free_energy(),
                        &block_2.calculate_helmholtz_free_energy()
                    );
                }
                #[test]
                fn positive()
                {
                    let mut block = get_block();
                    assert_eq!(block.calculate_helmholtz_free_energy(), 0.0);
                    block.set_current_nodal_coordinates(
                        get_current_coordinates_block()
                    );
                    assert!(block.calculate_helmholtz_free_energy() > 0.0);
                }
            }
            mod undeformed
            {
                use super::*;
                #[test]
                fn finite_difference()
                {
                    get_fd_helmholtz_free_energy(
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
                    let mut block = get_block();
                    let minimum = block.calculate_helmholtz_free_energy();
                    let mut perturbed_coordinates = get_reference_coordinates_block();
                    (0..D).for_each(|node|
                        (0..3).for_each(|i|{
                            perturbed_coordinates = get_reference_coordinates_block();
                            perturbed_coordinates[node][i] += 0.5 * EPSILON;
                            block.set_current_nodal_coordinates(
                                perturbed_coordinates.convert()
                            );
                            assert!(
                                block.calculate_helmholtz_free_energy() > minimum
                            );
                            perturbed_coordinates[node][i] -= EPSILON;
                            block.set_current_nodal_coordinates(
                                perturbed_coordinates.convert()
                            );
                            assert!(
                                block.calculate_helmholtz_free_energy() > minimum
                            );
                        })
                    )
                }
                #[test]
                fn objectivity()
                {
                    let mut block_1 = get_block();
                    let mut block_2 = get_block_transformed();
                    assert_eq_within_tols(
                        &block_1.calculate_helmholtz_free_energy(),
                        &block_2.calculate_helmholtz_free_energy()
                    );
                    block_1.set_current_nodal_coordinates(
                        get_reference_coordinates_block().convert()
                    );
                    block_2.set_current_nodal_coordinates(
                        get_reference_coordinates_transformed_block().convert()
                    );
                    assert_eq_within_tols(
                        &block_1.calculate_helmholtz_free_energy(),
                        &block_2.calculate_helmholtz_free_energy()
                    );
                }
                #[test]
                fn zero()
                {
                    let mut block = get_block();
                    assert_eq!(block.calculate_helmholtz_free_energy(), 0.0);
                    block.set_current_nodal_coordinates(
                        get_reference_coordinates_block().convert()
                    );
                    assert_eq!(block.calculate_helmholtz_free_energy(), 0.0);
                }
            }
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
                    let mut block = get_block();
                    block.set_current_nodal_coordinates(
                        get_current_coordinates_block()
                    );
                    block.calculate_nodal_stiffnesses().iter()
                    .zip(get_fd_nodal_forces(true).iter())
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
                    let mut block_1 = get_block();
                    let mut block_2 = get_block_transformed();
                    block_1.set_current_nodal_coordinates(
                        get_current_coordinates_block()
                    );
                    block_2.set_current_nodal_coordinates(
                        get_current_coordinates_transformed_block()
                    );
                    block_1.calculate_nodal_forces().iter()
                    .zip((
                        get_rotation_current_configuration().transpose() *
                        block_2.calculate_nodal_forces()
                    ).iter()
                    ).for_each(|(nodal_force, res_nodal_force)|
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
                    let mut block = get_block();
                    block.calculate_nodal_stiffnesses().iter()
                    .zip(get_fd_nodal_forces(false).iter())
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
                    );
                    block.set_current_nodal_coordinates(
                        get_reference_coordinates_block().convert()
                    );
                    block.calculate_nodal_stiffnesses().iter()
                    .zip(get_fd_nodal_forces(false).iter())
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
                    );
                }
                #[test]
                fn objectivity()
                {
                    let mut block = get_block_transformed();
                    block.calculate_nodal_forces().iter()
                    .for_each(|nodal_force|
                        nodal_force.iter()
                        .for_each(|nodal_force_i|
                            assert_eq_within_tols(
                                nodal_force_i, &0.0
                            )
                        )
                    );
                    block.set_current_nodal_coordinates(
                        get_reference_coordinates_transformed_block().convert()
                    );
                    block.calculate_nodal_forces().iter()
                    .for_each(|nodal_force|
                        nodal_force.iter()
                        .for_each(|nodal_force_i|
                            assert_eq_within_tols(
                                nodal_force_i, &0.0
                            )
                        )
                    );
                }
                #[test]
                fn zero()
                {
                    let mut block = get_block();
                    block.calculate_nodal_forces().iter()
                    .for_each(|nodal_force|
                        nodal_force.iter().for_each(|nodal_force_i|
                            assert_eq_within_tols(
                                nodal_force_i, &0.0
                            )
                        )
                    );
                    block.set_current_nodal_coordinates(
                        get_reference_coordinates_block().convert()
                    );
                    block.calculate_nodal_forces().iter()
                    .for_each(|nodal_force|
                        nodal_force.iter().for_each(|nodal_force_i|
                            assert_eq_within_tols(
                                nodal_force_i, &0.0
                            )
                        )
                    );
                }
            }
        }
        mod nodal_stiffness
        {
            use super::*;
            mod deformed
            {
                use super::*;
                #[test]
                fn objectivity()
                {
                    let mut block_1 = get_block();
                    let mut block_2 = get_block_transformed();
                    block_1.set_current_nodal_coordinates(
                        get_current_coordinates_block()
                    );
                    block_2.set_current_nodal_coordinates(
                        get_current_coordinates_transformed_block()
                    );
                    block_1.calculate_nodal_stiffnesses().iter()
                    .zip((
                        get_rotation_current_configuration().transpose() *
                        block_2.calculate_nodal_stiffnesses()
                        * get_rotation_current_configuration()
                    ).iter())
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
                    let mut block = get_block();
                    block.set_current_nodal_coordinates(
                        get_current_coordinates_block()
                    );
                    let nodal_stiffness = block.calculate_nodal_stiffnesses();
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
                    let converted: TensorRank2<3, 1, 1> = get_rotation_reference_configuration().convert();
                    let mut block_1 = get_block();
                    let mut block_2 = get_block_transformed();
                    block_1.calculate_nodal_stiffnesses().iter()
                    .zip((
                        converted.transpose() *
                        block_2.calculate_nodal_stiffnesses()
                        * converted
                    ).iter())
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
                    );
                    let converted: TensorRank2<3, 1, 1> = get_rotation_reference_configuration().convert();
                    block_1.set_current_nodal_coordinates(
                        get_reference_coordinates_block().convert()
                    );
                    block_2.set_current_nodal_coordinates(
                        get_reference_coordinates_transformed_block().convert()
                    );
                    block_1.calculate_nodal_stiffnesses().iter()
                    .zip((
                        converted.transpose() *
                        block_2.calculate_nodal_stiffnesses()
                        * converted
                    ).iter())
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
                    );
                }
                #[test]
                fn symmetry()
                {
                    let mut block = get_block();
                    let nodal_stiffness = block.calculate_nodal_stiffnesses();
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
                    );
                    block.set_current_nodal_coordinates(
                        get_reference_coordinates_block().convert()
                    );
                    let nodal_stiffness = block.calculate_nodal_stiffnesses();
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
                    );
                }
            }
        }
    }
}
pub(crate) use test_finite_element_block_with_constitutive_model;