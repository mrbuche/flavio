macro_rules! test_finite_element_block
{
    ($element: ident) =>
    {
        mod block
        {
            use crate::
            {
                EPSILON,
                fem::block::test::
                {
                    test_finite_element_block_with_elastic_constitutive_model,
                    test_finite_element_block_with_hyperelastic_constitutive_model,
                    test_finite_element_block_with_hyperviscoelastic_constitutive_model
                },
                math::
                {
                    Convert,
                    TensorRank2
                },
                mechanics::test::
                {
                    get_rotation_current_configuration,
                    get_rotation_rate_current_configuration,
                    get_rotation_reference_configuration,
                    get_translation_current_configuration,
                    get_translation_rate_current_configuration,
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
                    test_finite_element_block_with_elastic_constitutive_model!(
                        ElasticBlock, $element, AlmansiHamel, ALMANSIHAMELPARAMETERS
                    );
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
                    test_finite_element_block_with_hyperelastic_constitutive_model!(
                        ElasticBlock, $element, ArrudaBoyce, ARRUDABOYCEPARAMETERS
                    );
                }
                mod fung
                {
                    use super::*;
                    test_finite_element_block_with_hyperelastic_constitutive_model!(
                        ElasticBlock, $element, Fung, FUNGPARAMETERS
                    );
                }
                mod gent
                {
                    use super::*;
                    test_finite_element_block_with_hyperelastic_constitutive_model!(
                        ElasticBlock, $element, Gent, GENTPARAMETERS);
                }
                mod mooney_rivlin
                {
                    use super::*;
                    test_finite_element_block_with_hyperelastic_constitutive_model!(
                        ElasticBlock, $element, MooneyRivlin, MOONEYRIVLINPARAMETERS
                    );
                }
                mod neo_hookean
                {
                    use super::*;
                    test_finite_element_block_with_hyperelastic_constitutive_model!(
                        ElasticBlock, $element, NeoHookean, NEOHOOKEANPARAMETERS
                    );
                }
                mod saint_venant_kirchoff
                {
                    use super::*;
                    test_finite_element_block_with_hyperelastic_constitutive_model!(
                        ElasticBlock, $element, SaintVenantKirchoff, SAINTVENANTKIRCHOFFPARAMETERS
                    );
                }
                mod yeoh
                {
                    use super::*;
                    test_finite_element_block_with_hyperelastic_constitutive_model!(
                        ElasticBlock, $element, Yeoh, YEOHPARAMETERS
                    );
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
                    test_finite_element_block_with_hyperviscoelastic_constitutive_model!(
                        ViscoelasticBlock, $element, SaintVenantKirchoff, SAINTVENANTKIRCHOFFPARAMETERS
                    );
                }
            }
        }
    }
}
pub(crate) use test_finite_element_block;

macro_rules! test_nodal_forces_and_nodal_stiffnesses
{
    ($block: ident, $element: ident, $constitutive_model: ident, $constitutive_model_parameters: ident) =>
    {
        fn get_block<'a>() -> $block<D, E, $element<$constitutive_model<'a>>, G, N>
        {
            $block::new(
                $constitutive_model_parameters,
                get_connectivity(),
                get_reference_coordinates_block()
            )
        }
        fn get_block_transformed<'a>() -> $block<D, E, $element<$constitutive_model<'a>>, G, N>
        {
            $block::new(
                $constitutive_model_parameters,
                get_connectivity(),
                get_reference_coordinates_transformed_block()
            )
        }
        fn get_coordinates_transformed_block() -> NodalCoordinates<D>
        {
            get_coordinates_block().iter()
            .map(|coordinate|
                (get_rotation_current_configuration() * coordinate)
                + get_translation_current_configuration()
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
                    );
                }
                #[test]
                fn objectivity()
                {
                    get_nodal_forces(false, true).iter()
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
                    get_nodal_forces(false, false).iter()
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
                    );
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
                    );
                }
            }
        }
        #[test]
        fn size()
        {
            assert_eq!(
                std::mem::size_of::<ElasticBlock<D, E, $element::<$constitutive_model>, G, N>>(),
                std::mem::size_of::<Connectivity<E, N>>()
                + std::mem::size_of::<NodalCoordinates<D>>()
                + E * std::mem::size_of::<$element::<$constitutive_model>>()
            )
        }
    }
}
pub(crate) use test_nodal_forces_and_nodal_stiffnesses;

macro_rules! test_helmholtz_free_energy
{
    ($block: ident, $element: ident, $constitutive_model: ident, $constitutive_model_parameters: ident) =>
    {
        fn get_finite_difference_of_helmholtz_free_energy(is_deformed: bool) -> NodalForces<D>
        {
            let mut block = get_block();
            let mut finite_difference = 0.0;
            (0..D).map(|node|
                (0..3).map(|i|{
                    let mut nodal_coordinates = 
                    if is_deformed
                    {
                        get_coordinates_block()
                    }
                    else
                    {
                        get_reference_coordinates_block().convert()
                    };
                    nodal_coordinates[node][i] += 0.5 * EPSILON;
                    block.set_nodal_coordinates(nodal_coordinates);
                    finite_difference = block.calculate_helmholtz_free_energy();
                    let mut nodal_coordinates = 
                    if is_deformed
                    {
                        get_coordinates_block()
                    }
                    else
                    {
                        get_reference_coordinates_block().convert()
                    };
                    nodal_coordinates[node][i] -= 0.5 * EPSILON;
                    block.set_nodal_coordinates(nodal_coordinates);
                    finite_difference -= block.calculate_helmholtz_free_energy();
                    finite_difference/EPSILON
                }).collect()
            ).collect()
        }
        #[test]
        fn nodal_stiffnesses_deformed_symmetry()
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
                    block.set_nodal_coordinates(
                        get_coordinates_block()
                    );
                    block.calculate_nodal_forces().iter()
                    .zip(get_finite_difference_of_helmholtz_free_energy(true).iter())
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
                    block.set_nodal_coordinates(
                        get_coordinates_block()
                    );
                    let nodal_forces = block.calculate_nodal_forces();
                    let minimum = block.calculate_helmholtz_free_energy() - nodal_forces.dot(&get_coordinates_block());
                    let mut perturbed_coordinates = get_coordinates_block();
                    (0..D).for_each(|node|
                        (0..3).for_each(|i|{
                            perturbed_coordinates = get_coordinates_block();
                            perturbed_coordinates[node][i] += 0.5 * EPSILON;
                            block.set_nodal_coordinates(
                                &perturbed_coordinates * 1.0
                            );
                            assert!(
                                block.calculate_helmholtz_free_energy()
                                - nodal_forces.dot(
                                    &perturbed_coordinates
                                ) > minimum
                            );
                            perturbed_coordinates[node][i] -= EPSILON;
                            block.set_nodal_coordinates(
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
                    block_1.set_nodal_coordinates(
                        get_coordinates_block()
                    );
                    block_2.set_nodal_coordinates(
                        get_coordinates_transformed_block()
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
                    block.set_nodal_coordinates(
                        get_coordinates_block()
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
                    let mut block = get_block();
                    let minimum = block.calculate_helmholtz_free_energy();
                    let mut perturbed_coordinates = get_reference_coordinates_block();
                    (0..D).for_each(|node|
                        (0..3).for_each(|i|{
                            perturbed_coordinates = get_reference_coordinates_block();
                            perturbed_coordinates[node][i] += 0.5 * EPSILON;
                            block.set_nodal_coordinates(
                                perturbed_coordinates.convert()
                            );
                            assert!(
                                block.calculate_helmholtz_free_energy() > minimum
                            );
                            perturbed_coordinates[node][i] -= EPSILON;
                            block.set_nodal_coordinates(
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
                    block_1.set_nodal_coordinates(
                        get_reference_coordinates_block().convert()
                    );
                    block_2.set_nodal_coordinates(
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
                    block.set_nodal_coordinates(
                        get_reference_coordinates_block().convert()
                    );
                    assert_eq!(block.calculate_helmholtz_free_energy(), 0.0);
                }
            }
        }
    }
}
pub(crate) use test_helmholtz_free_energy;

macro_rules! test_finite_element_block_with_elastic_constitutive_model
{
    ($block: ident, $element: ident, $constitutive_model: ident, $constitutive_model_parameters: ident) =>
    {
        fn get_finite_difference_of_nodal_forces(is_deformed: bool) -> NodalStiffnesses<D>
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
                                get_coordinates_block()
                            }
                            else
                            {
                                get_reference_coordinates_block().convert()
                            };
                            nodal_coordinates[node_a][i] += 0.5 * EPSILON;
                            block.set_nodal_coordinates(nodal_coordinates);
                            finite_difference = block.calculate_nodal_forces()[node_b][j];
                            let mut nodal_coordinates = 
                            if is_deformed
                            {
                                get_coordinates_block()
                            }
                            else
                            {
                                get_reference_coordinates_block().convert()
                            };
                            nodal_coordinates[node_a][i] -= 0.5 * EPSILON;
                            block.set_nodal_coordinates(nodal_coordinates);
                            finite_difference -= block.calculate_nodal_forces()[node_b][j];
                            finite_difference/EPSILON
                        }).collect()
                    ).collect()
                ).collect()
            ).collect()
        }
        fn get_nodal_forces(is_deformed: bool, is_rotated: bool) -> NodalForces<D>
        {
            if is_rotated
            {
                if is_deformed
                {
                    let mut block = get_block_transformed();
                    block.set_nodal_coordinates(
                        get_coordinates_transformed_block()
                    );
                    get_rotation_current_configuration().transpose() *
                    block.calculate_nodal_forces()
                }
                else
                {
                    let converted: TensorRank2<3, 1, 1> = get_rotation_reference_configuration().convert();
                    converted.transpose() *
                    get_block_transformed().calculate_nodal_forces()
                }
            }
            else
            {
                if is_deformed
                {
                    let mut block = get_block();
                    block.set_nodal_coordinates(
                        get_coordinates_block()
                    );
                    block.calculate_nodal_forces()
                }
                else
                {
                    get_block().calculate_nodal_forces()
                }
            }
        }
        fn get_nodal_stiffnesses(is_deformed: bool, is_rotated: bool) -> NodalStiffnesses<D>
        {
            if is_rotated
            {
                if is_deformed
                {
                    let mut block = get_block_transformed();
                    block.set_nodal_coordinates(
                        get_coordinates_transformed_block()
                    );
                    get_rotation_current_configuration().transpose() *
                    block.calculate_nodal_stiffnesses()
                    * get_rotation_current_configuration()
                }
                else
                {
                    let converted: TensorRank2<3, 1, 1> = get_rotation_reference_configuration().convert();
                    converted.transpose() *
                    get_block_transformed().calculate_nodal_stiffnesses()
                    * converted
                }
            }
            else
            {
                if is_deformed
                {
                    let mut block = get_block();
                    block.set_nodal_coordinates(
                        get_coordinates_block()
                    );
                    block.calculate_nodal_stiffnesses()
                }
                else
                {
                    get_block().calculate_nodal_stiffnesses()
                }
            }
        }
        crate::fem::block::test::test_nodal_forces_and_nodal_stiffnesses!(
            $block, $element, $constitutive_model, $constitutive_model_parameters
        );
    }
}
pub(crate) use test_finite_element_block_with_elastic_constitutive_model;

macro_rules! test_finite_element_block_with_hyperelastic_constitutive_model
{
    ($block: ident, $element: ident, $constitutive_model: ident, $constitutive_model_parameters: ident) =>
    {
        crate::fem::block::test::test_finite_element_block_with_elastic_constitutive_model!(
            $block, $element, $constitutive_model, $constitutive_model_parameters
        );
        crate::fem::block::test::test_helmholtz_free_energy!(
            $block, $element, $constitutive_model, $constitutive_model_parameters
        );
    }
}
pub(crate) use test_finite_element_block_with_hyperelastic_constitutive_model;

macro_rules! test_finite_element_block_with_viscoelastic_constitutive_model
{
    ($block: ident, $element: ident, $constitutive_model: ident, $constitutive_model_parameters: ident) =>
    {
        fn get_finite_difference_of_nodal_forces(is_deformed: bool) -> NodalStiffnesses<D>
        {
            let mut block = get_block();
            if is_deformed
            {
                block.set_nodal_coordinates(get_coordinates_block());
            }
            else
            {
                block.set_nodal_coordinates(get_reference_coordinates_block().convert());
            }
            let mut finite_difference = 0.0;
            (0..D).map(|node_a|
                (0..D).map(|node_b|
                    (0..3).map(|i|
                        (0..3).map(|j|{
                            let mut nodal_velocities = 
                            if is_deformed
                            {
                                get_velocities_block()
                            }
                            else
                            {
                                NodalVelocities::zero()
                            };
                            nodal_velocities[node_a][i] += 0.5 * EPSILON;
                            block.set_nodal_velocities(nodal_velocities);
                            finite_difference = block.calculate_nodal_forces()[node_b][j];
                            let mut nodal_velocities = 
                            if is_deformed
                            {
                                get_velocities_block()
                            }
                            else
                            {
                                NodalVelocities::zero()
                            };
                            nodal_velocities[node_a][i] -= 0.5 * EPSILON;
                            block.set_nodal_velocities(nodal_velocities);
                            finite_difference -= block.calculate_nodal_forces()[node_b][j];
                            finite_difference/EPSILON
                        }).collect()
                    ).collect()
                ).collect()
            ).collect()
        }
        crate::fem::block::test::test_nodal_forces_and_nodal_stiffnesses!(
            $block, $element, $constitutive_model, $constitutive_model_parameters
        );
    }
}
pub(crate) use test_finite_element_block_with_viscoelastic_constitutive_model;

macro_rules! test_finite_element_block_with_hyperviscoelastic_constitutive_model
{
    ($block: ident, $element: ident, $constitutive_model: ident, $constitutive_model_parameters: ident) =>
    {
        fn get_velocities_transformed_block() -> NodalCoordinates<D>
        {
            get_coordinates_block().iter()
            .zip(get_velocities_block().iter())
            .map(|(coordinate, velocity)|
                get_rotation_current_configuration() * velocity
                + get_rotation_rate_current_configuration() * coordinate
                + get_translation_rate_current_configuration()
            ).collect()
        }
        crate::fem::block::test::test_finite_element_block_with_viscoelastic_constitutive_model!(
            $block, $element, $constitutive_model, $constitutive_model_parameters
        );
        fn get_finite_difference_of_viscous_dissipation(is_deformed: bool) -> NodalForces<D>
        {
            let mut block = get_block();
            if is_deformed
            {
                block.set_nodal_coordinates(get_coordinates_block());
            }
            else
            {
                block.set_nodal_coordinates(get_reference_coordinates_block().convert());
            }
            let mut finite_difference = 0.0;
            (0..D).map(|node|
                (0..3).map(|i|{
                    let mut nodal_velocities = 
                    if is_deformed
                    {
                        get_velocities_block()
                    }
                    else
                    {
                        NodalVelocities::zero()
                    };
                    nodal_velocities[node][i] += 0.5 * EPSILON;
                    block.set_nodal_velocities(nodal_velocities);
                    finite_difference = block.calculate_viscous_dissipation();
                    let mut nodal_velocities = 
                    if is_deformed
                    {
                        get_velocities_block()
                    }
                    else
                    {
                        NodalVelocities::zero()
                    };
                    nodal_velocities[node][i] -= 0.5 * EPSILON;
                    block.set_nodal_velocities(nodal_velocities);
                    finite_difference -= block.calculate_viscous_dissipation();
                    finite_difference/EPSILON
                }).collect()
            ).collect()
        }
        fn get_finite_difference_of_dissipation_potential(is_deformed: bool) -> NodalForces<D>
        {
            let mut block = get_block();
            if is_deformed
            {
                block.set_nodal_coordinates(get_coordinates_block());
            }
            else
            {
                block.set_nodal_coordinates(get_reference_coordinates_block().convert());
            }
            let mut finite_difference = 0.0;
            (0..D).map(|node|
                (0..3).map(|i|{
                    let mut nodal_velocities = 
                    if is_deformed
                    {
                        get_velocities_block()
                    }
                    else
                    {
                        NodalVelocities::zero()
                    };
                    nodal_velocities[node][i] += 0.5 * EPSILON;
                    block.set_nodal_velocities(nodal_velocities);
                    finite_difference = block.calculate_dissipation_potential();
                    let mut nodal_velocities = 
                    if is_deformed
                    {
                        get_velocities_block()
                    }
                    else
                    {
                        NodalVelocities::zero()
                    };
                    nodal_velocities[node][i] -= 0.5 * EPSILON;
                    block.set_nodal_velocities(nodal_velocities);
                    finite_difference -= block.calculate_dissipation_potential();
                    finite_difference/EPSILON
                }).collect()
            ).collect()
        }
        fn get_nodal_forces(is_deformed: bool, is_rotated: bool) -> NodalForces<D>
        {
            if is_rotated
            {
                if is_deformed
                {
                    let mut block = get_block_transformed();
                    block.set_nodal_coordinates(
                        get_coordinates_transformed_block()
                    );
                    block.set_nodal_velocities(
                        get_velocities_transformed_block()
                    );
                    get_rotation_current_configuration().transpose() *
                    block.calculate_nodal_forces()
                }
                else
                {
                    let converted: TensorRank2<3, 1, 1> = get_rotation_reference_configuration().convert();
                    converted.transpose() *
                    get_block_transformed().calculate_nodal_forces()
                }
            }
            else
            {
                if is_deformed
                {
                    let mut block = get_block();
                    block.set_nodal_coordinates(
                        get_coordinates_block()
                    );
                    block.set_nodal_velocities(
                        get_velocities_block()
                    );
                    block.calculate_nodal_forces()
                }
                else
                {
                    get_block().calculate_nodal_forces()
                }
            }
        }
        fn get_nodal_stiffnesses(is_deformed: bool, is_rotated: bool) -> NodalStiffnesses<D>
        {
            if is_rotated
            {
                if is_deformed
                {
                    let mut block = get_block_transformed();
                    block.set_nodal_coordinates(
                        get_coordinates_transformed_block()
                    );
                    block.set_nodal_velocities(
                        get_velocities_transformed_block()
                    );
                    get_rotation_current_configuration().transpose() *
                    block.calculate_nodal_stiffnesses()
                    * get_rotation_current_configuration()
                }
                else
                {
                    let converted: TensorRank2<3, 1, 1> = get_rotation_reference_configuration().convert();
                    converted.transpose() *
                    get_block_transformed().calculate_nodal_stiffnesses()
                    * converted
                }
            }
            else
            {
                if is_deformed
                {
                    let mut block = get_block();
                    block.set_nodal_coordinates(
                        get_coordinates_block()
                    );
                    block.set_nodal_velocities(
                        get_velocities_block()
                    );
                    block.calculate_nodal_stiffnesses()
                }
                else
                {
                    get_block().calculate_nodal_stiffnesses()
                }
            }
        }
        crate::fem::block::test::test_helmholtz_free_energy!(
            $block, $element, $constitutive_model, $constitutive_model_parameters
        );
        mod viscous_dissipation
        {
            use super::*;
            mod deformed
            {
                use super::*;
                #[test]
                fn finite_difference()
                {
                    let mut block = get_block();
                    block.set_nodal_coordinates(
                        get_coordinates_block()
                    );
                    let nodal_forces_0 = block.calculate_nodal_forces();
                    block.set_nodal_velocities(
                        get_velocities_block()
                    );
                    (block.calculate_nodal_forces() - nodal_forces_0).iter()
                    .zip(get_finite_difference_of_viscous_dissipation(true).iter())
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
                    block.set_nodal_coordinates(
                        get_coordinates_block()
                    );
                    let nodal_forces_0 = block.calculate_nodal_forces();
                    block.set_nodal_velocities(
                        get_velocities_block()
                    );
                    let nodal_forces = block.calculate_nodal_forces() - nodal_forces_0;
                    let minimum = block.calculate_viscous_dissipation() - nodal_forces.dot(&get_velocities_block());
                    let mut perturbed_velocities = get_velocities_block();
                    (0..D).for_each(|node|
                        (0..3).for_each(|i|{
                            perturbed_velocities = get_velocities_block();
                            perturbed_velocities[node][i] += 0.5 * EPSILON;
                            block.set_nodal_velocities(
                                &perturbed_velocities * 1.0
                            );
                            assert!(
                                block.calculate_viscous_dissipation()
                                - nodal_forces.dot(
                                    &perturbed_velocities
                                ) > minimum
                            );
                            perturbed_velocities[node][i] -= EPSILON;
                            block.set_nodal_velocities(
                                &perturbed_velocities * 1.0
                            );
                            assert!(
                                block.calculate_viscous_dissipation()
                                - nodal_forces.dot(
                                    &perturbed_velocities
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
                    block_1.set_nodal_coordinates(
                        get_coordinates_block()
                    );
                    block_1.set_nodal_velocities(
                        get_velocities_block()
                    );
                    block_2.set_nodal_coordinates(
                        get_coordinates_transformed_block()
                    );
                    block_2.set_nodal_velocities(
                        get_velocities_transformed_block()
                    );
                    assert_eq_within_tols(
                        &block_1.calculate_viscous_dissipation(),
                        &block_2.calculate_viscous_dissipation()
                    );
                }
                #[test]
                fn positive()
                {
                    let mut block = get_block();
                    assert_eq!(block.calculate_viscous_dissipation(), 0.0);
                    block.set_nodal_coordinates(
                        get_coordinates_block()
                    );
                    block.set_nodal_velocities(
                        get_velocities_block()
                    );
                    assert!(block.calculate_viscous_dissipation() > 0.0);
                }
            }
            mod undeformed
            {
                use super::*;
                #[test]
                fn finite_difference()
                {
                    get_finite_difference_of_viscous_dissipation(
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
                    let minimum = block.calculate_viscous_dissipation();
                    let mut perturbed_velocities = NodalVelocities::zero();
                    (0..D).for_each(|node|
                        (0..3).for_each(|i|{
                            perturbed_velocities = NodalVelocities::zero();
                            perturbed_velocities[node][i] += 0.5 * EPSILON;
                            block.set_nodal_velocities(
                                perturbed_velocities.convert()
                            );
                            assert!(
                                block.calculate_viscous_dissipation() > minimum
                            );
                            perturbed_velocities[node][i] -= EPSILON;
                            block.set_nodal_velocities(
                                perturbed_velocities.convert()
                            );
                            assert!(
                                block.calculate_viscous_dissipation() > minimum
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
                        &block_1.calculate_viscous_dissipation(),
                        &block_2.calculate_viscous_dissipation()
                    );
                    block_1.set_nodal_coordinates(
                        get_reference_coordinates_block().convert()
                    );
                    block_1.set_nodal_velocities(
                        NodalVelocities::zero()
                    );
                    block_2.set_nodal_coordinates(
                        get_reference_coordinates_transformed_block().convert()
                    );
                    block_2.set_nodal_velocities(
                        NodalVelocities::zero()
                    );
                    assert_eq_within_tols(
                        &block_1.calculate_viscous_dissipation(),
                        &block_2.calculate_viscous_dissipation()
                    );
                }
                #[test]
                fn zero()
                {
                    let mut block = get_block();
                    assert_eq!(block.calculate_viscous_dissipation(), 0.0);
                    block.set_nodal_coordinates(
                        get_reference_coordinates_block().convert()
                    );
                    assert_eq!(block.calculate_viscous_dissipation(), 0.0);
                }
            }
        }
        mod dissipation_potential
        {
            use super::*;
            mod deformed
            {
                use super::*;
                #[test]
                fn finite_difference()
                {
                    let mut block = get_block();
                    block.set_nodal_coordinates(
                        get_coordinates_block()
                    );
                    block.set_nodal_velocities(
                        get_velocities_block()
                    );
                    block.calculate_nodal_forces().iter()
                    .zip(get_finite_difference_of_dissipation_potential(true).iter())
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
                    block.set_nodal_coordinates(
                        get_coordinates_block()
                    );
                    block.set_nodal_velocities(
                        get_velocities_block()
                    );
                    let nodal_forces = block.calculate_nodal_forces();
                    let minimum = block.calculate_dissipation_potential() - nodal_forces.dot(&get_velocities_block());
                    let mut perturbed_velocities = get_velocities_block();
                    (0..D).for_each(|node|
                        (0..3).for_each(|i|{
                            perturbed_velocities = get_velocities_block();
                            perturbed_velocities[node][i] += 0.5 * EPSILON;
                            block.set_nodal_velocities(
                                &perturbed_velocities * 1.0
                            );
                            assert!(
                                block.calculate_dissipation_potential()
                                - nodal_forces.dot(
                                    &perturbed_velocities
                                ) > minimum
                            );
                            perturbed_velocities[node][i] -= EPSILON;
                            block.set_nodal_velocities(
                                &perturbed_velocities * 1.0
                            );
                            assert!(
                                block.calculate_dissipation_potential()
                                - nodal_forces.dot(
                                    &perturbed_velocities
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
                    block_1.set_nodal_coordinates(
                        get_coordinates_block()
                    );
                    block_1.set_nodal_velocities(
                        get_velocities_block()
                    );
                    block_2.set_nodal_coordinates(
                        get_coordinates_transformed_block()
                    );
                    block_2.set_nodal_velocities(
                        get_velocities_transformed_block()
                    );
                    assert_eq_within_tols(
                        &block_1.calculate_dissipation_potential(),
                        &block_2.calculate_dissipation_potential()
                    );
                }
                #[test]
                fn positive()
                {
                    let mut block = get_block();
                    assert_eq!(block.calculate_dissipation_potential(), 0.0);
                    block.set_nodal_coordinates(
                        get_coordinates_block()
                    );
                    block.set_nodal_velocities(
                        get_velocities_block()
                    );
                    assert!(block.calculate_dissipation_potential() > 0.0);
                }
            }
            mod undeformed
            {
                use super::*;
                #[test]
                fn finite_difference()
                {
                    get_finite_difference_of_dissipation_potential(
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
                    let minimum = block.calculate_dissipation_potential();
                    let mut perturbed_velocities = NodalVelocities::zero();
                    (0..D).for_each(|node|
                        (0..3).for_each(|i|{
                            perturbed_velocities = NodalVelocities::zero();
                            perturbed_velocities[node][i] += 0.5 * EPSILON;
                            block.set_nodal_velocities(
                                perturbed_velocities.convert()
                            );
                            assert!(
                                block.calculate_dissipation_potential() > minimum
                            );
                            perturbed_velocities[node][i] -= EPSILON;
                            block.set_nodal_velocities(
                                perturbed_velocities.convert()
                            );
                            assert!(
                                block.calculate_dissipation_potential() > minimum
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
                        &block_1.calculate_dissipation_potential(),
                        &block_2.calculate_dissipation_potential()
                    );
                    block_1.set_nodal_coordinates(
                        get_reference_coordinates_block().convert()
                    );
                    block_1.set_nodal_velocities(
                        NodalVelocities::zero()
                    );
                    block_2.set_nodal_coordinates(
                        get_reference_coordinates_transformed_block().convert()
                    );
                    block_2.set_nodal_velocities(
                        NodalVelocities::zero()
                    );
                    assert_eq_within_tols(
                        &block_1.calculate_dissipation_potential(),
                        &block_2.calculate_dissipation_potential()
                    );
                }
                #[test]
                fn zero()
                {
                    let mut block = get_block();
                    assert_eq!(block.calculate_dissipation_potential(), 0.0);
                    block.set_nodal_coordinates(
                        get_reference_coordinates_block().convert()
                    );
                    assert_eq!(block.calculate_dissipation_potential(), 0.0);
                }
            }
        }
    }
}
pub(crate) use test_finite_element_block_with_hyperviscoelastic_constitutive_model;