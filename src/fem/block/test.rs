macro_rules! test_finite_element_block {
    ($element: ident) => {
        mod block {
            use super::*;
            use crate::{
                fem::block::test::{
                    test_finite_element_block_with_elastic_constitutive_model,
                    test_finite_element_block_with_elastic_hyperviscous_constitutive_model,
                    test_finite_element_block_with_hyperelastic_constitutive_model,
                    test_finite_element_block_with_hyperviscoelastic_constitutive_model,
                },
                math::{
                    test::{assert_eq, assert_eq_from_fd, assert_eq_within_tols, TestError},
                    Convert, TensorRank2,
                },
                mechanics::test::{
                    get_rotation_current_configuration, get_rotation_rate_current_configuration,
                    get_rotation_reference_configuration, get_translation_current_configuration,
                    get_translation_rate_current_configuration,
                    get_translation_reference_configuration,
                },
                EPSILON,
            };
            mod elastic {
                use super::*;
                use crate::constitutive::solid::elastic::{
                    test::ALMANSIHAMELPARAMETERS, AlmansiHamel,
                };
                mod almansi_hamel {
                    use super::*;
                    test_finite_element_block_with_elastic_constitutive_model!(
                        ElasticBlock,
                        $element,
                        AlmansiHamel,
                        ALMANSIHAMELPARAMETERS
                    );
                }
            }
            mod hyperelastic {
                use super::*;
                use crate::constitutive::solid::hyperelastic::{
                    test::{
                        ARRUDABOYCEPARAMETERS, FUNGPARAMETERS, GENTPARAMETERS,
                        MOONEYRIVLINPARAMETERS, NEOHOOKEANPARAMETERS,
                        SAINTVENANTKIRCHOFFPARAMETERS, YEOHPARAMETERS,
                    },
                    ArrudaBoyce, Fung, Gent, MooneyRivlin, NeoHookean, SaintVenantKirchoff, Yeoh,
                };
                mod arruda_boyce {
                    use super::*;
                    test_finite_element_block_with_hyperelastic_constitutive_model!(
                        ElasticBlock,
                        $element,
                        ArrudaBoyce,
                        ARRUDABOYCEPARAMETERS
                    );
                }
                mod fung {
                    use super::*;
                    test_finite_element_block_with_hyperelastic_constitutive_model!(
                        ElasticBlock,
                        $element,
                        Fung,
                        FUNGPARAMETERS
                    );
                }
                mod gent {
                    use super::*;
                    test_finite_element_block_with_hyperelastic_constitutive_model!(
                        ElasticBlock,
                        $element,
                        Gent,
                        GENTPARAMETERS
                    );
                }
                mod mooney_rivlin {
                    use super::*;
                    test_finite_element_block_with_hyperelastic_constitutive_model!(
                        ElasticBlock,
                        $element,
                        MooneyRivlin,
                        MOONEYRIVLINPARAMETERS
                    );
                }
                mod neo_hookean {
                    use super::*;
                    test_finite_element_block_with_hyperelastic_constitutive_model!(
                        ElasticBlock,
                        $element,
                        NeoHookean,
                        NEOHOOKEANPARAMETERS
                    );
                }
                mod saint_venant_kirchoff {
                    use super::*;
                    test_finite_element_block_with_hyperelastic_constitutive_model!(
                        ElasticBlock,
                        $element,
                        SaintVenantKirchoff,
                        SAINTVENANTKIRCHOFFPARAMETERS
                    );
                }
                mod yeoh {
                    use super::*;
                    test_finite_element_block_with_hyperelastic_constitutive_model!(
                        ElasticBlock,
                        $element,
                        Yeoh,
                        YEOHPARAMETERS
                    );
                }
            }
            mod elastic_hyperviscous {
                use super::*;
                use crate::constitutive::solid::elastic_hyperviscous::{
                    test::ALMANSIHAMELPARAMETERS, AlmansiHamel,
                };
                mod almansi_hamel {
                    use super::*;
                    test_finite_element_block_with_elastic_hyperviscous_constitutive_model!(
                        ViscoelasticBlock,
                        $element,
                        AlmansiHamel,
                        ALMANSIHAMELPARAMETERS
                    );
                }
            }
            mod hyperviscoelastic {
                use super::*;
                use crate::constitutive::solid::hyperviscoelastic::{
                    test::SAINTVENANTKIRCHOFFPARAMETERS, SaintVenantKirchoff,
                };
                mod saint_venant_kirchoff {
                    use super::*;
                    test_finite_element_block_with_hyperviscoelastic_constitutive_model!(
                        ViscoelasticBlock,
                        $element,
                        SaintVenantKirchoff,
                        SAINTVENANTKIRCHOFFPARAMETERS
                    );
                }
            }
        }
    };
}
pub(crate) use test_finite_element_block;

macro_rules! setup_for_test_finite_element_block_with_elastic_constitutive_model {
    ($block: ident, $element: ident, $constitutive_model: ident, $constitutive_model_parameters: ident) => {
        fn get_block<'a>() -> $block<D, E, $element<$constitutive_model<'a>>, G, N> {
            $block::<D, E, $element<$constitutive_model<'a>>, G, N>::new(
                $constitutive_model_parameters,
                get_connectivity(),
                get_reference_coordinates_block(),
            )
        }
        fn get_block_transformed<'a>() -> $block<D, E, $element<$constitutive_model<'a>>, G, N> {
            $block::<D, E, $element<$constitutive_model<'a>>, G, N>::new(
                $constitutive_model_parameters,
                get_connectivity(),
                get_reference_coordinates_transformed_block(),
            )
        }
    };
}
pub(crate) use setup_for_test_finite_element_block_with_elastic_constitutive_model;

macro_rules! test_nodal_forces_and_nodal_stiffnesses {
    ($block: ident, $element: ident, $constitutive_model: ident, $constitutive_model_parameters: ident) => {
        setup_for_test_finite_element_block_with_elastic_constitutive_model!(
            $block,
            $element,
            $constitutive_model,
            $constitutive_model_parameters
        );
        fn get_coordinates_transformed_block() -> NodalCoordinates<D> {
            get_coordinates_block()
                .iter()
                .map(|coordinate| {
                    (get_rotation_current_configuration() * coordinate)
                        + get_translation_current_configuration()
                })
                .collect()
        }
        fn get_reference_coordinates_transformed_block() -> ReferenceNodalCoordinates<D> {
            get_reference_coordinates_block()
                .iter()
                .map(|reference_coordinate| {
                    (get_rotation_reference_configuration() * reference_coordinate)
                        + get_translation_reference_configuration()
                })
                .collect()
        }
        mod nodal_forces {
            use super::*;
            mod deformed {
                use super::*;
                #[test]
                fn finite_difference() -> Result<(), TestError> {
                    assert_eq_from_fd(
                        &get_nodal_stiffnesses(true, false)?,
                        &get_finite_difference_of_nodal_forces(true)?,
                    )
                }
                #[test]
                fn objectivity() -> Result<(), TestError> {
                    assert_eq_within_tols(
                        &get_nodal_forces(true, false)?,
                        &get_nodal_forces(true, true)?,
                    )
                }
            }
            mod undeformed {
                use super::*;
                #[test]
                fn finite_difference() -> Result<(), TestError> {
                    assert_eq_from_fd(
                        &get_nodal_stiffnesses(false, false)?,
                        &get_finite_difference_of_nodal_forces(false)?,
                    )
                }
                #[test]
                fn objectivity() -> Result<(), TestError> {
                    assert_eq_within_tols(&get_nodal_forces(false, true)?, &NodalForces::zero())
                }
                #[test]
                fn zero() -> Result<(), TestError> {
                    assert_eq_within_tols(&get_nodal_forces(false, false)?, &NodalForces::zero())
                }
            }
        }
        mod nodal_stiffnesses {
            use super::*;
            mod deformed {
                use super::*;
                #[test]
                fn objectivity() -> Result<(), TestError> {
                    assert_eq_within_tols(
                        &get_nodal_stiffnesses(true, false)?,
                        &get_nodal_stiffnesses(true, true)?,
                    )
                }
            }
            mod undeformed {
                use super::*;
                #[test]
                fn objectivity() -> Result<(), TestError> {
                    assert_eq_within_tols(
                        &get_nodal_stiffnesses(false, false)?,
                        &get_nodal_stiffnesses(false, true)?,
                    )
                }
            }
        }
        #[test]
        fn size() {
            assert_eq!(
                std::mem::size_of::<ElasticBlock<D, E, $element::<$constitutive_model>, G, N>>(),
                std::mem::size_of::<Connectivity<E, N>>()
                    + std::mem::size_of::<NodalCoordinates<D>>()
                    + E * std::mem::size_of::<$element::<$constitutive_model>>()
            )
        }
    };
}
pub(crate) use test_nodal_forces_and_nodal_stiffnesses;

macro_rules! test_helmholtz_free_energy {
    ($block: ident, $element: ident, $constitutive_model: ident, $constitutive_model_parameters: ident) => {
        fn get_finite_difference_of_helmholtz_free_energy(
            is_deformed: bool,
        ) -> Result<NodalForces<D>, TestError> {
            let mut block = get_block();
            let mut finite_difference = 0.0;
            (0..D)
                .map(|node| {
                    (0..3)
                        .map(|i| {
                            let mut nodal_coordinates = if is_deformed {
                                get_coordinates_block()
                            } else {
                                get_reference_coordinates_block().into()
                            };
                            nodal_coordinates[node][i] += 0.5 * EPSILON;
                            block.set_nodal_coordinates(nodal_coordinates);
                            finite_difference = block.calculate_helmholtz_free_energy()?;
                            let mut nodal_coordinates = if is_deformed {
                                get_coordinates_block()
                            } else {
                                get_reference_coordinates_block().into()
                            };
                            nodal_coordinates[node][i] -= 0.5 * EPSILON;
                            block.set_nodal_coordinates(nodal_coordinates);
                            finite_difference -= block.calculate_helmholtz_free_energy()?;
                            Ok(finite_difference / EPSILON)
                        })
                        .collect()
                })
                .collect()
        }
        mod helmholtz_free_energy {
            use super::*;
            mod deformed {
                use super::*;
                #[test]
                fn finite_difference() -> Result<(), TestError> {
                    let mut block = get_block();
                    block.set_nodal_coordinates(get_coordinates_block());
                    assert_eq_from_fd(
                        &block.calculate_nodal_forces()?,
                        &get_finite_difference_of_helmholtz_free_energy(true)?,
                    )
                }
                #[test]
                #[should_panic(expected = "Invalid Jacobian")]
                fn invalid_jacobian() {
                    let mut block = get_block();
                    let mut deformation_gradient = DeformationGradient::identity();
                    deformation_gradient[0][0] = 0.0;
                    let coordinates_block = get_reference_coordinates_block()
                        .iter()
                        .map(|reference_coordinates| &deformation_gradient * reference_coordinates)
                        .collect();
                    block.set_nodal_coordinates(coordinates_block);
                    block.calculate_helmholtz_free_energy().unwrap();
                }
                #[test]
                fn minimized() -> Result<(), TestError> {
                    let mut block = get_block();
                    block.set_nodal_coordinates(get_coordinates_block());
                    let nodal_forces = block.calculate_nodal_forces()?;
                    let minimum = block.calculate_helmholtz_free_energy()?
                        - nodal_forces.dot(&get_coordinates_block());
                    let mut perturbed = 0.0;
                    let mut perturbed_coordinates = get_coordinates_block();
                    (0..D).try_for_each(|node| {
                        (0..3).try_for_each(|i| {
                            perturbed_coordinates = get_coordinates_block();
                            perturbed_coordinates[node][i] += 0.5 * EPSILON;
                            block.set_nodal_coordinates(&perturbed_coordinates * 1.0);
                            perturbed = block.calculate_helmholtz_free_energy()?
                                - nodal_forces.dot(&perturbed_coordinates);
                            if assert_eq_within_tols(&perturbed, &minimum).is_err() {
                                assert!(perturbed > minimum)
                            }
                            perturbed_coordinates[node][i] -= EPSILON;
                            block.set_nodal_coordinates(&perturbed_coordinates * 1.0);
                            perturbed = block.calculate_helmholtz_free_energy()?
                                - nodal_forces.dot(&perturbed_coordinates);
                            if assert_eq_within_tols(&perturbed, &minimum).is_err() {
                                assert!(perturbed > minimum)
                            }
                            Ok(())
                        })
                    })
                }
                #[test]
                fn objectivity() -> Result<(), TestError> {
                    let mut block_1 = get_block();
                    let mut block_2 = get_block_transformed();
                    block_1.set_nodal_coordinates(get_coordinates_block());
                    block_2.set_nodal_coordinates(get_coordinates_transformed_block());
                    assert_eq_within_tols(
                        &block_1.calculate_helmholtz_free_energy()?,
                        &block_2.calculate_helmholtz_free_energy()?,
                    )
                }
                #[test]
                fn positive() -> Result<(), TestError> {
                    let mut block = get_block();
                    assert_eq_within_tols(&block.calculate_helmholtz_free_energy()?.abs(), &0.0)?;
                    block.set_nodal_coordinates(get_coordinates_block());
                    assert!(block.calculate_helmholtz_free_energy()? > 0.0);
                    Ok(())
                }
            }
            mod undeformed {
                use super::*;
                #[test]
                fn finite_difference() -> Result<(), TestError> {
                    assert_eq_from_fd(
                        &get_finite_difference_of_helmholtz_free_energy(false)?,
                        &NodalForces::zero(),
                    )
                }
                #[test]
                fn minimized() -> Result<(), TestError> {
                    let mut block = get_block();
                    let minimum = block.calculate_helmholtz_free_energy()?;
                    let mut perturbed = 0.0;
                    let mut perturbed_coordinates = get_reference_coordinates_block();
                    (0..D).try_for_each(|node| {
                        (0..3).try_for_each(|i| {
                            perturbed_coordinates = get_reference_coordinates_block();
                            perturbed_coordinates[node][i] += 0.5 * EPSILON;
                            block.set_nodal_coordinates(perturbed_coordinates.convert());
                            perturbed = block.calculate_helmholtz_free_energy()?;
                            if assert_eq_within_tols(&perturbed, &minimum).is_err() {
                                assert!(perturbed > minimum)
                            }
                            perturbed_coordinates[node][i] -= EPSILON;
                            block.set_nodal_coordinates(perturbed_coordinates.convert());
                            perturbed = block.calculate_helmholtz_free_energy()?;
                            if assert_eq_within_tols(&perturbed, &minimum).is_err() {
                                assert!(perturbed > minimum)
                            }
                            Ok(())
                        })
                    })
                }
                #[test]
                fn objectivity() -> Result<(), TestError> {
                    let mut block_1 = get_block();
                    let mut block_2 = get_block_transformed();
                    assert_eq_within_tols(
                        &block_1.calculate_helmholtz_free_energy()?,
                        &block_2.calculate_helmholtz_free_energy()?,
                    )?;
                    block_1.set_nodal_coordinates(get_reference_coordinates_block().into());
                    block_2.set_nodal_coordinates(
                        get_reference_coordinates_transformed_block().into(),
                    );
                    assert_eq_within_tols(
                        &block_1.calculate_helmholtz_free_energy()?,
                        &block_2.calculate_helmholtz_free_energy()?,
                    )
                }
                #[test]
                fn zero() -> Result<(), TestError> {
                    let mut block = get_block();
                    assert_eq_within_tols(&block.calculate_helmholtz_free_energy()?.abs(), &0.0)?;
                    block.set_nodal_coordinates(get_reference_coordinates_block().into());
                    assert_eq_within_tols(&block.calculate_helmholtz_free_energy()?.abs(), &0.0)
                }
            }
        }
        #[test]
        fn nodal_stiffnesses_deformed_symmetry() -> Result<(), TestError> {
            let nodal_stiffness = get_nodal_stiffnesses(true, false)?;
            let result =
                nodal_stiffness
                    .iter()
                    .enumerate()
                    .try_for_each(|(a, nodal_stiffness_a)| {
                        nodal_stiffness_a.iter().enumerate().try_for_each(
                            |(b, nodal_stiffness_ab)| {
                                nodal_stiffness_ab
                                    .iter()
                                    .enumerate()
                                    .zip(nodal_stiffness_ab.transpose().iter())
                                    .try_for_each(
                                        |((i, nodal_stiffness_ab_i), nodal_stiffness_ab_j)| {
                                            nodal_stiffness_ab_i
                                                .iter()
                                                .enumerate()
                                                .zip(nodal_stiffness_ab_j.iter())
                                                .try_for_each(
                                                    |(
                                                        (j, nodal_stiffness_ab_ij),
                                                        nodal_stiffness_ab_ji,
                                                    )| {
                                                        if a == b {
                                                            assert_eq_within_tols(
                                                                nodal_stiffness_ab_ij,
                                                                &nodal_stiffness_ab_ji,
                                                            )
                                                        } else if i == j {
                                                            assert_eq_within_tols(
                                                                nodal_stiffness_ab_ij,
                                                                &nodal_stiffness[b][a][i][j],
                                                            )
                                                        } else {
                                                            Ok(())
                                                        }
                                                    },
                                                )
                                        },
                                    )
                            },
                        )
                    });
            result
        }
        #[test]
        fn nodal_stiffnesses_undeformed_symmetry() -> Result<(), TestError> {
            let nodal_stiffness = get_nodal_stiffnesses(false, false)?;
            let result =
                nodal_stiffness
                    .iter()
                    .enumerate()
                    .try_for_each(|(a, nodal_stiffness_a)| {
                        nodal_stiffness_a.iter().enumerate().try_for_each(
                            |(b, nodal_stiffness_ab)| {
                                nodal_stiffness_ab
                                    .iter()
                                    .enumerate()
                                    .zip(nodal_stiffness_ab.transpose().iter())
                                    .try_for_each(
                                        |((i, nodal_stiffness_ab_i), nodal_stiffness_ab_j)| {
                                            nodal_stiffness_ab_i
                                                .iter()
                                                .enumerate()
                                                .zip(nodal_stiffness_ab_j.iter())
                                                .try_for_each(
                                                    |(
                                                        (j, nodal_stiffness_ab_ij),
                                                        nodal_stiffness_ab_ji,
                                                    )| {
                                                        if a == b {
                                                            assert_eq_within_tols(
                                                                nodal_stiffness_ab_ij,
                                                                &nodal_stiffness_ab_ji,
                                                            )
                                                        } else if i == j {
                                                            assert_eq_within_tols(
                                                                nodal_stiffness_ab_ij,
                                                                &nodal_stiffness[b][a][i][j],
                                                            )
                                                        } else {
                                                            Ok(())
                                                        }
                                                    },
                                                )
                                        },
                                    )
                            },
                        )
                    });
            result
        }
    };
}
pub(crate) use test_helmholtz_free_energy;

macro_rules! test_finite_element_block_with_elastic_constitutive_model {
    ($block: ident, $element: ident, $constitutive_model: ident, $constitutive_model_parameters: ident) => {
        fn get_finite_difference_of_nodal_forces(
            is_deformed: bool,
        ) -> Result<NodalStiffnesses<D>, TestError> {
            let mut block = get_block();
            let mut finite_difference = 0.0;
            (0..D)
                .map(|node_a| {
                    (0..D)
                        .map(|node_b| {
                            (0..3)
                                .map(|i| {
                                    (0..3)
                                        .map(|j| {
                                            let mut nodal_coordinates = if is_deformed {
                                                get_coordinates_block()
                                            } else {
                                                get_reference_coordinates_block().into()
                                            };
                                            nodal_coordinates[node_b][j] += 0.5 * EPSILON;
                                            block.set_nodal_coordinates(nodal_coordinates);
                                            finite_difference =
                                                block.calculate_nodal_forces()?[node_a][i];
                                            let mut nodal_coordinates = if is_deformed {
                                                get_coordinates_block()
                                            } else {
                                                get_reference_coordinates_block().into()
                                            };
                                            nodal_coordinates[node_b][j] -= 0.5 * EPSILON;
                                            block.set_nodal_coordinates(nodal_coordinates);
                                            finite_difference -=
                                                block.calculate_nodal_forces()?[node_a][i];
                                            Ok(finite_difference / EPSILON)
                                        })
                                        .collect()
                                })
                                .collect()
                        })
                        .collect()
                })
                .collect()
        }
        fn get_nodal_forces(
            is_deformed: bool,
            is_rotated: bool,
        ) -> Result<NodalForces<D>, TestError> {
            if is_rotated {
                if is_deformed {
                    let mut block = get_block_transformed();
                    block.set_nodal_coordinates(get_coordinates_transformed_block());
                    Ok(get_rotation_current_configuration().transpose()
                        * block.calculate_nodal_forces()?)
                } else {
                    let converted: TensorRank2<3, 1, 1> =
                        get_rotation_reference_configuration().into();
                    Ok(converted.transpose() * get_block_transformed().calculate_nodal_forces()?)
                }
            } else {
                if is_deformed {
                    let mut block = get_block();
                    block.set_nodal_coordinates(get_coordinates_block());
                    Ok(block.calculate_nodal_forces()?)
                } else {
                    Ok(get_block().calculate_nodal_forces()?)
                }
            }
        }
        fn get_nodal_stiffnesses(
            is_deformed: bool,
            is_rotated: bool,
        ) -> Result<NodalStiffnesses<D>, TestError> {
            if is_rotated {
                if is_deformed {
                    let mut block = get_block_transformed();
                    block.set_nodal_coordinates(get_coordinates_transformed_block());
                    Ok(get_rotation_current_configuration().transpose()
                        * block.calculate_nodal_stiffnesses()?
                        * get_rotation_current_configuration())
                } else {
                    let converted: TensorRank2<3, 1, 1> =
                        get_rotation_reference_configuration().into();
                    Ok(converted.transpose()
                        * get_block_transformed().calculate_nodal_stiffnesses()?
                        * converted)
                }
            } else {
                if is_deformed {
                    let mut block = get_block();
                    block.set_nodal_coordinates(get_coordinates_block());
                    Ok(block.calculate_nodal_stiffnesses()?)
                } else {
                    Ok(get_block().calculate_nodal_stiffnesses()?)
                }
            }
        }
        crate::fem::block::test::test_nodal_forces_and_nodal_stiffnesses!(
            $block,
            $element,
            $constitutive_model,
            $constitutive_model_parameters
        );
    };
}
pub(crate) use test_finite_element_block_with_elastic_constitutive_model;

macro_rules! test_finite_element_block_with_hyperelastic_constitutive_model {
    ($block: ident, $element: ident, $constitutive_model: ident, $constitutive_model_parameters: ident) => {
        crate::fem::block::test::test_finite_element_block_with_elastic_constitutive_model!(
            $block,
            $element,
            $constitutive_model,
            $constitutive_model_parameters
        );
        crate::fem::block::test::test_helmholtz_free_energy!(
            $block,
            $element,
            $constitutive_model,
            $constitutive_model_parameters
        );
    };
}
pub(crate) use test_finite_element_block_with_hyperelastic_constitutive_model;

macro_rules! test_finite_element_block_with_viscoelastic_constitutive_model {
    ($block: ident, $element: ident, $constitutive_model: ident, $constitutive_model_parameters: ident) => {
        fn get_velocities_transformed_block() -> NodalCoordinates<D> {
            get_coordinates_block()
                .iter()
                .zip(get_velocities_block().iter())
                .map(|(coordinate, velocity)| {
                    get_rotation_current_configuration() * velocity
                        + get_rotation_rate_current_configuration() * coordinate
                        + get_translation_rate_current_configuration()
                })
                .collect()
        }
        fn get_nodal_forces(
            is_deformed: bool,
            is_rotated: bool,
        ) -> Result<NodalForces<D>, TestError> {
            if is_rotated {
                if is_deformed {
                    let mut block = get_block_transformed();
                    block.set_nodal_coordinates(get_coordinates_transformed_block());
                    block.set_nodal_velocities(get_velocities_transformed_block());
                    Ok(get_rotation_current_configuration().transpose()
                        * block.calculate_nodal_forces()?)
                } else {
                    let converted: TensorRank2<3, 1, 1> =
                        get_rotation_reference_configuration().into();
                    Ok(converted.transpose() * get_block_transformed().calculate_nodal_forces()?)
                }
            } else {
                if is_deformed {
                    let mut block = get_block();
                    block.set_nodal_coordinates(get_coordinates_block());
                    block.set_nodal_velocities(get_velocities_block());
                    Ok(block.calculate_nodal_forces()?)
                } else {
                    Ok(get_block().calculate_nodal_forces()?)
                }
            }
        }
        fn get_nodal_stiffnesses(
            is_deformed: bool,
            is_rotated: bool,
        ) -> Result<NodalStiffnesses<D>, TestError> {
            if is_rotated {
                if is_deformed {
                    let mut block = get_block_transformed();
                    block.set_nodal_coordinates(get_coordinates_transformed_block());
                    block.set_nodal_velocities(get_velocities_transformed_block());
                    Ok(get_rotation_current_configuration().transpose()
                        * block.calculate_nodal_stiffnesses()?
                        * get_rotation_current_configuration())
                } else {
                    let converted: TensorRank2<3, 1, 1> =
                        get_rotation_reference_configuration().into();
                    Ok(converted.transpose()
                        * get_block_transformed().calculate_nodal_stiffnesses()?
                        * converted)
                }
            } else {
                if is_deformed {
                    let mut block = get_block();
                    block.set_nodal_coordinates(get_coordinates_block());
                    block.set_nodal_velocities(get_velocities_block());
                    Ok(block.calculate_nodal_stiffnesses()?)
                } else {
                    Ok(get_block().calculate_nodal_stiffnesses()?)
                }
            }
        }
        fn get_finite_difference_of_nodal_forces(
            is_deformed: bool,
        ) -> Result<NodalStiffnesses<D>, TestError> {
            let mut block = get_block();
            if is_deformed {
                block.set_nodal_coordinates(get_coordinates_block());
            } else {
                block.set_nodal_coordinates(get_reference_coordinates_block().into());
            }
            let mut finite_difference = 0.0;
            (0..D)
                .map(|node_a| {
                    (0..D)
                        .map(|node_b| {
                            (0..3)
                                .map(|i| {
                                    (0..3)
                                        .map(|j| {
                                            let mut nodal_velocities = if is_deformed {
                                                get_velocities_block()
                                            } else {
                                                NodalVelocities::zero()
                                            };
                                            nodal_velocities[node_a][i] += 0.5 * EPSILON;
                                            block.set_nodal_velocities(nodal_velocities);
                                            finite_difference =
                                                block.calculate_nodal_forces()?[node_b][j];
                                            let mut nodal_velocities = if is_deformed {
                                                get_velocities_block()
                                            } else {
                                                NodalVelocities::zero()
                                            };
                                            nodal_velocities[node_a][i] -= 0.5 * EPSILON;
                                            block.set_nodal_velocities(nodal_velocities);
                                            finite_difference -=
                                                block.calculate_nodal_forces()?[node_b][j];
                                            Ok(finite_difference / EPSILON)
                                        })
                                        .collect()
                                })
                                .collect()
                        })
                        .collect()
                })
                .collect()
        }
        crate::fem::block::test::test_nodal_forces_and_nodal_stiffnesses!(
            $block,
            $element,
            $constitutive_model,
            $constitutive_model_parameters
        );
    };
}
pub(crate) use test_finite_element_block_with_viscoelastic_constitutive_model;

macro_rules! test_finite_element_block_with_elastic_hyperviscous_constitutive_model {
    ($block: ident, $element: ident, $constitutive_model: ident, $constitutive_model_parameters: ident) => {
        crate::fem::block::test::test_finite_element_block_with_viscoelastic_constitutive_model!(
            $block,
            $element,
            $constitutive_model,
            $constitutive_model_parameters
        );
        fn get_finite_difference_of_viscous_dissipation(
            is_deformed: bool,
        ) -> Result<NodalForces<D>, TestError> {
            let mut block = get_block();
            if is_deformed {
                block.set_nodal_coordinates(get_coordinates_block());
            } else {
                block.set_nodal_coordinates(get_reference_coordinates_block().into());
            }
            let mut finite_difference = 0.0;
            (0..D)
                .map(|node| {
                    (0..3)
                        .map(|i| {
                            let mut nodal_velocities = if is_deformed {
                                get_velocities_block()
                            } else {
                                NodalVelocities::zero()
                            };
                            nodal_velocities[node][i] += 0.5 * EPSILON;
                            block.set_nodal_velocities(nodal_velocities);
                            finite_difference = block.calculate_viscous_dissipation()?;
                            let mut nodal_velocities = if is_deformed {
                                get_velocities_block()
                            } else {
                                NodalVelocities::zero()
                            };
                            nodal_velocities[node][i] -= 0.5 * EPSILON;
                            block.set_nodal_velocities(nodal_velocities);
                            finite_difference -= block.calculate_viscous_dissipation()?;
                            Ok(finite_difference / EPSILON)
                        })
                        .collect()
                })
                .collect()
        }
        fn get_finite_difference_of_dissipation_potential(
            is_deformed: bool,
        ) -> Result<NodalForces<D>, TestError> {
            let mut block = get_block();
            if is_deformed {
                block.set_nodal_coordinates(get_coordinates_block());
            } else {
                block.set_nodal_coordinates(get_reference_coordinates_block().into());
            }
            let mut finite_difference = 0.0;
            (0..D)
                .map(|node| {
                    (0..3)
                        .map(|i| {
                            let mut nodal_velocities = if is_deformed {
                                get_velocities_block()
                            } else {
                                NodalVelocities::zero()
                            };
                            nodal_velocities[node][i] += 0.5 * EPSILON;
                            block.set_nodal_velocities(nodal_velocities);
                            finite_difference = block.calculate_dissipation_potential()?;
                            let mut nodal_velocities = if is_deformed {
                                get_velocities_block()
                            } else {
                                NodalVelocities::zero()
                            };
                            nodal_velocities[node][i] -= 0.5 * EPSILON;
                            block.set_nodal_velocities(nodal_velocities);
                            finite_difference -= block.calculate_dissipation_potential()?;
                            Ok(finite_difference / EPSILON)
                        })
                        .collect()
                })
                .collect()
        }
        mod viscous_dissipation {
            use super::*;
            mod deformed {
                use super::*;
                #[test]
                fn finite_difference() -> Result<(), TestError> {
                    let mut block = get_block();
                    block.set_nodal_coordinates(get_coordinates_block());
                    let nodal_forces_0 = block.calculate_nodal_forces()?;
                    block.set_nodal_velocities(get_velocities_block());
                    assert_eq_from_fd(
                        &(block.calculate_nodal_forces()? - nodal_forces_0),
                        &get_finite_difference_of_viscous_dissipation(true)?,
                    )
                }
                #[test]
                fn minimized() -> Result<(), TestError> {
                    let mut block = get_block();
                    block.set_nodal_coordinates(get_coordinates_block());
                    let nodal_forces_0 = block.calculate_nodal_forces()?;
                    block.set_nodal_velocities(get_velocities_block());
                    let nodal_forces = block.calculate_nodal_forces()? - nodal_forces_0;
                    let minimum = block.calculate_viscous_dissipation()?
                        - nodal_forces.dot(&get_velocities_block());
                    let mut perturbed_velocities = get_velocities_block();
                    (0..D).try_for_each(|node| {
                        (0..3).try_for_each(|i| {
                            perturbed_velocities = get_velocities_block();
                            perturbed_velocities[node][i] += 0.5 * EPSILON;
                            block.set_nodal_velocities(&perturbed_velocities * 1.0);
                            assert!(
                                block.calculate_viscous_dissipation()?
                                    - nodal_forces.dot(&perturbed_velocities)
                                    >= minimum
                            );
                            perturbed_velocities[node][i] -= EPSILON;
                            block.set_nodal_velocities(&perturbed_velocities * 1.0);
                            assert!(
                                block.calculate_viscous_dissipation()?
                                    - nodal_forces.dot(&perturbed_velocities)
                                    >= minimum
                            );
                            Ok(())
                        })
                    })
                }
                #[test]
                fn objectivity() -> Result<(), TestError> {
                    let mut block_1 = get_block();
                    let mut block_2 = get_block_transformed();
                    block_1.set_nodal_coordinates(get_coordinates_block());
                    block_1.set_nodal_velocities(get_velocities_block());
                    block_2.set_nodal_coordinates(get_coordinates_transformed_block());
                    block_2.set_nodal_velocities(get_velocities_transformed_block());
                    assert_eq_within_tols(
                        &block_1.calculate_viscous_dissipation()?,
                        &block_2.calculate_viscous_dissipation()?,
                    )
                }
                #[test]
                fn positive() -> Result<(), TestError> {
                    let mut block = get_block();
                    assert_eq(&block.calculate_viscous_dissipation()?, &0.0)?;
                    block.set_nodal_coordinates(get_coordinates_block());
                    block.set_nodal_velocities(get_velocities_block());
                    assert!(block.calculate_viscous_dissipation()? > 0.0);
                    Ok(())
                }
            }
            mod undeformed {
                use super::*;
                #[test]
                fn finite_difference() -> Result<(), TestError> {
                    assert_eq_from_fd(
                        &get_finite_difference_of_viscous_dissipation(false)?,
                        &NodalForces::zero(),
                    )
                }
                #[test]
                fn minimized() -> Result<(), TestError> {
                    let mut block = get_block();
                    let minimum = block.calculate_viscous_dissipation()?;
                    let mut perturbed_velocities = NodalVelocities::zero();
                    (0..D).try_for_each(|node| {
                        (0..3).try_for_each(|i| {
                            perturbed_velocities = NodalVelocities::zero();
                            perturbed_velocities[node][i] += 0.5 * EPSILON;
                            block.set_nodal_velocities(perturbed_velocities.convert());
                            assert!(block.calculate_viscous_dissipation()? >= minimum);
                            perturbed_velocities[node][i] -= EPSILON;
                            block.set_nodal_velocities(perturbed_velocities.convert());
                            assert!(block.calculate_viscous_dissipation()? >= minimum);
                            Ok(())
                        })
                    })
                }
                #[test]
                fn objectivity() -> Result<(), TestError> {
                    let mut block_1 = get_block();
                    let mut block_2 = get_block_transformed();
                    assert_eq_within_tols(
                        &block_1.calculate_viscous_dissipation()?,
                        &block_2.calculate_viscous_dissipation()?,
                    )?;
                    block_1.set_nodal_coordinates(get_reference_coordinates_block().into());
                    block_1.set_nodal_velocities(NodalVelocities::zero());
                    block_2.set_nodal_coordinates(
                        get_reference_coordinates_transformed_block().into(),
                    );
                    block_2.set_nodal_velocities(NodalVelocities::zero());
                    assert_eq_within_tols(
                        &block_1.calculate_viscous_dissipation()?,
                        &block_2.calculate_viscous_dissipation()?,
                    )
                }
                #[test]
                fn zero() -> Result<(), TestError> {
                    let mut block = get_block();
                    assert_eq(&block.calculate_viscous_dissipation()?, &0.0)?;
                    block.set_nodal_coordinates(get_reference_coordinates_block().into());
                    assert_eq(&block.calculate_viscous_dissipation()?, &0.0)
                }
            }
        }
        mod dissipation_potential {
            use super::*;
            mod deformed {
                use super::*;
                #[test]
                fn finite_difference() -> Result<(), TestError> {
                    let mut block = get_block();
                    block.set_nodal_coordinates(get_coordinates_block());
                    block.set_nodal_velocities(get_velocities_block());
                    assert_eq_from_fd(
                        &block.calculate_nodal_forces()?,
                        &get_finite_difference_of_dissipation_potential(true)?,
                    )
                }
                #[test]
                fn minimized() -> Result<(), TestError> {
                    let mut block = get_block();
                    block.set_nodal_coordinates(get_coordinates_block());
                    block.set_nodal_velocities(get_velocities_block());
                    let nodal_forces = block.calculate_nodal_forces()?;
                    let minimum = block.calculate_dissipation_potential()?
                        - nodal_forces.dot(&get_velocities_block());
                    let mut perturbed_velocities = get_velocities_block();
                    (0..D).try_for_each(|node| {
                        (0..3).try_for_each(|i| {
                            perturbed_velocities = get_velocities_block();
                            perturbed_velocities[node][i] += 0.5 * EPSILON;
                            block.set_nodal_velocities(&perturbed_velocities * 1.0);
                            assert!(
                                block.calculate_dissipation_potential()?
                                    - nodal_forces.dot(&perturbed_velocities)
                                    >= minimum
                            );
                            perturbed_velocities[node][i] -= EPSILON;
                            block.set_nodal_velocities(&perturbed_velocities * 1.0);
                            assert!(
                                block.calculate_dissipation_potential()?
                                    - nodal_forces.dot(&perturbed_velocities)
                                    >= minimum
                            );
                            Ok(())
                        })
                    })
                }
                #[test]
                fn objectivity() -> Result<(), TestError> {
                    let mut block_1 = get_block();
                    let mut block_2 = get_block_transformed();
                    block_1.set_nodal_coordinates(get_coordinates_block());
                    block_1.set_nodal_velocities(get_velocities_block());
                    block_2.set_nodal_coordinates(get_coordinates_transformed_block());
                    block_2.set_nodal_velocities(get_velocities_transformed_block());
                    assert_eq_within_tols(
                        &block_1.calculate_dissipation_potential()?,
                        &block_2.calculate_dissipation_potential()?,
                    )
                }
            }
            mod undeformed {
                use super::*;
                #[test]
                fn finite_difference() -> Result<(), TestError> {
                    assert_eq_from_fd(
                        &get_finite_difference_of_dissipation_potential(false)?,
                        &NodalForces::zero(),
                    )
                }
                #[test]
                fn minimized() -> Result<(), TestError> {
                    let mut block = get_block();
                    let minimum = block.calculate_dissipation_potential()?;
                    let mut perturbed_velocities = NodalVelocities::zero();
                    (0..D).try_for_each(|node| {
                        (0..3).try_for_each(|i| {
                            perturbed_velocities = NodalVelocities::zero();
                            perturbed_velocities[node][i] += 0.5 * EPSILON;
                            block.set_nodal_velocities(perturbed_velocities.convert());
                            assert!(block.calculate_dissipation_potential()? >= minimum);
                            perturbed_velocities[node][i] -= EPSILON;
                            block.set_nodal_velocities(perturbed_velocities.convert());
                            assert!(block.calculate_dissipation_potential()? >= minimum);
                            Ok(())
                        })
                    })
                }
                #[test]
                fn objectivity() -> Result<(), TestError> {
                    let mut block_1 = get_block();
                    let mut block_2 = get_block_transformed();
                    assert_eq_within_tols(
                        &block_1.calculate_dissipation_potential()?,
                        &block_2.calculate_dissipation_potential()?,
                    )?;
                    block_1.set_nodal_coordinates(get_reference_coordinates_block().into());
                    block_1.set_nodal_velocities(NodalVelocities::zero());
                    block_2.set_nodal_coordinates(
                        get_reference_coordinates_transformed_block().into(),
                    );
                    block_2.set_nodal_velocities(NodalVelocities::zero());
                    assert_eq_within_tols(
                        &block_1.calculate_dissipation_potential()?,
                        &block_2.calculate_dissipation_potential()?,
                    )
                }
                #[test]
                fn zero() -> Result<(), TestError> {
                    let mut block = get_block();
                    assert_eq(&block.calculate_dissipation_potential()?, &0.0)?;
                    block.set_nodal_coordinates(get_reference_coordinates_block().into());
                    assert_eq(&block.calculate_dissipation_potential()?, &0.0)
                }
            }
        }
    };
}
pub(crate) use test_finite_element_block_with_elastic_hyperviscous_constitutive_model;

macro_rules! test_finite_element_block_with_hyperviscoelastic_constitutive_model
{
    ($block: ident, $element: ident, $constitutive_model: ident, $constitutive_model_parameters: ident) =>
    {
        crate::fem::block::test::test_finite_element_block_with_elastic_hyperviscous_constitutive_model!(
            $block, $element, $constitutive_model, $constitutive_model_parameters
        );
        crate::fem::block::test::test_helmholtz_free_energy!(
            $block, $element, $constitutive_model, $constitutive_model_parameters
        );
        #[test]
        fn dissipation_potential_deformed_positive() -> Result<(), TestError>
        {
            let mut block = get_block();
            block.set_nodal_coordinates(
                get_coordinates_block()
            );
            block.set_nodal_velocities(
                get_velocities_block()
            );
            assert!(block.calculate_dissipation_potential()? > 0.0);
            Ok(())
        }
    }
}
pub(crate) use test_finite_element_block_with_hyperviscoelastic_constitutive_model;
