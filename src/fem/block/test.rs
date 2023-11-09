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
                        GentModel,
                        MooneyRivlinModel,
                        NeoHookeanModel,
                        YeohModel,
                    },
                    test::
                    {
                        GENTPARAMETERS,
                        MOONEYRIVLINPARAMETERS,
                        NEOHOOKEANPARAMETERS,
                        YEOHPARAMETERS
                    }
                },
                fem::block::test::test_finite_element_block_with_constitutive_model,
                math::Convert,
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
            get_current_coordinates_block().iter().map(|current_coordinate|
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
        fn get_reference_coordinates_transformed_block() -> ReferenceNodalCoordinates<D>
        {
            get_reference_coordinates_block().iter().map(|reference_coordinate|
                (get_rotation_reference_configuration() * reference_coordinate)
                + get_translation_reference_configuration()
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
                    todo!()
                }
                #[test]
                fn objectivity()
                {
                    todo!()
                }
            }
            mod undeformed
            {
                use super::*;
                #[test]
                fn finite_difference()
                {
                    todo!()
                }
                #[test]
                fn objectivity()
                {
                    todo!()
                }
                #[test]
                fn zero()
                {
                    todo!()
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
                fn symmetry()
                {
                    todo!()
                }
                #[test]
                fn objectivity()
                {
                    todo!()
                }
            }
            mod undeformed
            {
                use super::*;
                #[test]
                fn symmetry()
                {
                    todo!()
                }
                #[test]
                fn objectivity()
                {
                    todo!()
                }
            }
        }
    }
}
pub(crate) use test_finite_element_block_with_constitutive_model;