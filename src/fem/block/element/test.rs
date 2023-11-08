macro_rules! test_finite_element
{
    ($element: ident) =>
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
            fem::block::element::test::test_finite_element_for_constitutive_model,
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
        mod gent
        {
            use super::*;
            test_finite_element_for_constitutive_model!($element, GentModel, GENTPARAMETERS);
        }
        mod mooney_rivlin
        {
            use super::*;
            test_finite_element_for_constitutive_model!($element, MooneyRivlinModel, MOONEYRIVLINPARAMETERS);
        }
        mod neo_hookean
        {
            use super::*;
            test_finite_element_for_constitutive_model!($element, NeoHookeanModel, NEOHOOKEANPARAMETERS);
        }
        mod yeoh
        {
            use super::*;
            test_finite_element_for_constitutive_model!($element, YeohModel, YEOHPARAMETERS);
        }
    }
}
pub(crate) use test_finite_element;
macro_rules! test_finite_element_for_constitutive_model
{
    ($element: ident, $constitutive_model: ident, $constitutive_model_parameters: ident) =>
    {
        fn get_current_coordinates() -> CurrentNodalCoordinates<N>
        {
            get_reference_coordinates().iter()
            .map(|reference_coordinate|
                get_deformation_gradient() * reference_coordinate
            ).collect()
        }
        fn get_current_coordinates_transformed() -> CurrentNodalCoordinates<N>
        {
            get_current_coordinates().iter()
            .map(|current_coordinate|
                get_rotation_current_configuration() * current_coordinate
                + get_translation_current_configuration()
            ).collect()
        }
        fn get_element<'a>() -> $element<'a, $constitutive_model<'a>>
        {
            $element::new(
                $constitutive_model_parameters,
                get_reference_coordinates()
            )
        }
        fn get_fd_helmholtz_free_energy(is_deformed: bool) -> NodalForces<N>
        {
            let element = get_element();
            let mut finite_difference = 0.0;
            (0..N).map(|node|
                (0..3).map(|index|{
                    let mut nodal_coordinates = 
                    if is_deformed
                    {
                        get_current_coordinates()
                    }
                    else
                    {
                        get_reference_coordinates()
                        .convert()
                    };
                    nodal_coordinates[node][index] += 0.5 * EPSILON;
                    finite_difference = element.calculate_helmholtz_free_energy(&nodal_coordinates);
                    nodal_coordinates[node][index] -= EPSILON;
                    finite_difference -= element.calculate_helmholtz_free_energy(&nodal_coordinates);
                    finite_difference/EPSILON
                }).collect()
            ).collect()
        }
        fn get_reference_coordinates_transformed() -> ReferenceNodalCoordinates<N>
        {
            get_reference_coordinates().iter()
            .map(|reference_coordinate|
                get_rotation_reference_configuration() * reference_coordinate
                + get_translation_reference_configuration()
            ).collect()
        }
        fn get_element_transformed<'a>() -> $element<'a, $constitutive_model<'a>>
        {
            $element::<$constitutive_model>::new
            (
                $constitutive_model_parameters,
                get_reference_coordinates_transformed()
            )
        }
        #[test]
        fn integration_weights_sum_to_one()
        {
            assert_eq!(get_element().get_integration_weights().iter().sum::<Scalar>(), 1.0)
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
                        &get_current_coordinates()
                    ).iter()
                    .zip(get_fd_helmholtz_free_energy(
                        true
                    ).iter())
                    .for_each(|(nodal_force_i, fd_nodal_force_i)|
                        nodal_force_i.iter()
                        .zip(fd_nodal_force_i.iter())
                        .for_each(|(nodal_force_ij, fd_nodal_force_ij)|
                            assert!((nodal_force_ij/fd_nodal_force_ij - 1.0).abs() < EPSILON)
                        )
                    )
                }
                #[test]
                fn objectivity()
                {
                    assert_eq_within_tols(
                        &get_element()
                        .calculate_helmholtz_free_energy(
                            &get_current_coordinates()
                        ),
                        &get_element_transformed()
                        .calculate_helmholtz_free_energy(
                            &get_current_coordinates_transformed()
                        )
                    )
                }
                #[test]
                fn positive()
                {
                    assert!(
                        get_element()
                        .calculate_helmholtz_free_energy(
                            &get_current_coordinates()
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
                    get_element().calculate_nodal_forces(
                        &get_reference_coordinates()
                        .convert()
                    ).iter()
                    .zip(get_fd_helmholtz_free_energy(
                        false
                    ).iter())
                    .for_each(|(nodal_force_i, fd_nodal_force_i)|
                        nodal_force_i.iter()
                        .zip(fd_nodal_force_i.iter())
                        .for_each(|(nodal_force_ij, fd_nodal_force_ij)|
                            assert!(fd_nodal_force_ij.abs() < EPSILON)
                        )
                    )
                }
                #[test]
                fn objectivity()
                {
                    assert_eq_within_tols(
                        &get_element()
                        .calculate_helmholtz_free_energy(
                            &get_reference_coordinates()
                            .convert()
                        ),
                        &get_element_transformed()
                        .calculate_helmholtz_free_energy(
                            &get_reference_coordinates_transformed()
                            .convert()
                        )
                    )
                }
                #[test]
                fn zero()
                {
                    assert_eq!(
                        get_element()
                        .calculate_helmholtz_free_energy(
                            &get_reference_coordinates()
                            .convert()
                        ), 0.0
                    )
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
                    get_element().calculate_nodal_forces(
                        &get_current_coordinates()
                    ).iter().zip((
                        get_rotation_current_configuration().transpose() *
                        get_element().calculate_nodal_forces(
                            &get_current_coordinates_transformed()
                        )
                    ).iter()
                    ).for_each(|(nodal_force, res_nodal_force)|
                        nodal_force.iter()
                        .zip(res_nodal_force.iter())
                        .for_each(|(nodal_force_i, res_nodal_force_i)|
                            assert_eq_within_tols(nodal_force_i, res_nodal_force_i)
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
                    todo!()
                }
                #[test]
                fn objectivity()
                {
                    get_element().calculate_nodal_forces(
                        &get_reference_coordinates_transformed()
                        .convert()
                    ).iter()
                    .for_each(|nodal_force|
                        nodal_force.iter()
                        .for_each(|nodal_force_i|
                        assert_eq_within_tols(nodal_force_i, &0.0)
                        )
                    )
                }
                #[test]
                fn zero()
                {
                    get_element().calculate_nodal_forces(
                        &get_reference_coordinates()
                        .convert()
                    ).iter().for_each(|nodal_force|
                        nodal_force.iter().for_each(|nodal_force_i|
                            assert_eq_within_tols(nodal_force_i, &0.0)
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
                    todo!()
                }
                #[test]
                fn symmetry()
                {
                    todo!()
                }
            }
            mod undeformed
            {
                use super::*;
                #[test]
                fn objectivity()
                {
                    todo!()
                }
                #[test]
                fn symmetry()
                {
                    todo!()
                }
            }
        }
    }
}
pub(crate) use test_finite_element_for_constitutive_model;