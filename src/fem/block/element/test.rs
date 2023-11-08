macro_rules! test_finite_element
{
    ($element: ident) =>
    {
        use crate::
        {
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
                    todo!()
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
                    todo!()
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