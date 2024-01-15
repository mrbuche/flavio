macro_rules! test_linear_finite_element
{
    ($element: ident) =>
    {
        mod linear_finite_element
        {
            use crate::
            {
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
                fem::block::element::linear::test::test_linear_finite_element_with_constitutive_model,
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
            pub mod arruda_boyce
            {
                use super::*;
                test_linear_finite_element_with_constitutive_model!($element, ArrudaBoyceModel, ARRUDABOYCEPARAMETERS);
            }
            mod gent
            {
                use super::*;
                test_linear_finite_element_with_constitutive_model!($element, GentModel, GENTPARAMETERS);
            }
            mod mooney_rivlin
            {
                use super::*;
                test_linear_finite_element_with_constitutive_model!($element, MooneyRivlinModel, MOONEYRIVLINPARAMETERS);
            }
            mod neo_hookean
            {
                use super::*;
                test_linear_finite_element_with_constitutive_model!($element, NeoHookeanModel, NEOHOOKEANPARAMETERS);
            }
            mod yeoh
            {
                use super::*;
                test_linear_finite_element_with_constitutive_model!($element, YeohModel, YEOHPARAMETERS);
            }
        }
    }
}
pub(crate) use test_linear_finite_element;
macro_rules! test_linear_finite_element_with_constitutive_model
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
        fn get_element_transformed<'a>() -> $element<'a, $constitutive_model<'a>>
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
        mod deformation_gradient
        {
            use super::*;
            mod deformed
            {
                use super::*;
                #[test]
                fn objectivity()
                {
                    get_element().calculate_deformation_gradient(
                        &get_current_coordinates()
                    ).iter().zip((
                        get_rotation_current_configuration().transpose() *
                        get_element_transformed().calculate_deformation_gradient(
                            &get_current_coordinates_transformed()
                        ) * get_rotation_reference_configuration()
                    ).iter()
                    ).for_each(|(deformation_gradient_i, res_deformation_gradient_i)|
                        deformation_gradient_i.iter()
                        .zip(res_deformation_gradient_i.iter())
                        .for_each(|(deformation_gradient_ij, res_deformation_gradient_ij)|
                            assert_eq_within_tols(
                                deformation_gradient_ij, res_deformation_gradient_ij
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
                    get_element_transformed().calculate_deformation_gradient(
                        &get_reference_coordinates_transformed()
                        .convert()
                    ).iter()
                    .enumerate()
                    .for_each(|(i, deformation_gradient_i)|
                        deformation_gradient_i.iter()
                        .enumerate()
                        .for_each(|(j, deformation_gradient_ij)|
                            if i == j
                            {
                                assert_eq_within_tols(deformation_gradient_ij, &1.0)
                            }
                            else
                            {
                                assert_eq_within_tols(deformation_gradient_ij, &0.0)
                            }
                        )
                    )
                }
            }
        }
        mod gradient_vectors
        {
            use super::*;
            #[test]
            fn get<'a>()
            {
                $element::<'a, $constitutive_model<'a>>::calculate_gradient_vectors(
                    &get_reference_coordinates()
                ).iter().zip((
                    get_element().get_gradient_vectors()
                ).iter())
                .for_each(|(gradient_vector, res_gradient_vector)|
                    gradient_vector.iter()
                    .zip(res_gradient_vector.iter())
                    .for_each(|(gradient_vector_i, res_gradient_vector_i)|
                        assert_eq_within_tols(gradient_vector_i, res_gradient_vector_i)
                    )
                )
            }
            #[test]
            fn objectivity<'a>()
            {
                $element::<'a, $constitutive_model<'a>>::calculate_gradient_vectors(
                    &get_reference_coordinates()
                ).iter().zip((
                    get_rotation_reference_configuration().transpose() *
                    $element::<'a, $constitutive_model<'a>>::calculate_gradient_vectors(
                        &get_reference_coordinates_transformed()
                    )
                ).iter())
                .for_each(|(gradient_vector, res_gradient_vector)|
                    gradient_vector.iter()
                    .zip(res_gradient_vector.iter())
                    .for_each(|(gradient_vector_i, res_gradient_vector_i)|
                        assert_eq_within_tols(gradient_vector_i, res_gradient_vector_i)
                    )
                )
            }
        }
        mod standard_gradient_operator
        {
            use super::*;
            #[test]
            fn partition_of_unity<'a>()
            {
                let mut sum = [0.0_f64; 3];
                $element::<'a, $constitutive_model<'a>>::calculate_standard_gradient_operator()
                .iter()
                .for_each(|row|
                    row.iter()
                    .zip(sum.iter_mut())
                    .for_each(|(entry, sum_i)|
                        *sum_i += entry
                    )
                );
                sum.iter()
                .for_each(|sum_i|
                    assert_eq!(sum_i, &0.0)
                )
            }
        }
        #[test]
        fn size()
        {
            // really only for hyperelastic constitutive models (no state variables)
            assert_eq!(
                std::mem::size_of::<$element::<$constitutive_model>>(),
                std::mem::size_of::<$constitutive_model>()
                + std::mem::size_of::<GradientVectors<N>>()
            )
        }
    }
}
pub(crate) use test_linear_finite_element_with_constitutive_model;