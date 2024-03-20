macro_rules! test_composite_element
{
    ($element: ident) =>
    {
        crate::fem::block::element::composite::test::setup_for_composite_elements!($element);
        crate::fem::block::element::test::test_finite_element!($element);
    }
}
pub(crate) use test_composite_element;

macro_rules! setup_for_composite_elements
{
    ($element: ident) =>
    {
        use crate::
        {
            constitutive::solid::elastic::AlmansiHamel,
            mechanics::test::
            {
                get_deformation_gradient,
                get_deformation_gradient_rate
            },
            test::assert_eq_within_tols
        };
        fn get_coordinates() -> NodalCoordinates<N>
        {
            get_deformation_gradient() * get_reference_coordinates()
        }
        fn get_velocities() -> NodalVelocities<N>
        {
            get_deformation_gradient_rate() * get_reference_coordinates()
        }
        mod jacobians
        {
            use super::*;
            mod undeformed
            {
                use super::*;
                #[test]
                fn jacobians()
                {
                    $element::<AlmansiHamel>::calculate_jacobians_and_parametric_gradient_operators(
                        &get_reference_coordinates()
                    ).0.iter()
                    .for_each(|jacobian|
                        assert_eq!(jacobian, &1.0)
                    )
                }
                #[test]
                fn scaled_composite_jacobians()
                {
                    $element::<AlmansiHamel>::calculate_scaled_composite_jacobian_at_integration_points(
                        &get_reference_coordinates()
                    ).iter()
                    .for_each(|scaled_composite_jacobian|
                        assert_eq!(scaled_composite_jacobian, &INTEGRATION_WEIGHT)
                    )
                }
            }
        }
        mod partition_of_unity
        {
            use super::*;
            #[test]
            fn shape_functions()
            {
                $element::<AlmansiHamel>::calculate_shape_functions_at_integration_points().iter()
                .for_each(|shape_functions|
                    assert_eq!(shape_functions.iter().sum::<Scalar>(), 1.0)
                )
            }
            #[test]
            fn standard_gradient_operators()
            {
                let mut sum = [0.0_f64; 3];
                $element::<AlmansiHamel>::calculate_standard_gradient_operators().iter()
                .for_each(|standard_gradient_operator|{
                    standard_gradient_operator.iter()
                    .for_each(|row|
                        row.iter()
                        .zip(sum.iter_mut())
                        .for_each(|(entry, sum_i)|
                            *sum_i += entry
                        )
                    );
                    sum.iter()
                    .for_each(|sum_i|
                        assert_eq_within_tols(sum_i, &0.0)
                    )
                })
            }
        }
        #[test]
        fn normalized_projection_matrix()
        {
            $element::<AlmansiHamel>::calculate_shape_function_integrals_products()
            .iter().map(|dummy| dummy * 1.0).sum::<TensorRank2<4, 9, 9>>().iter()
            .zip($element::<AlmansiHamel>::calculate_inverse_normalized_projection_matrix()
            .inverse().iter())
            .for_each(|(sum_i, projection_matrix_i)|
                sum_i.iter()
                .zip(projection_matrix_i.iter())
                .for_each(|(sum_ij, projection_matrix_ij)|
                    assert_eq_within_tols(sum_ij, projection_matrix_ij)
                )
            )
        }
        #[test]
        fn size()
        {
            assert_eq!(
                std::mem::size_of::<$element::<AlmansiHamel>>(),
                std::mem::size_of::<[AlmansiHamel; G]>()
                + std::mem::size_of::<ProjectedGradientVectors<G, N>>()
                + std::mem::size_of::<Scalars<G>>()
            )
        }
        #[test]
        fn todo()
        {
            todo!("Test this across constitutive models, and add similar tests to linear/test.rs.")
        }
        crate::fem::block::element::test::setup_for_element_tests_any_element!($element);
    }
}
pub(crate) use setup_for_composite_elements;
