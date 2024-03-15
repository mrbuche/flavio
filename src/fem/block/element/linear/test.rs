macro_rules! test_linear_element
{
    ($element: ident) =>
    {
        crate::fem::block::element::test::setup_for_elements!($element);
        crate::fem::block::element::linear::test::test_linear_element_inner!($element);
    }
}
pub(crate) use test_linear_element;

macro_rules! test_linear_element_inner
{
    ($element: ident) =>
    {
        crate::fem::block::element::test::test_finite_element!($element);
        mod linear_element
        {
            use crate::
            {
                fem::block::element::linear::test::
                {
                    test_linear_element_with_constitutive_model
                },
                math::Convert,
                test::assert_eq_within_tols
            };
            use super::*;
            mod elastic
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
                mod almansi_hamel
                {
                    use super::*;
                    test_linear_element_with_constitutive_model!($element, AlmansiHamel, ALMANSIHAMELPARAMETERS);
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
                    test_linear_element_with_constitutive_model!($element, ArrudaBoyce, ARRUDABOYCEPARAMETERS);
                }
                mod fung
                {
                    use super::*;
                    test_linear_element_with_constitutive_model!($element, Fung, FUNGPARAMETERS);
                }
                mod gent
                {
                    use super::*;
                    test_linear_element_with_constitutive_model!($element, Gent, GENTPARAMETERS);
                }
                mod mooney_rivlin
                {
                    use super::*;
                    test_linear_element_with_constitutive_model!($element, MooneyRivlin, MOONEYRIVLINPARAMETERS);
                }
                mod neo_hookean
                {
                    use super::*;
                    test_linear_element_with_constitutive_model!($element, NeoHookean, NEOHOOKEANPARAMETERS);
                }
                mod saint_venant_kirchoff
                {
                    use super::*;
                    test_linear_element_with_constitutive_model!($element, SaintVenantKirchoff, SAINTVENANTKIRCHOFFPARAMETERS);
                }
                mod yeoh
                {
                    use super::*;
                    test_linear_element_with_constitutive_model!($element, Yeoh, YEOHPARAMETERS);
                }
            }
            mod elastic_hyperviscous
            {
                use crate::
                {
                    constitutive::solid::elastic_hyperviscous::
                    {
                        AlmansiHamel,
                        test::ALMANSIHAMELPARAMETERS
                    }
                };
                use super::*;
                mod almansi_hamel
                {
                    use super::*;
                    test_linear_element_with_constitutive_model!($element, AlmansiHamel, ALMANSIHAMELPARAMETERS);
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
                    test_linear_element_with_constitutive_model!($element, SaintVenantKirchoff, SAINTVENANTKIRCHOFFPARAMETERS);
                }
            }
        }
    }
}
pub(crate) use test_linear_element_inner;

macro_rules! test_linear_element_with_constitutive_model
{
    ($element: ident, $constitutive_model: ident, $constitutive_model_parameters: ident) =>
    {
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
        mod deformation_gradient
        {
            use super::*;
            mod deformed
            {
                use super::*;
                #[test]
                fn calculate()
                {
                    get_element().calculate_deformation_gradient(
                        &get_coordinates()
                    ).iter().zip(get_deformation_gradient().iter())
                    .for_each(|(calculated_deformation_gradient_i, deformation_gradient_i)|
                        calculated_deformation_gradient_i.iter()
                        .zip(deformation_gradient_i.iter())
                        .for_each(|(calculated_deformation_gradient_ij, deformation_gradient_ij)|
                            assert_eq_within_tols(
                                calculated_deformation_gradient_ij, deformation_gradient_ij
                            )
                        )
                    )
                }
                #[test]
                fn objectivity()
                {
                    get_element().calculate_deformation_gradient(
                        &get_coordinates()
                    ).iter().zip((
                        get_rotation_current_configuration().transpose() *
                        get_element_transformed().calculate_deformation_gradient(
                            &get_coordinates_transformed()
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
                fn calculate()
                {
                    get_element().calculate_deformation_gradient(
                        &get_reference_coordinates().convert()
                    ).iter().enumerate()
                    .for_each(|(i, calculated_deformation_gradient_i)|
                        calculated_deformation_gradient_i.iter().enumerate()
                        .for_each(|(j, calculated_deformation_gradient_ij)|
                            if i == j
                            {
                                assert_eq_within_tols(calculated_deformation_gradient_ij, &1.0)
                            }
                            else
                            {
                                assert_eq_within_tols(calculated_deformation_gradient_ij, &0.0)
                            }
                        )
                    )
                }
                #[test]
                fn objectivity()
                {
                    get_element_transformed().calculate_deformation_gradient(
                        &get_reference_coordinates_transformed().convert()
                    ).iter().enumerate()
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
        mod deformation_gradient_rate
        {
            use super::*;
            mod deformed
            {
                use super::*;
                #[test]
                fn calculate()
                {
                    get_element().calculate_deformation_gradient_rate(
                        &get_coordinates(), &get_velocities()
                    ).iter().zip(get_deformation_gradient_rate().iter())
                    .for_each(|(calculated_deformation_gradient_rate_i, deformation_gradient_rate_i)|
                        calculated_deformation_gradient_rate_i.iter()
                        .zip(deformation_gradient_rate_i.iter())
                        .for_each(|(calculated_deformation_gradient_rate_ij, deformation_gradient_rate_ij)|
                            assert_eq_within_tols(
                                calculated_deformation_gradient_rate_ij, deformation_gradient_rate_ij
                            )
                        )
                    )
                }
                #[test]
                fn objectivity()
                {
                    get_element().calculate_deformation_gradient_rate(
                        &get_coordinates(), &get_velocities()
                    ).iter().zip((
                        get_rotation_current_configuration().transpose() *(
                            get_element_transformed().calculate_deformation_gradient_rate(
                                &get_coordinates_transformed(), &get_velocities_transformed()
                            ) * get_rotation_reference_configuration() -
                            get_rotation_rate_current_configuration() *
                            get_element().calculate_deformation_gradient(
                                &get_coordinates()
                            )
                        )
                    ).iter()
                    ).for_each(|(deformation_gradient_rate_i, res_deformation_gradient_rate_i)|
                        deformation_gradient_rate_i.iter()
                        .zip(res_deformation_gradient_rate_i.iter())
                        .for_each(|(deformation_gradient_rate_ij, res_deformation_gradient_rate_ij)|
                            assert_eq_within_tols(
                                deformation_gradient_rate_ij, res_deformation_gradient_rate_ij
                            )
                        )
                    )
                }
            }
            mod undeformed
            {
                use super::*;
                #[test]
                fn calculate()
                {
                    get_element().calculate_deformation_gradient_rate(
                        &get_reference_coordinates().convert(),
                        &NodalVelocities::zero().convert()
                    ).iter()
                    .for_each(|calculated_deformation_gradient_rate_i|
                        calculated_deformation_gradient_rate_i.iter()
                        .for_each(|calculated_deformation_gradient_rate_ij|
                            assert_eq_within_tols(calculated_deformation_gradient_rate_ij, &0.0)
                        )
                    )
                }
                #[test]
                fn objectivity()
                {
                    get_element_transformed().calculate_deformation_gradient_rate(
                        &get_reference_coordinates_transformed().convert(),
                        &NodalVelocities::zero().convert()
                    ).iter()
                    .for_each(|deformation_gradient_rate_i|
                        deformation_gradient_rate_i.iter()
                        .for_each(|deformation_gradient_rate_ij|
                            assert_eq_within_tols(deformation_gradient_rate_ij, &0.0)
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
                $element::<$constitutive_model<'a>>::calculate_gradient_vectors(
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
                $element::<$constitutive_model<'a>>::calculate_gradient_vectors(
                    &get_reference_coordinates()
                ).iter().zip((
                    get_rotation_reference_configuration().transpose() *
                    $element::<$constitutive_model<'a>>::calculate_gradient_vectors(
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
                $element::<$constitutive_model<'a>>::calculate_standard_gradient_operator()
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
    }
}
pub(crate) use test_linear_element_with_constitutive_model;