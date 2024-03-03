macro_rules! test_linear_surface_element
{
    ($element: ident) =>
    {
        crate::fem::block::element::test::setup_for_surface_elements!($element);
        crate::fem::block::element::linear::test::test_linear_element_inner!($element);
        crate::fem::block::element::linear::surface::test::test_linear_surface_element_inner!($element);
    }
}
pub(crate) use test_linear_surface_element;

macro_rules! test_linear_surface_element_inner
{
    ($element: ident) =>
    {
        mod linear_surface_element
        {
            use crate::
            {
                fem::block::element::linear::surface::test::
                {
                    test_linear_surface_element_with_constitutive_model
                }
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
                    test_linear_surface_element_with_constitutive_model!($element, AlmansiHamel, ALMANSIHAMELPARAMETERS);
                }
                #[test]
                fn other_constitutive_models()
                {
                    todo!()
                }
            }
        }
    }
}
pub(crate) use test_linear_surface_element_inner;

macro_rules! test_linear_surface_element_with_constitutive_model
{
    ($element: ident, $constitutive_model: ident, $constitutive_model_parameters: ident) =>
    {
        mod normal
        {
            use super::*;
            mod deformed
            {
                use super::*;
                #[test]
                fn finite_difference()
                {
                    todo!("test the normal gradients")
                }
                #[test]
                fn normal()
                {
                    todo!("normal to surface basis")
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
                fn normal()
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
        mod normal_gradients
        {
            use super::*;
            mod deformed
            {
                #[test]
                fn objectivity()
                {
                    todo!()
                }
            }
            mod undeformed
            {
                #[test]
                fn objectivity()
                {
                    todo!()
                }
            }
        }
        mod normal_rate
        {
            use super::*;
            mod deformed
            {
                use super::*;
                #[test]
                fn finite_difference()
                {
                    get_normal_rate(true).iter()
                    .zip(get_normal_rate_from_finite_difference().iter())
                    .for_each(|(normal_rate_i, fd_normal_rate_i)|
                        assert!(
                            (normal_rate_i/fd_normal_rate_i - 1.0).abs() < EPSILON
                        )
                    )
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
                    get_normal_rate(false).iter()
                    .zip(get_normal_rate_from_finite_difference().iter())
                    .for_each(|(normal_rate_i, fd_normal_rate_i)|
                        assert!(
                            (normal_rate_i/fd_normal_rate_i - 1.0).abs() < EPSILON ||
                            normal_rate_i.abs() < EPSILON
                        )
                    )
                }
                #[test]
                fn objectivity()
                {
                    todo!()
                }
            }
        }
        mod reference_normal
        {
            use super::*;
            mod deformed
            {
                use super::*;
                #[test]
                fn normal()
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
                fn normal()
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
pub(crate) use test_linear_surface_element_with_constitutive_model;