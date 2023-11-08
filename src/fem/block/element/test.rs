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
            fem::block::element::test::test_finite_element_for_constitutive_model
        };
        mod gent
        {
            use super::*;
            test_finite_element_for_constitutive_model!($element, GentModel, GENTPARAMETERS);
        }
        mod mooney_rivlin
        {
            use super::*;
            test_finite_element_for_constitutive_model!($element, MooneyRivlinModel, NEOHOOKEANPARAMETERS);
        }
        mod neo_hookean
        {
            use super::*;
            test_finite_element_for_constitutive_model!($element, NeoHookeanModel, MOONEYRIVLINPARAMETERS);
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
        fn get_element<'a>() -> $element<'a, $constitutive_model<'a>>
        {
            $element::new(
                $constitutive_model_parameters,
                get_element_standard_reference_coordinates()
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
                    todo!()
                }
                #[test]
                fn positive()
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