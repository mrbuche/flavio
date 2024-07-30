macro_rules! test_hybrid_elastic_constitutive_models
{
    ($hybrid_type: ident) =>
    {
        use crate::
        {
            constitutive::
            {
                Constitutive,
                hybrid::Hybrid,
                solid::
                {
                    elastic::
                    {
                        AlmansiHamel,
                        Elastic,
                        test::*
                    },
                    hyperelastic::
                    {
                        NeoHookean,
                        test::NEOHOOKEANPARAMETERS
                    }
                }
            },
            math::TensorRank2Trait,
            mechanics::
            {
                CauchyTangentStiffness,
                DeformationGradient,
                FirstPiolaKirchoffTangentStiffness,
                SecondPiolaKirchoffTangentStiffness
            }
        };
        mod hybrid_1
        {
            use super::*;
            test_constructed_solid_constitutive_model!(
                $hybrid_type::construct(
                    AlmansiHamel::new(ALMANSIHAMELPARAMETERS),
                    NeoHookean::new(NEOHOOKEANPARAMETERS)
                )
            );
        }
        mod hybrid_nested_1
        {
            use super::*;
            test_constructed_solid_constitutive_model!(
                $hybrid_type::construct(
                    AlmansiHamel::new(ALMANSIHAMELPARAMETERS),
                    $hybrid_type::construct(
                        AlmansiHamel::new(ALMANSIHAMELPARAMETERS),
                        NeoHookean::new(NEOHOOKEANPARAMETERS)
                    )
                )
            );
        }
        mod hybrid_nested_2
        {
            use super::*;
            test_constructed_solid_constitutive_model!(
                $hybrid_type::construct(
                    $hybrid_type::construct(
                        NeoHookean::new(NEOHOOKEANPARAMETERS),
                        AlmansiHamel::new(ALMANSIHAMELPARAMETERS)
                    ),
                    $hybrid_type::construct(
                        NeoHookean::new(NEOHOOKEANPARAMETERS),
                        $hybrid_type::construct(
                            AlmansiHamel::new(ALMANSIHAMELPARAMETERS),
                            NeoHookean::new(NEOHOOKEANPARAMETERS)
                        )
                    )
                )
            );
        }
    }
}
pub(crate) use test_hybrid_elastic_constitutive_models;

macro_rules! test_hybrid_elastic_constitutive_models_no_tangents
{
    ($hybrid_type: ident) =>
    {
        use crate::
        {
            constitutive::
            {
                Constitutive,
                hybrid::Hybrid,
                solid::
                {
                    elastic::
                    {
                        AlmansiHamel,
                        Elastic,
                        test::*
                    },
                    hyperelastic::
                    {
                        NeoHookean,
                        test::NEOHOOKEANPARAMETERS
                    }
                }
            },
            math::TensorRank2Trait,
            mechanics::DeformationGradient
        };
        mod hybrid_1
        {
            use super::*;
            test_solid_constitutive_model_no_tangents!(
                $hybrid_type::construct(
                    AlmansiHamel::new(ALMANSIHAMELPARAMETERS),
                    NeoHookean::new(NEOHOOKEANPARAMETERS)
                )
            );
        }
    }
}
pub(crate) use test_hybrid_elastic_constitutive_models_no_tangents;
