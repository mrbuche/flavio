macro_rules! test_hybrid_hyperelastic_constitutive_models
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
                    elastic::Elastic,
                    hyperelastic::
                    {
                        ArrudaBoyce,
                        Fung,
                        Gent,
                        Hyperelastic,
                        MooneyRivlin,
                        NeoHookean,
                        SaintVenantKirchoff,
                        Yeoh,
                        test::*
                    }
                }
            },
            math::TensorRank2Trait,
            mechanics::
            {
                CauchyTangentStiffness,
                DeformationGradient,
                FirstPiolaKirchoffStress,
                FirstPiolaKirchoffTangentStiffness,
                SecondPiolaKirchoffTangentStiffness
            }
        };
        use_elastic_macros!();
        mod hybrid_1
        {
            use super::*;
            test_constructed_solid_hyperelastic_constitutive_model!(
                $hybrid_type::construct(
                    ArrudaBoyce::new(ARRUDABOYCEPARAMETERS),
                    Fung::new(FUNGPARAMETERS)
                )
            );
        }
        mod hybrid_2
        {
            use super::*;
            test_constructed_solid_hyperelastic_constitutive_model!(
                $hybrid_type::construct(
                    Gent::new(GENTPARAMETERS),
                    MooneyRivlin::new(MOONEYRIVLINPARAMETERS)
                )
            );
        }
        mod hybrid_nested_1
        {
            use super::*;
            test_constructed_solid_hyperelastic_constitutive_model!(
                $hybrid_type::construct(
                    NeoHookean::new(NEOHOOKEANPARAMETERS),
                    $hybrid_type::construct(
                        SaintVenantKirchoff::new(SAINTVENANTKIRCHOFFPARAMETERS),
                        Yeoh::new(YEOHPARAMETERS)
                    )
                )
            );
        }
        mod hybrid_nested_2
        {
            use super::*;
            test_constructed_solid_hyperelastic_constitutive_model!(
                $hybrid_type::construct(
                    $hybrid_type::construct(
                        Gent::new(GENTPARAMETERS),
                        MooneyRivlin::new(MOONEYRIVLINPARAMETERS)
                    ),
                    $hybrid_type::construct(
                        NeoHookean::new(NEOHOOKEANPARAMETERS),
                        $hybrid_type::construct(
                            SaintVenantKirchoff::new(SAINTVENANTKIRCHOFFPARAMETERS),
                            Yeoh::new(YEOHPARAMETERS)
                        )
                    )
                )
            );
        }
    }
}
pub(crate) use test_hybrid_hyperelastic_constitutive_models;
