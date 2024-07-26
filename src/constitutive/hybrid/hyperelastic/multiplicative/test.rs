use crate::constitutive::hybrid::
{
    Multiplicative,
    hyperelastic::test::test_hybrid_hyperelastic_constitutive_models
};

test_hybrid_hyperelastic_constitutive_models!(Multiplicative);
