use super::
{
    NeoHookean,
    super::test::
    {
        NEOHOOKEANPARAMETERS,
        test_hyperelastic_constitutive_model
    }
};

test_hyperelastic_constitutive_model!(
    NeoHookean,
    NEOHOOKEANPARAMETERS,
    NeoHookean::new(NEOHOOKEANPARAMETERS)
);
