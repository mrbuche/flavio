use super::
{
    NeoHookeanModel,
    super::test::
    {
        NEOHOOKEANPARAMETERS,
        test_hyperelastic_constitutive_model
    }
};

test_hyperelastic_constitutive_model!(
    NeoHookeanModel,
    NEOHOOKEANPARAMETERS,
    NeoHookeanModel::new(NEOHOOKEANPARAMETERS)
);
