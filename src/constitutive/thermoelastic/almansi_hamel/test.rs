use super::
{
    AlmansiHamelModel,
    ThermoelasticConstitutiveModel,
    super::test::
    {
        ALMANSIHAMELPARAMETERS,
        test_thermoelastic_constitutive_model,
        test_thermoelastic_only_constitutive_model_constructed
    }
};

test_thermoelastic_constitutive_model!(
    AlmansiHamelModel,
    ALMANSIHAMELPARAMETERS,
    AlmansiHamelModel::new(ALMANSIHAMELPARAMETERS)
);

test_thermoelastic_only_constitutive_model_constructed!(
    AlmansiHamelModel::new(ALMANSIHAMELPARAMETERS)
);