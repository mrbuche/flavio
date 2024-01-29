use super::
{
    SaintVenantKirchoffModel,
    super::test::
    {
        SAINTVENANTKIRCHOFFPARAMETERS,
        test_thermohyperelastic_constitutive_model
    }
};

test_thermohyperelastic_constitutive_model!(
    SaintVenantKirchoffModel,
    SAINTVENANTKIRCHOFFPARAMETERS,
    SaintVenantKirchoffModel::new(SAINTVENANTKIRCHOFFPARAMETERS)
);
