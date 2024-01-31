use super::
{
    SaintVenantKirchoffModel,
    super::test::
    {
        SAINTVENANTKIRCHOFFPARAMETERS,
        test_hyperelastic_constitutive_model
    }
};

test_hyperelastic_constitutive_model!(
    SaintVenantKirchoffModel,
    SAINTVENANTKIRCHOFFPARAMETERS,
    SaintVenantKirchoffModel::new(SAINTVENANTKIRCHOFFPARAMETERS)
);
