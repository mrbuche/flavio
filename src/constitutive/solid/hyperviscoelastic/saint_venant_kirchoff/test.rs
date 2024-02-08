use super::
{
    SaintVenantKirchoff,
    super::test::
    {
        SAINTVENANTKIRCHOFFPARAMETERS,
        test_hyperviscoelastic_constitutive_model
    }
};

test_hyperviscoelastic_constitutive_model!(
    SaintVenantKirchoff,
    SAINTVENANTKIRCHOFFPARAMETERS,
    SaintVenantKirchoff::new(SAINTVENANTKIRCHOFFPARAMETERS)
);