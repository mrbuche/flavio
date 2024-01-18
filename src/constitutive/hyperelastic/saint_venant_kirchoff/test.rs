use crate::constitutive::
{
    test::SAINTVENANTKIRCHOFFPARAMETERS,
    hyperelastic::test::test_hyperelastic_constitutive_model
};
use super::SaintVenantKirchoffModel;

test_hyperelastic_constitutive_model!(
    SaintVenantKirchoffModel,
    SAINTVENANTKIRCHOFFPARAMETERS,
    SaintVenantKirchoffModel::new(SAINTVENANTKIRCHOFFPARAMETERS)
);
