use crate::constitutive::
{
    test::ALMANSIHAMELPARAMETERS,
    hyperelastic::test::test_hyperelastic_constitutive_model
};
use super::AlmansiHamelModel;

test_hyperelastic_constitutive_model!(
    AlmansiHamelModel,
    ALMANSIHAMELPARAMETERS,
    AlmansiHamelModel::new(ALMANSIHAMELPARAMETERS)
);
