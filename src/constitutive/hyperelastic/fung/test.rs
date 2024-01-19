use crate::constitutive::
{
    test::FUNGPARAMETERS,
    hyperelastic::test::test_hyperelastic_constitutive_model
};
use super::FungModel;

test_hyperelastic_constitutive_model!(
    FungModel,
    FUNGPARAMETERS,
    FungModel::new(FUNGPARAMETERS)
);
