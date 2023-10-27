use crate::constitutive::
{
    ConstitutiveModel,
    test::NEOHOOKEANPARAMETERS,
    hyperelastic::test::test_hyperelastic_constitutive_model
};
use super::NeoHookeanModel;

test_hyperelastic_constitutive_model!(
    NeoHookeanModel,
    NEOHOOKEANPARAMETERS
);
