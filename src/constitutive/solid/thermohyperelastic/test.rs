use crate::
{
    constitutive::solid::thermoelastic::test::ALMANSIHAMELPARAMETERS,
    mechanics::Scalar
};
pub const SAINTVENANTKIRCHOFFPARAMETERS: &[Scalar; 4] = &[ALMANSIHAMELPARAMETERS[0], ALMANSIHAMELPARAMETERS[1], ALMANSIHAMELPARAMETERS[2], ALMANSIHAMELPARAMETERS[3]];

macro_rules! test_thermohyperelastic_constitutive_model
{
    ($thermohyperelastic_constitutive_model: ident, $thermohyperelastic_constitutive_model_parameters: expr, $thermohyperelastic_constitutive_model_constructed: expr) =>
    {
        use crate::constitutive::solid::
        {
            thermoelastic::
            {
                ThermoelasticConstitutiveModel,
                test::test_thermoelastic_constitutive_model,
            },
            thermohyperelastic::
            {
                test::test_thermohyperelastic_constitutive_model_constructed
            }
        };
        test_thermoelastic_constitutive_model!($thermohyperelastic_constitutive_model, $thermohyperelastic_constitutive_model_parameters, $thermohyperelastic_constitutive_model_constructed);
        test_thermohyperelastic_constitutive_model_constructed!($thermohyperelastic_constitutive_model_constructed);
    }
}
pub(crate) use test_thermohyperelastic_constitutive_model;
macro_rules! test_thermohyperelastic_constitutive_model_constructed
{
    ($thermohyperelastic_constitutive_model_constructed: expr) =>
    {
    }
}
pub(crate) use test_thermohyperelastic_constitutive_model_constructed;