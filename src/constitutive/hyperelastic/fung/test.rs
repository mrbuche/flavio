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

#[test]
fn get_extra_modulus()
{
    assert_eq!(
        &FUNGPARAMETERS[2],
        FungModel::new(FUNGPARAMETERS).get_extra_modulus()
    )
}

#[test]
fn get_exponent()
{
    assert_eq!(
        &FUNGPARAMETERS[3],
        FungModel::new(FUNGPARAMETERS).get_exponent()
    )
}
