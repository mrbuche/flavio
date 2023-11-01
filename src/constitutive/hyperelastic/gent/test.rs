use crate::constitutive::
{
    test::GENTPARAMETERS,
    hyperelastic::test::test_hyperelastic_constitutive_model
};
use super::GentModel;

test_hyperelastic_constitutive_model!(
    GentModel,
    GENTPARAMETERS,
    GentModel::new(GENTPARAMETERS)
);

#[test]
fn get_extensibility()
{
    assert_eq!(&GENTPARAMETERS[2], GentModel::new(GENTPARAMETERS).get_extensibility())
}
