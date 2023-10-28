use crate::constitutive::
{
    test::MOONEYRIVLINPARAMETERS,
    hyperelastic::test::test_hyperelastic_constitutive_model
};
use super::MooneyRivlinModel;

test_hyperelastic_constitutive_model!(
    MooneyRivlinModel,
    MOONEYRIVLINPARAMETERS
);

#[test]
fn get_extra_modulus()
{
    assert_eq!(&MOONEYRIVLINPARAMETERS[2], MooneyRivlinModel::new(MOONEYRIVLINPARAMETERS).get_extra_modulus())
}
