use super::
{
    MooneyRivlinModel,
    super::test::
    {
        MOONEYRIVLINPARAMETERS,
        test_hyperelastic_constitutive_model
    }
};

test_hyperelastic_constitutive_model!(
    MooneyRivlinModel,
    MOONEYRIVLINPARAMETERS,
    MooneyRivlinModel::new(MOONEYRIVLINPARAMETERS)
);

#[test]
fn get_extra_modulus()
{
    assert_eq!(
        &MOONEYRIVLINPARAMETERS[2],
        MooneyRivlinModel::new(MOONEYRIVLINPARAMETERS).get_extra_modulus()
    )
}
