use super::
{
    MooneyRivlin,
    super::test::
    {
        MOONEYRIVLINPARAMETERS,
        test_hyperelastic_constitutive_model
    }
};

test_hyperelastic_constitutive_model!(
    MooneyRivlin,
    MOONEYRIVLINPARAMETERS,
    MooneyRivlin::new(MOONEYRIVLINPARAMETERS)
);

#[test]
fn get_extra_modulus()
{
    assert_eq!(
        &MOONEYRIVLINPARAMETERS[2],
        MooneyRivlin::new(MOONEYRIVLINPARAMETERS).get_extra_modulus()
    )
}
