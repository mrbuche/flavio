use super::
{
    Gent,
    super::test::
    {
        GENTPARAMETERS,
        test_hyperelastic_constitutive_model
    }
};

test_hyperelastic_constitutive_model!(
    Gent,
    GENTPARAMETERS,
    Gent::new(GENTPARAMETERS)
);

#[test]
fn get_extensibility()
{
    assert_eq!(
        &GENTPARAMETERS[2],
        Gent::new(GENTPARAMETERS).get_extensibility()
    )
}
