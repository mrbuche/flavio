use super::
{
    ArrudaBoyce,
    super::test::
    {
        ARRUDABOYCEPARAMETERS,
        test_hyperelastic_constitutive_model
    }
};

test_hyperelastic_constitutive_model!(
    ArrudaBoyce,
    ARRUDABOYCEPARAMETERS,
    ArrudaBoyce::new(ARRUDABOYCEPARAMETERS)
);

#[test]
fn get_number_of_links()
{
    assert_eq!(
        &ARRUDABOYCEPARAMETERS[2],
        ArrudaBoyce::new(ARRUDABOYCEPARAMETERS).get_number_of_links()
    )
}
