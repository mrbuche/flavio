use super::
{
    ArrudaBoyceModel,
    super::test::
    {
        ARRUDABOYCEPARAMETERS,
        test_hyperelastic_constitutive_model
    }
};

test_hyperelastic_constitutive_model!(
    ArrudaBoyceModel,
    ARRUDABOYCEPARAMETERS,
    ArrudaBoyceModel::new(ARRUDABOYCEPARAMETERS)
);

#[test]
fn get_number_of_links()
{
    assert_eq!(
        &ARRUDABOYCEPARAMETERS[2],
        ArrudaBoyceModel::new(ARRUDABOYCEPARAMETERS).get_number_of_links()
    )
}
