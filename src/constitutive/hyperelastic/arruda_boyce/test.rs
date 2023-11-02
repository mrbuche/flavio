use crate::constitutive::
{
    test::ARRUDABOYCEPARAMETERS,
    hyperelastic::test::test_hyperelastic_constitutive_model
};
use super::ArrudaBoyceModel;

test_hyperelastic_constitutive_model!(
    ArrudaBoyceModel,
    ARRUDABOYCEPARAMETERS,
    ArrudaBoyceModel::new(ARRUDABOYCEPARAMETERS)
);

#[test]
fn get_number_density()
{
    assert_eq!(
        &ARRUDABOYCEPARAMETERS[2],
        ArrudaBoyceModel::new(ARRUDABOYCEPARAMETERS).get_number_density()
    )
}

#[test]
fn get_number_of_links()
{
    assert_eq!(
        &ARRUDABOYCEPARAMETERS[3],
        ArrudaBoyceModel::new(ARRUDABOYCEPARAMETERS).get_number_of_links()
    )
}
