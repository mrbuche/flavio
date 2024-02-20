use super::*;
use super::super::test::*;

use_elastic_macros!();

test_solid_hyperelastic_constitutive_model!(
    ArrudaBoyce, ARRUDABOYCEPARAMETERS,
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