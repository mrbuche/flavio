use super::*;
use super::super::test::*;

use_elastic_macros!();

test_solid_hyperelastic_constitutive_model!(
    Gent, GENTPARAMETERS,
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