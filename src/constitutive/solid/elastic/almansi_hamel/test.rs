use super::*;
use super::super::test::*;

test_solid_elastic_constitutive_model!(
    AlmansiHamel, ALMANSIHAMELPARAMETERS,
    AlmansiHamel::new(ALMANSIHAMELPARAMETERS)
);
