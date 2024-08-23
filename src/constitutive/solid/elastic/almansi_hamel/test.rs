use super::super::test::*;
use super::*;

test_solid_elastic_constitutive_model!(
    AlmansiHamel,
    ALMANSIHAMELPARAMETERS,
    AlmansiHamel::new(ALMANSIHAMELPARAMETERS)
);

test_solve_uniaxial_tension!(AlmansiHamel::new(ALMANSIHAMELPARAMETERS));
