use super::*;
use super::super::test::*;

use_elastic_macros!();

test_solid_hyperelastic_constitutive_model!(
    SaintVenantKirchoff, SAINTVENANTKIRCHOFFPARAMETERS,
    SaintVenantKirchoff::new(SAINTVENANTKIRCHOFFPARAMETERS)
);
