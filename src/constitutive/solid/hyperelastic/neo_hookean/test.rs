use super::*;
use super::super::test::*;

use_elastic_macros!();

test_solid_hyperelastic_constitutive_model!(
    NeoHookean, NEOHOOKEANPARAMETERS,
    NeoHookean::new(NEOHOOKEANPARAMETERS)
);
