use super::super::test::*;
use super::*;

use_elastic_macros!();

test_solid_hyperelastic_constitutive_model!(
    NeoHookean,
    NEOHOOKEANPARAMETERS,
    NeoHookean::new(NEOHOOKEANPARAMETERS)
);

test_solve_uniaxial!(NeoHookean::new(NEOHOOKEANPARAMETERS));
