#[cfg(test)]
mod test;

use super::*;

pub struct CompositeHyperelasticConstitutiveModelMultiplicativeDecomposition<C1, C2>
{
    hyperelastic_constitutive_model_1: C1,
    hyperelastic_constitutive_model_2: C2
}
