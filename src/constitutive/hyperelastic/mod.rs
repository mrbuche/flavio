#[cfg(test)]
mod test;

mod neo_hookean;

use super::*;

pub trait HyperelasticConstitutiveModel<'a>: ConstitutiveModel<'a> {}
