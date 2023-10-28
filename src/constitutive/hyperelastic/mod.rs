#[cfg(test)]
mod test;

pub mod gent;
pub mod neo_hookean;

use super::*;

pub trait HyperelasticConstitutiveModel<'a>: ConstitutiveModel<'a> {}
