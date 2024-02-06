//! Constitutive model library.

#[cfg(test)]
pub mod test;

pub mod multiphysics;
pub mod thermal;
pub mod solid;

use crate::mechanics::Scalar;

/// Array of constitutive model parameters.
pub type Parameters<'a> = &'a [Scalar];

/// Required methods for constitutive models.
pub trait Constitutive<'a>
{
    /// Constructs and returns a new constitutive model.
    fn new(parameters: Parameters<'a>) -> Self;
}
