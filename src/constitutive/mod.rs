//! Constitutive model library.

#[cfg(test)]
pub mod test;

pub mod thermal;
pub mod solid;

use crate::mechanics::Scalar;

/// Array of constitutive model parameters.
pub type ConstitutiveModelParameters<'a> = &'a [Scalar];

/// Required methods for constitutive models.
pub trait ConstitutiveModel<'a>
{
    /// Constructs and returns a new constitutive model.
    fn new(parameters: ConstitutiveModelParameters<'a>) -> Self;
}
