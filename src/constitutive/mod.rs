//! Constitutive model library.

#[cfg(test)]
pub mod test;

pub mod cohesive;
pub mod fluid;
pub mod hybrid;
pub mod multiphysics;
pub mod solid;
pub mod thermal;

use crate::
{
    get_error_message,
    mechanics::{DeformationGradient, Scalar}
};
use std::fmt;

/// Possible errors encountered in constitutive models.
pub enum ConstitutiveError {
    InvalidJacobian(Scalar, DeformationGradient, String),
}

impl fmt::Display for ConstitutiveError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let error = match self {
            Self::InvalidJacobian(jacobian, deformation_gradient, model) => format!(
                "Invalid Jacobian: {}.\n\
                 From deformation gradient: {}.\n\
                 In constitutive model: {}.",
                jacobian, deformation_gradient, model
            )
        };
        write!(f, "\x1b[91m{}\n\x1b[0;2;31m{}\x1b[0m\n", error, get_error_message())
    }
}

/// Array of constitutive model parameters.
pub type Parameters<'a> = &'a [Scalar];

/// Required methods for constitutive models.
pub trait Constitutive<'a>
{
    /// Constructs and returns a new constitutive model.
    fn new(parameters: Parameters<'a>) -> Self;
}
