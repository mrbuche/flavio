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

/// Array of constitutive model parameters.
pub type Parameters<'a> = &'a [Scalar];

/// Required methods for constitutive models.
pub trait Constitutive<'a>
{
    /// Constructs and returns a new constitutive model.
    fn new(parameters: Parameters<'a>) -> Self;
}

/// Possible errors encountered in constitutive models.
pub enum ConstitutiveError<'a> {
    Custom(String, &'a DeformationGradient, String),
    InvalidJacobianElastic(Scalar, &'a DeformationGradient, String),
    InvalidJacobianThermoelastic(Scalar, &'a DeformationGradient, &'a Scalar, String),
}

/// Display implementation for constitutive model errors.
impl<'a> fmt::Display for ConstitutiveError<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let error = match self {
            Self::Custom(message, deformation_gradient, model) => format!(
                "{}\n\
                 From deformation gradient: {}.\n\
                 In constitutive model: {}.",
                message, deformation_gradient, model
            ),
            Self::InvalidJacobianElastic(jacobian, deformation_gradient, model) => format!(
                "Invalid Jacobian: {:.6e}.\n\
                 From deformation gradient: {}.\n\
                 In constitutive model: {}.",
                jacobian, deformation_gradient, model
            ),
            Self::InvalidJacobianThermoelastic(jacobian, deformation_gradient, temperature, model) => format!(
                "Invalid Jacobian: {:.6e}.\n\
                 From deformation gradient: {}.\n\
                 For temperature: {:.6e}.\n\
                 In constitutive model: {}.",
                jacobian, deformation_gradient, temperature, model
            ),
        };
        write!(f, "\x1b[91m{}\n\x1b[0;2;31m{}\x1b[0m\n", error, get_error_message())
    }
}
