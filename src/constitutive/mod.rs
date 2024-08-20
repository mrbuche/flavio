//! Constitutive model library.

#[cfg(test)]
pub mod test;

pub mod cohesive;
pub mod fluid;
pub mod hybrid;
pub mod multiphysics;
pub mod solid;
pub mod thermal;

use crate::{
    get_defeat_message,
    mechanics::{DeformationGradient, Scalar},
};
use std::fmt;

/// Array of constitutive model parameters.
pub type Parameters<'a> = &'a [Scalar];

/// Required methods for constitutive models.
pub trait Constitutive<'a> {
    /// Constructs and returns a new constitutive model.
    fn new(parameters: Parameters<'a>) -> Self;
}

/// Error message prefix for constitutive model errors.
pub const CONSTITUTIVE_MODEL_ERROR: &str = "\x1b[91mConstitutive model error.\x1b[0\n";

/// Possible errors encountered in constitutive models.
pub enum ConstitutiveError {
    Custom(String, DeformationGradient, String),
    InvalidJacobian(Scalar, DeformationGradient, String),
}

/// Debug implementation for constitutive model errors.
impl fmt::Debug for ConstitutiveError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let error = match self {
            Self::Custom(message, deformation_gradient, constitutive_model) => format!(
                "{}\n\
                 From deformation gradient: {}.\n\
                 In constitutive model: {}.",
                message, deformation_gradient, constitutive_model
            ),
            Self::InvalidJacobian(jacobian, deformation_gradient, constitutive_model) => {
                format!(
                    "Invalid Jacobian: {:.6e}.\n\
                     From deformation gradient: {}.\n\
                     In constitutive model: {}.",
                    jacobian, deformation_gradient, constitutive_model
                )
            }
        };
        write!(
            f,
            "\x1b[91m{}\n\x1b[0;2;31m{}\x1b[0m\n",
            error,
            get_defeat_message()
        )
    }
}

/// Display implementation for constitutive model errors.
impl fmt::Display for ConstitutiveError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let error = match self {
            Self::Custom(message, deformation_gradient, constitutive_model) => format!(
                "{}\n\
                 From deformation gradient: {}.\n\
                 In constitutive model: {}.",
                message, deformation_gradient, constitutive_model
            ),
            Self::InvalidJacobian(jacobian, deformation_gradient, constitutive_model) => {
                format!(
                    "Invalid Jacobian: {:.6e}.\n\
                     From deformation gradient: {}.\n\
                     In constitutive model: {}.",
                    jacobian, deformation_gradient, constitutive_model
                )
            }
        };
        write!(
            f,
            "\x1b[91m{}\n\x1b[0;2;31m{}\x1b[0m\n",
            error,
            get_defeat_message()
        )
    }
}
