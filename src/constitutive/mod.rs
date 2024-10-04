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

#[doc(hidden)]
pub const CONSTITUTIVE_MODEL_ERROR: &str = "\n\x1b[91mConstitutive model error.\x1b[0";

/// Possible errors encountered in constitutive models.
pub enum ConstitutiveError {
    Custom(String, DeformationGradient, String),
    InvalidJacobian(Scalar, DeformationGradient, String),
    SolveError,
}

impl fmt::Debug for ConstitutiveError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let error = match self {
            Self::Custom(message, deformation_gradient, constitutive_model) => format!(
                "\x1b[1;91m{}\x1b[0;91m\n\
                 From deformation gradient: {}.\n\
                 In constitutive model: {}.",
                message, deformation_gradient, constitutive_model
            ),
            Self::InvalidJacobian(jacobian, deformation_gradient, constitutive_model) => {
                format!(
                    "\x1b[1;91mInvalid Jacobian: {:.6e}.\x1b[0;91m\n\
                     From deformation gradient: {}.\n\
                     In constitutive model: {}.",
                    jacobian, deformation_gradient, constitutive_model
                )
            },
            Self::SolveError => "\x1b[1;91mThe maximum number of steps was reached before the tolerance was satisfied.\x1b[0;91m".to_string()
        };
        write!(
            f,
            "\n{}\n\x1b[0;2;31m{}\x1b[0m\n",
            error,
            get_defeat_message()
        )
    }
}

impl PartialEq for ConstitutiveError {
    fn eq(&self, other: &Self) -> bool {
        match self {
            Self::Custom(a, b, c) => match other {
                Self::Custom(d, e, f) => a == d && c == f,
                _ => false,
            },
            Self::InvalidJacobian(a, b, c) => match other {
                Self::InvalidJacobian(d, e, f) => a == d && c == f,
                _ => false,
            },
            Self::SolveError => todo!(),
        }
    }
}
