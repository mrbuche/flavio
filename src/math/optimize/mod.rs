#[cfg(test)]
mod test;

mod newton;

use super::Tensor;
use crate::get_defeat_message;
use std::{fmt, ops::Div};

pub use newton::Newton;

/// Base trait for optimization algorithms.
pub trait Optimize<H, J, X>
where
    H: Tensor,
    J: Tensor + Div<H, Output = X>,
    X: Tensor,
{
    fn minimize(
        &self,
        jacobian: impl Fn(&X) -> J,
        hessian: impl Fn(&X) -> H,
        initial_guess: X,
    ) -> Result<X, OptimizeError>;
}

/// Possible optimization algorithms.
#[derive(Debug)]
pub enum Optimization {
    Newton(Newton),
}

/// Possible errors encountered when optimizing.
pub enum OptimizeError {
    MaximumStepsReached(usize, String),
    NotMinimum(String, String),
}

impl fmt::Debug for OptimizeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let error = match self {
            Self::MaximumStepsReached(steps, optimizer) => {
                format!(
                    "\x1b[1;91mMaximum number of steps ({}) reached.\x1b[0;91m\n\
                     In optimizer: {}.",
                    steps, optimizer
                )
            }
            Self::NotMinimum(solution, optimizer) => {
                format!(
                    "\x1b[1;91mThe obtained solution is not a minimum.\x1b[0;91m\n\
                     For solution: {}.\n\
                     In optimizer: {}.",
                    solution, optimizer
                )
            }
        };
        write!(
            f,
            "\n{}\n\x1b[0;2;31m{}\x1b[0m\n",
            error,
            get_defeat_message()
        )
    }
}
