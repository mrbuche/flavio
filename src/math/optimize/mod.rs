#[cfg(test)]
mod test;

mod newton;

use super::Tensor;
use crate::get_defeat_message;
use std::{fmt, ops::Div};

// Newton doesn't need the objective function,
// so whether it's actually an optimization or root finding problem,
// the i/o is the same, so it can live in optimize safely.

/// ???
pub trait Optimize<X, J, H>
where
    Self: fmt::Debug,
    X: Tensor,
    J: Tensor + Div<H, Output = X>,
    H: Tensor,
{
    fn minimize(
        &self,
        jacobian: impl Fn(&X) -> J,
        hessian: impl Fn(&X) -> H,
        initial_guess: X,
    ) -> Result<X, OptimizeError>;
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
