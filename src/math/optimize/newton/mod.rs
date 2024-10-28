#[cfg(test)]
mod test;

use super::{
    super::{Tensor, TensorRank0},
    Optimize, OptimizeError,
};
use crate::ABS_TOL;
use std::{fmt, ops::Div};

/// Newton's method.
#[derive(Debug)]
pub struct Newton {
    /// Absolute error tolerance.
    pub abs_tol: TensorRank0,
    /// Maximum number of steps.
    pub max_steps: usize,
}

impl Default for Newton {
    fn default() -> Self {
        Self {
            abs_tol: ABS_TOL,
            max_steps: 100,
        }
    }
}

impl<X, J, H> Optimize<X, J, H> for Newton
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
    ) -> Result<X, OptimizeError> {
        let mut residual = jacobian(&initial_guess);
        let mut residual_norm = residual.norm();
        let mut solution = initial_guess;
        let mut steps = 0;
        while residual_norm >= self.abs_tol {
            if steps >= self.max_steps {
                return Err(OptimizeError::MaximumStepsReached(
                    steps,
                    format!("{:?}", &self),
                ));
            } else {
                residual = jacobian(&solution);
                residual_norm = residual.norm();
                solution -= residual / hessian(&solution);
                steps += 1;
            }
        }
        if hessian(&solution).is_positive_definite() {
            Ok(solution)
        } else {
            Err(OptimizeError::NotMinimum(
                format!("{}", solution),
                format!("{:?}", &self),
            ))
        }
    }
}
