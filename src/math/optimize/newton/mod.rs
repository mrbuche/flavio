#[cfg(test)]
mod test;

use super::{
    super::{Tensor, TensorRank0},
    Optimize, OptimizeError,
};
use crate::ABS_TOL;
use std::ops::Div;

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
            max_steps: 1_000,
        }
    }
}

impl<H, J, X> Optimize<H, J, X> for Newton
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
    ) -> Result<X, OptimizeError> {
        let mut residual;
        let mut solution = initial_guess;
        for _ in 0..self.max_steps {
            residual = jacobian(&solution);
            if residual.norm() < self.abs_tol {
                if hessian(&solution).is_positive_definite() {
                    return Ok(solution);
                } else {
                    return Err(OptimizeError::NotMinimum(
                        format!("{}", solution),
                        format!("{:?}", &self),
                    ));
                }
            } else {
                solution -= residual / hessian(&solution);
            }
        }
        Err(OptimizeError::MaximumStepsReached(
            self.max_steps,
            format!("{:?}", &self),
        ))
    }
}
