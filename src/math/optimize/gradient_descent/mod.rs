#[cfg(test)]
mod test;

use super::{
    super::{Tensor, TensorRank0},
    FirstOrder, OptimizeError,
};
use crate::ABS_TOL;

/// The method of gradient descent.
#[derive(Debug)]
pub struct GradientDescent {
    /// Absolute error tolerance.
    pub abs_tol: TensorRank0,
    /// Maximum number of steps.
    pub max_steps: usize,
}

impl Default for GradientDescent {
    fn default() -> Self {
        Self {
            abs_tol: ABS_TOL,
            max_steps: 100,
        }
    }
}

impl<X> FirstOrder<X> for GradientDescent
where
    X: Tensor,
{
    fn minimize(&self, jacobian: impl Fn(&X) -> X, initial_guess: X) -> Result<X, OptimizeError> {
        //
        // How to choose short (below, dx*dg/dg*dg) or long (dx*dx/dx*dg) steps?
        // Or even allow different options for calculating step size?
        //
        let mut residual;
        let mut residual_change = X::zero();
        let mut solution = initial_guess;
        let mut solution_previous = X::zero();
        let mut step_size = 0.1;
        let mut step_trial;
        for _ in 0..self.max_steps {
            residual = jacobian(&solution);
            if residual.norm() < self.abs_tol {
                return Ok(solution);
            } else {
                residual_change -= &residual;
                step_trial = residual_change
                    .full_contraction(&(solution_previous - &solution))
                    .abs()
                    / residual_change.norm_squared();
                if step_trial > 0.0 && !step_trial.is_nan() {
                    step_size = step_trial
                } else {
                    step_size *= 1.1
                }
                residual_change = residual.copy();
                solution_previous = solution.copy();
                solution -= residual * step_size;
            }
        }
        Err(OptimizeError::MaximumStepsReached(
            self.max_steps,
            format!("{:?}", &self),
        ))
    }
}
