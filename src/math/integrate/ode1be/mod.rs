#[cfg(test)]
mod test;

use super::{
    super::{Tensor, TensorRank0, TensorRank0List, Tensors},
    Implicit, IntegrationError, OdeSolver,
};
use crate::{ABS_TOL, REL_TOL};

/// Implicit, single-stage, first-order, fixed-step, Runge-Kutta method (the backward Euler method).
#[derive(Debug)]
pub struct Ode1be {
    /// Absolute error tolerance.
    pub abs_tol: TensorRank0,
    /// Relative error tolerance.
    pub rel_tol: TensorRank0,
}

impl Default for Ode1be {
    fn default() -> Self {
        Self {
            abs_tol: ABS_TOL,
            rel_tol: REL_TOL,
        }
    }
}

impl<Y, U, const W: usize> OdeSolver<Y, U, W> for Ode1be
where
    Y: Tensor,
    U: Tensors<Item = Y>,
{
}

impl<Y, J, U, const W: usize> Implicit<Y, J, U, W> for Ode1be
where
    Y: Tensor,
    for<'a> &'a Y: std::ops::Mul<TensorRank0, Output = Y>,
    U: Tensors<Item = Y>,
{
    fn integrate(
        &self,
        _function: impl Fn(&TensorRank0, &Y) -> Y,
        _jacobian: impl Fn(&TensorRank0, &Y) -> J,
        initial_time: TensorRank0,
        initial_condition: Y,
        evaluation_times: &TensorRank0List<W>,
    ) -> Result<U, IntegrationError<W>> {
        //
        // The timestep would be an input argument.
        // Do you even want to bother with fixed timestep methods?
        // Are there adaptive timestep BE methods out there?
        // Still call it ode1be?
        // Explicit Heun with Euler for error would be ode12, or Fehlberg RK1(2), try to be MATLAB consistent.
        // Would the Implicit flavor be ode12i or something?
        // Nah, probably ode12t since it ends up being trapezoidal.
        // What would be an implicit ode23? is that the ode23t in MATLAB?
        //
        // Can you even have embedded error for an ode1 ?
        // There are no lower order methods!
        // Seems like you need ode12 or higher.
        //
        // maybe skip ode1x
        // try to do what Julia does for Trapezoid() -> ode23t()
        //
        // do ode12 as Heun
        // do implicit as Trapezoidal (directly analogous)
        // what to call it? ode12t?
        // do you want to do that generalized alpha thing?
        //
        let mut solution = U::zero();
        {
            let _ = self.setup(
                // let (mut eval_times, mut dt, mut t, mut y, mut y_sol) = self.setup(
                initial_time,
                initial_condition,
                evaluation_times,
                &mut solution,
            )?;
        }
        Ok(solution)
    }
}
