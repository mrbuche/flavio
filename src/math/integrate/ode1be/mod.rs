#[cfg(test)]
mod test;

use super::{
    super::{Tensor, TensorRank0, TensorRank0List, Tensors},
    Implicit, IntegrationError, OdeSolver,
};
use crate::{ABS_TOL, REL_TOL};
use std::ops::{Div, Mul, Sub};

/// Implicit, single-stage, first-order, variable-step, Runge-Kutta method (the backward Euler method).
#[derive(Debug)]
pub struct Ode1be {
    /// Absolute error tolerance.
    pub abs_tol: TensorRank0,
    /// Multiplying factor when decreasing time steps.
    pub dec_fac: TensorRank0,
    /// Multiplying factor when increasing time steps.
    pub inc_fac: TensorRank0,
    /// Maximum number of Newton steps MOVE THIS FIELD AND ANY ERROR HANDLING TO DEDICATED NEWTON SOLVER.
    pub max_steps: usize,
    /// Relative error tolerance.
    pub rel_tol: TensorRank0,
}

impl Default for Ode1be {
    fn default() -> Self {
        Self {
            abs_tol: ABS_TOL,
            dec_fac: 0.5,
            inc_fac: 1.1,
            max_steps: 1_000,
            rel_tol: REL_TOL,
        }
    }
}

impl<Y, J, U, const W: usize> Implicit<Y, J, U, W> for Ode1be
where
    Y: Tensor + Div<J, Output = Y>,
    for<'a> &'a Y: Mul<TensorRank0, Output = Y> + Sub<&'a Y, Output = Y>,
    J: Tensor,
    U: Tensors<Item = Y>,
{
    fn integrate(
        &self,
        function: impl Fn(&TensorRank0, &Y) -> Y,
        jacobian: impl Fn(&TensorRank0, &Y) -> J,
        initial_time: TensorRank0,
        initial_condition: Y,
        evaluation_times: &TensorRank0List<W>,
    ) -> Result<U, IntegrationError<W>> {
        let mut e;
        let mut k_1 = function(&initial_time, &initial_condition);
        let mut k_2;
        let mut solution = U::zero();
        let mut t_trial;
        let mut y_trial;
        let mut residual;
        let mut residual_norm;
        let mut steps;
        let identity = J::identity();
        {
            let (mut eval_times, mut dt, mut t, mut y, mut y_sol) = self.setup(
                initial_time,
                initial_condition,
                evaluation_times,
                &mut solution,
            )?;
            while eval_times.peek().is_some() {
                //
                // what about making a dedicated Newton solver implementation?
                // with max steps, tolerances (ABS_TOL below), etc. as an attributes?
                //
                k_2 = Y::zero();
                residual_norm = 1.0;
                steps = 0;
                t_trial = t + dt;
                y_trial = y.copy();
                while residual_norm >= ABS_TOL {
                    if steps >= self.max_steps {
                        return Err(IntegrationError::MaximumStepsReached(
                            steps,
                            format!("{:?}", &self),
                        ));
                    } else {
                        println!("ADFGHGFDDFGHGF {:?}", steps);
                        k_2 = function(&t_trial, &y_trial);
                        residual = &y_trial - &y - &(&k_2 * dt);
                        residual_norm = residual.norm();
                        y_trial += residual / (jacobian(&t_trial, &y_trial) * dt - &identity);
                        //
                        // (1) you dont want to reset the trial guess after evaluating the residual, in case it's below the tolerance
                        // (2) why arent we taking more than one netwon step?
                        //
                        steps += 1;
                    }
                }
                e = ((&k_2 - &k_1) * (dt / 2.0)).norm();
                if e < self.abs_tol || e / y_trial.norm() < self.rel_tol {
                    while let Some(eval_time) = eval_times.next_if(|&eval_time| t > eval_time) {
                        *y_sol.next().ok_or("not ok")? =
                            (&y_trial - &y) / dt * (eval_time - t) + &y;
                    }
                    k_1 = k_2;
                    t += dt;
                    dt *= self.inc_fac;
                    y = y_trial;
                } else {
                    dt *= self.dec_fac;
                }
            }
        }
        Ok(solution)
    }
}
