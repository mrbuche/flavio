#[cfg(test)]
mod test;

use super::{
    super::{Tensor, TensorRank0, TensorRank0List, Tensors},
    Explicit, IntegrationError,
};
use crate::{ABS_TOL, REL_TOL};

/// Explicit, three-stage, third-order, variable-step, Runge-Kutta method ([Bogacki and Shampine, 1989](https://doi.org/10.1016/0893-9659(89)90079-7)).
///
/// ```math
/// \frac{dy}{dt} = f(t, y)
/// ```
/// ```math
/// t_{n+1} = t_n + h
/// ```
/// ```math
/// k_1 = f(t_n, y_n)
/// ```
/// ```math
/// k_2 = f(t_n + \tfrac{1}{2} h, y_n + \tfrac{1}{2} h k_1)
/// ```
/// ```math
/// k_3 = f(t_n + \tfrac{3}{4} h, y_n + \tfrac{3}{4} h k_2)
/// ```
/// ```math
/// y_{n+1} = y_n + \frac{h}{9}\left(2k_1 + 3k_2 + 4k_3\right)
/// ```
/// ```math
/// k_4 = f(t_{n+1}, y_{n+1})
/// ```
/// ```math
/// e_{n+1} = \frac{h}{72}\left(-5k_1 + 6k_2 + 8k_3 - 9k_4\right)
/// ```
pub struct Ode23 {
    /// Absolute error tolerance.
    pub abs_tol: TensorRank0,
    /// Multiplying factor when decreasing timesteps.
    pub dec_fac: TensorRank0,
    /// Multiplying factor when increasing timesteps.
    pub inc_fac: TensorRank0,
    /// Relative error tolerance.
    pub rel_tol: TensorRank0,
}

impl Default for Ode23 {
    fn default() -> Self {
        Self {
            abs_tol: ABS_TOL,
            dec_fac: 0.5,
            inc_fac: 1.1,
            rel_tol: REL_TOL,
        }
    }
}

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

/// TODO: Should start time using t_0, leaving option for evaluation times not to include the initial time, but assert eval_times >= t_0.
///       Dont need t_end since what is important is last evaluation time.
///       But do need to enforce evaluation times are all greater than t_0.
///       And also that they are monotonic.
/// TODO: Need to do something to prevent large timesteps from going past more than one evaluation time.
impl Explicit for Ode23 {
    fn integrate<Y, U, const W: usize>(
        &self,
        function: impl Fn(&TensorRank0, &Y) -> Y,
        _t_0: TensorRank0,
        y_0: Y,
        evaluation_times: &TensorRank0List<W>,
    ) -> Result<U, IntegrationError>
    where
        Y: Tensor,
        for<'a> &'a Y: std::ops::Mul<TensorRank0, Output = Y>,
        U: Tensors<Item = Y>,
    {
        let mut error;
        let mut error_norm;
        let mut evaluation_time = evaluation_times.0.into_iter().peekable();
        let mut k_1;
        let mut k_2;
        let mut k_3;
        let mut k_4;
        let mut time = evaluation_time.next().ok_or("not ok")?;
        let mut timestep = evaluation_time.peek().ok_or("not ok")? - time;
        let mut solution = U::zero();
        {
            let mut solution_iter_mut = solution.iter_mut();
            *solution_iter_mut.next().ok_or("not ok")? = y_0.copy();
            let mut y = y_0;
            let mut y_trial;
            k_4 = function(&time, &y);
            while let Some(next_evaluation_time) = evaluation_time.peek() {
                k_1 = k_4;
                k_2 = function(&(time + 0.5 * timestep), &(&k_1 * (0.5 * timestep) + &y));
                k_3 = function(&(time + 0.75 * timestep), &(&k_2 * (0.75 * timestep) + &y));
                y_trial = (&k_1 * 2.0 + &k_2 * 3.0 + &k_3 * 4.0) * (timestep / 9.0) + &y;
                k_4 = function(&(time + timestep), &y_trial);
                error = (k_1 * -5.0 + k_2 * 6.0 + k_3 * 8.0 + &k_4 * -9.0) * (timestep / 72.0);
                error_norm = error.norm();
                if error_norm < self.abs_tol || error_norm / y_trial.norm() < self.rel_tol {
                    if &time > next_evaluation_time {
                        *solution_iter_mut.next().ok_or("not ok")? = (y_trial.copy() - &y)
                            / timestep
                            * (evaluation_time.next().ok_or("not ok")? - time)
                            + &y;
                    }
                    time += timestep;
                    timestep *= self.inc_fac;
                    y = y_trial;
                } else {
                    timestep *= self.dec_fac;
                }
            }
        }
        Ok(solution)
    }
}
