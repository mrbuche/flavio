#[cfg(test)]
mod test;

use super::{
    super::{Tensor, TensorRank0, TensorRank0List, Tensors},
    SolveError,
};
use crate::{ABS_TOL, REL_TOL};

fn ode23<const W: usize, T, U>(
    function: impl Fn(&TensorRank0, &T) -> T,
    evaluation_times: &TensorRank0List<W>,
    y_0: T,
) -> Result<U, SolveError>
where
    T: Tensor,
    for<'a> &'a T: std::ops::Mul<TensorRank0, Output = T>,
    U: Tensors<Item = T>,
{
    let mut error;
    let mut error_norm;
    let mut evaluation_time = evaluation_times.0.into_iter().peekable();
    let mut s_1;
    let mut s_2;
    let mut s_3;
    let mut s_4;
    let mut time = evaluation_time.next().ok_or("not ok")?;
    let mut timestep = evaluation_time.peek().ok_or("not ok")? - time;
    let mut solution = U::zero();
    {
        let mut solution_iter_mut = solution.iter_mut();
        *solution_iter_mut.next().ok_or("not ok")? = y_0.copy();
        let mut y = y_0;
        let mut y_trial;
        while let Some(next_evaluation_time) = evaluation_time.peek() {
            s_1 = function(&time, &y);
            s_2 = function(&(time + 0.5 * timestep), &(&s_1 * (0.5 * timestep) + &y));
            s_3 = function(&(time + 0.75 * timestep), &(&s_2 * (0.75 * timestep) + &y));
            y_trial = (&s_1 * 2.0 + &s_2 * 3.0 + &s_3 * 4.0) * (timestep / 9.0) + &y;
            s_4 = function(&(time + timestep), &y_trial);
            error = (s_1 * -5.0 + s_2 * 6.0 + s_3 * 8.0 + s_4 * -9.0) * (timestep / 72.0);
            error_norm = error.norm();
            if error_norm < ABS_TOL || error_norm / y_trial.norm() < REL_TOL {
                if &time > next_evaluation_time {
                    *solution_iter_mut.next().ok_or("not ok")? = (y_trial.copy() - &y) / timestep
                        * (evaluation_time.next().ok_or("not ok")? - time)
                        + &y;
                }
                time += timestep;
                timestep *= 1.1;
                y = y_trial;
            } else {
                timestep *= 0.5;
            }
        }
    }
    Ok(solution)
}
