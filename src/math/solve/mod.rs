#[cfg(test)]
mod test;

use crate::{ABS_TOL, REL_TOL, math::{Tensor, Tensors, TensorRank0}};
use std::ops::{Add, Mul};

fn ode_23<const W: usize, T, U>(function: impl Fn(&TensorRank0, &T) -> T, evaluation_times: [TensorRank0; W], y_0: T) -> U
where
    T: Add<T, Output = T> + for <'a> Add<&'a T, Output = T> + Mul<TensorRank0, Output = T> + Tensor

+ std::fmt::Debug
,
    for <'a> &'a T: Mul<TensorRank0, Output = T>,
    U: Iterator + Tensors
{
    // CHANGE UNWRAPS INTO ? ONCE GET THIS TO RETURN A RESULT
    let mut error;
    let mut error_norm;
    let mut evaluation_time = evaluation_times.into_iter().peekable();
    let mut output = U::zeroy().iter_muty();
    // let mut output = U::zeroy().into_iter();
    let mut s_1;
    let mut s_2;
    let mut s_3;
    let mut s_4;
    let mut time = evaluation_time.next().unwrap();
    let mut timestep = evaluation_time.next().unwrap() - time;
    println!("{}", time);
    println!("{}", timestep);
    println!("{}", evaluation_time.peek().unwrap());
    println!("{}", evaluation_time.next().unwrap());
    println!("{}", evaluation_time.peek().unwrap());
    println!("{}", evaluation_time.peek().unwrap());
    panic!();
    let mut y = y_0;
    let mut y_trial;
    while time <= evaluation_times.last().copied().unwrap() {
        s_1 = function(&time, &y);
        s_2 = function(&(time + 0.5 * timestep), &(&s_1 * (0.5 * timestep) + &y));
        s_3 = function(&(time + 0.75 * timestep), &(&s_2 * (0.75 * timestep) + &y));
        y_trial = (&s_1 * 2.0 + &s_2 * 3.0 + &s_3 * 4.0) * (timestep / 9.0) + &y;
        s_4 = function(&(time + timestep), &y_trial);
        error = (s_1 * -5.0 + s_2 * 6.0 + s_3 * 8.0 + s_4 * -9.0) * (timestep / 72.0);
        error_norm = error.normy();
        if error_norm < 1e-5 || error_norm / y_trial.normy() < 1e-5 {
            time += timestep;
            timestep *= 1.2;
            y = y_trial;
            if &time > evaluation_time.peek().unwrap() {
                *output.next().unwrap() = y_trial;
            // maybe need to say that Item of U is T
                // let _: u8 = output.next().unwrap();
            }
        } else {
            timestep *= 0.5;
        }
    }
    output
}
