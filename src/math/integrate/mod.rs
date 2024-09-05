#[cfg(test)]
mod test;

mod ode1be;
mod ode23;
mod ode45;

pub use ode1be::Ode1be;
pub use ode23::Ode23;
pub use ode45::Ode45;

use super::{Tensor, TensorRank0, TensorRank0List, Tensors};
use crate::get_defeat_message;
use std::fmt;

/// Base trait for explicit ordinary different equation solvers.
pub trait Explicit {
    /// Solves an initial value problem by explicitly integrating a system of ordinary different equations.
    ///
    /// ```math
    /// \frac{dy}{dt} = f(t, y),\quad y(t_0) = y_0
    /// ```
    fn integrate<Y, U, const W: usize>(
        &self,
        function: impl Fn(&TensorRank0, &Y) -> Y,
        initial_time: TensorRank0,
        initial_condition: Y,
        evaluation_times: &TensorRank0List<W>,
    ) -> Result<U, IntegrationError<W>>
    where
        Y: Tensor,
        for<'a> &'a Y: std::ops::Mul<TensorRank0, Output = Y>,
        U: Tensors<Item = Y>;
}

/// Base trait for implicit ordinary different equation solvers.
pub trait Implicit {
    /// Solves an initial value problem by implicitly integrating a system of ordinary different equations.
    ///
    /// ```math
    /// \frac{dy}{dt} = f(t, y),\quad y(t_0) = y_0,\quad \frac{\partial f}{\partial y} = J(t, y)
    /// ```
    fn integrate<Y, J, U, const W: usize>(
        &self,
        function: impl Fn(&TensorRank0, &Y) -> Y,
        jacobian: impl Fn(&TensorRank0, &Y) -> J,
        initial_time: TensorRank0,
        initial_condition: Y,
        evaluation_times: &TensorRank0List<W>,
    ) -> Result<U, IntegrationError<W>>
    where
        Y: Tensor,
        for<'a> &'a Y: std::ops::Mul<TensorRank0, Output = Y>,
        U: Tensors<Item = Y>;
}

/// Possible errors encountered when integrating.
pub enum IntegrationError<const W: usize> {
    EvaluationTimesNoFinalTime(TensorRank0List<W>, String),
    EvaluationTimesNotStrictlyIncreasing(TensorRank0List<W>, String),
    EvaluationTimesPreceedInitialTime(TensorRank0List<W>, TensorRank0, String),
}

impl<const W: usize> From<&str> for IntegrationError<W> {
    fn from(string: &str) -> Self {
        todo!("{}", string)
    }
}

impl<const W: usize> fmt::Debug for IntegrationError<W> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let error = match self {
            Self::EvaluationTimesNoFinalTime(evaluation_times, integrator) => {
                format!(
                    "\x1b[1;91mEvaluation times must include a final time.\x1b[0;91m\n\
                     From evaluation times: {}.\n\
                     In integrator: {}.",
                    evaluation_times, integrator
                )
            }
            Self::EvaluationTimesNotStrictlyIncreasing(evaluation_times, integrator) => {
                format!(
                    "\x1b[1;91mEvaluation times must be strictly increasing.\x1b[0;91m\n\
                     From evaluation times: {}.\n\
                     In integrator: {}.",
                    evaluation_times, integrator
                )
            }
            Self::EvaluationTimesPreceedInitialTime(evaluation_times, initial_time, integrator) => {
                format!(
                    "\x1b[1;91mEvaluation times preceed the initial time.\x1b[0;91m\n\
                     From evaluation times: {}.\n\
                     With initial time: {}.\n\
                     In integrator: {}.",
                    evaluation_times, initial_time, integrator
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

impl<const W: usize> fmt::Display for IntegrationError<W> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let error = match self {
            Self::EvaluationTimesNoFinalTime(evaluation_times, integrator) => {
                format!(
                    "\x1b[1;91mEvaluation times must include a final time.\x1b[0;91m\n\
                     From evaluation times: {}.\n\
                     In integrator: {}.",
                    evaluation_times, integrator
                )
            }
            Self::EvaluationTimesNotStrictlyIncreasing(evaluation_times, integrator) => {
                format!(
                    "\x1b[1;91mEvaluation times must be strictly increasing.\x1b[0;91m\n\
                     From evaluation times: {}.\n\
                     In integrator: {}.",
                    evaluation_times, integrator
                )
            }
            Self::EvaluationTimesPreceedInitialTime(evaluation_times, initial_time, integrator) => {
                format!(
                    "\x1b[1;91mEvaluation times preceed the initial time.\x1b[0;91m\n\
                     From evaluation times: {}.\n\
                     With initial time: {}.\n\
                     In integrator: {}.",
                    evaluation_times, initial_time, integrator
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
