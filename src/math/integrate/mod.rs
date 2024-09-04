#[cfg(test)]
mod test;

mod ode1be;
mod ode23;
mod ode45;

pub use ode45::ode45;

pub use ode1be::Ode1be;
pub use ode23::Ode23;

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
        t_0: TensorRank0,
        y_0: Y,
        evaluation_times: &TensorRank0List<W>,
    ) -> Result<U, IntegrationError>
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
        t_0: TensorRank0,
        y_0: Y,
        jac: J,
        evaluation_times: &TensorRank0List<W>,
    ) -> Result<U, IntegrationError>
    where
        Y: Tensor,
        for<'a> &'a Y: std::ops::Mul<TensorRank0, Output = Y>,
        U: Tensors<Item = Y>;
}

/// Possible errors encountered when integrating.
pub enum IntegrationError {
    GeneralError,
}

impl fmt::Debug for IntegrationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let error = match self {
            Self::GeneralError => "\x1b[1;91m???.\x1b[0;91m".to_string(),
        };
        write!(
            f,
            "\n{}\n\x1b[0;2;31m{}\x1b[0m\n",
            error,
            get_defeat_message()
        )
    }
}

impl fmt::Display for IntegrationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let error = match self {
            Self::GeneralError => "\x1b[1;91m???.\x1b[0;91m".to_string(),
        };
        write!(
            f,
            "\n{}\n\x1b[0;2;31m{}\x1b[0m\n",
            error,
            get_defeat_message()
        )
    }
}

impl From<&str> for IntegrationError {
    fn from(string: &str) -> Self {
        todo!("{}", string)
    }
}
