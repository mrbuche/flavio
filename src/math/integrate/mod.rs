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

/// ???
pub enum Ivp<F: Fn(&TensorRank0, &Y) -> Y, Y> {
    A(F, Y), // give fun, y0 (SHOULD GIVE t0 too, THAT IS PART OF THE IVP STATEMENT)
    B(F, Y), // give that plus Jacobian: if no Jacobian given in Implicit methods, use finite difference?
}

/// Base trait for ordinary different equation solvers.
pub trait OdeSolver {
    /// Solves an initial value problem by numerically integrating a system of ordinary different equations.
    ///
    /// ```math
    /// \frac{dy}{dt} = f(t, y),\quad y(t_0) = y_0
    /// ```
    fn integrate<F: Fn(&TensorRank0, &T) -> T, T, U, const W: usize>(
        &self,
        ivp: Ivp<F, T>,
        // function: impl Fn(&TensorRank0, &T) -> T,
        evaluation_times: &TensorRank0List<W>,
    ) -> Result<U, IntegrationError>
    where
        T: Tensor,
        for<'a> &'a T: std::ops::Mul<TensorRank0, Output = T>,
        U: Tensors<Item = T>;
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
