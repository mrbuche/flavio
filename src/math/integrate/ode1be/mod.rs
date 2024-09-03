#[cfg(test)]
mod test;

use super::{
    super::{Tensor, TensorRank0, TensorRank0List, Tensors},
    IntegrationError,
};
// use crate::{ABS_TOL, REL_TOL};

/// Implicit first-order Runge-Kutta method.
///
/// [`ode1be`] is an implicit, single-stage, first-order, fixed-step Runge-Kutta method (the backward Euler method).
///
pub fn ode1be<const W: usize, T, U>(
    _function: impl Fn(&TensorRank0, &T) -> T,
    _evaluation_times: &TensorRank0List<W>,
    _y_0: T,
) -> Result<U, IntegrationError>
where
    T: Tensor,
    for<'a> &'a T: std::ops::Mul<TensorRank0, Output = T>,
    U: Tensors<Item = T>,
{
    todo!()
}
