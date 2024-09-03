#[cfg(test)]
mod test;

use super::{
    super::{Tensor, TensorRank0, TensorRank0List, Tensors},
    IntegrationError,
};
// use crate::{ABS_TOL, REL_TOL};

/// Explicit fifth-order Runge-Kutta method.
///
/// [`ode45`] is an adaptive, explicit, six-stage, fifth-order Runge-Kutta method ([Dormand and Prince, 1980](https://doi.org/10.1016/0771-050X(80)90013-3)).
pub fn ode45<const W: usize, T, U>(
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
