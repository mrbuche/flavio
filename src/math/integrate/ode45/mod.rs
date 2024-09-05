#[cfg(test)]
mod test;

use super::super::TensorRank0;
use crate::{ABS_TOL, REL_TOL};

/// Explicit, six-stage, fifth-order, variable-step, Runge-Kutta method ([Dormand and Prince, 1980](https://doi.org/10.1016/0771-050X(80)90013-3)).
pub struct Ode45 {
    /// Absolute error tolerance.
    pub abs_tol: TensorRank0,
    /// Multiplying factor when decreasing time steps.
    pub dec_fac: TensorRank0,
    /// Multiplying factor when increasing time steps.
    pub inc_fac: TensorRank0,
    /// Relative error tolerance.
    pub rel_tol: TensorRank0,
}

impl Default for Ode45 {
    fn default() -> Self {
        Self {
            abs_tol: ABS_TOL,
            dec_fac: 0.5,
            inc_fac: 1.1,
            rel_tol: REL_TOL,
        }
    }
}
