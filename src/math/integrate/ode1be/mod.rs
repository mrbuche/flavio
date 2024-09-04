#[cfg(test)]
mod test;

use super::super::TensorRank0;
use crate::{ABS_TOL, REL_TOL};

/// Implicit, single-stage, first-order, fixed-step, Runge-Kutta method (the backward Euler method).
pub struct Ode1be {
    /// Absolute error tolerance.
    pub abs_tol: TensorRank0,
    /// Relative error tolerance.
    pub rel_tol: TensorRank0,
}

impl Default for Ode1be {
    fn default() -> Self {
        Self {
            abs_tol: ABS_TOL,
            rel_tol: REL_TOL,
        }
    }
}
