mod newton;

use super::Tensor;
use std::{fmt, ops::Div};

// Newton doesn't need the objective function,
// so whether it's actually an optimization or root finding problem,
// the i/o is the same, so it can live in optimize safely.

/// ???
pub trait Optimize<X, J, H>
where
    Self: fmt::Debug,
    X: Tensor,
    J: Tensor + Div<H, Output = X>,
    H: Tensor,
{
    fn minimize(
        &self,
        jacobian: impl Fn(&X) -> J,
        hessian: impl Fn(&X) -> H,
        initial_guess: X,
    ) -> Result<X, OptimizeError>;
}

/// ???
pub enum OptimizeError {
    MaximumStepsReached(usize, String),
}
