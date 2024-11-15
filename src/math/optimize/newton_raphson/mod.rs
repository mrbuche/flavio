#[cfg(test)]
mod test;

use super::{
    super::{Tensor, TensorRank0},
    Dirichlet, Neumann, OptimizeError, SecondOrder,
};
use crate::ABS_TOL;
use std::ops::Div;

/// The Newton-Raphson method.
#[derive(Debug)]
pub struct NewtonRaphson {
    /// Absolute error tolerance.
    pub abs_tol: TensorRank0,
    /// Whether to check if solution is minimum.
    pub check_minimum: bool,
    /// Maximum number of steps.
    pub max_steps: usize,
}

impl Default for NewtonRaphson {
    fn default() -> Self {
        Self {
            abs_tol: ABS_TOL,
            check_minimum: true,
            max_steps: 100,
        }
    }
}

impl<H: Tensor, J: Tensor, X: Tensor> SecondOrder<H, J, X> for NewtonRaphson
where
    J: Div<H, Output = X>,
{
    fn minimize(
        &self,
        jacobian: impl Fn(&X) -> Result<J, OptimizeError>,
        hessian: impl Fn(&X) -> Result<H, OptimizeError>,
        initial_guess: X,
        dirichlet: Option<Dirichlet>,
        neumann: Option<Neumann>,
    ) -> Result<X, OptimizeError> {
        let mut residual;
        let mut solution = initial_guess;
        if let Some(ref bc) = dirichlet {
            bc.places
                .iter()
                .zip(bc.values.iter())
                .for_each(|(place, value)| *solution.get_at_mut(place) = *value)
        }
        let mut tangent;
        for step in 0..self.max_steps {
            residual = jacobian(&solution)?;
            if let Some(ref bc) = neumann {
                bc.places
                    .iter()
                    .zip(bc.values.iter())
                    .for_each(|(place, value)| *residual.get_at_mut(place) -= value)
            }
            if let Some(ref bc) = dirichlet {
                bc.places
                    .iter()
                    .for_each(|place| *residual.get_at_mut(place) = 0.0)
            }
            tangent = hessian(&solution)?;

            // how to get rid of the rigid body modes causing low-ass eigenvalues?
            // you might have to remove rows/colums of the Dirichlet BCs
            // but then you no longer know the size unless BC slices become sized
            // and you will further need the generic_const_exprs to subtract
            println!("{:?}", (step, residual.norm()));
            println!("{}", residual);
            // println!("{}", tangent);
            println!("{}", residual.copy() / tangent.copy());

            if residual.norm() < self.abs_tol {
                if self.check_minimum && !tangent.is_positive_definite() {
                    return Err(OptimizeError::NotMinimum(
                        format!("{}", solution),
                        format!("{:?}", &self),
                    ));
                } else {
                    return Ok(solution);
                }
            } else {
                solution -= residual / tangent;
            }
        }
        Err(OptimizeError::MaximumStepsReached(
            self.max_steps,
            format!("{:?}", &self),
        ))
    }
}
