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
    /// Maximum steps per degree of freedom.
    pub max_steps: usize,
}

impl Default for NewtonRaphson {
    fn default() -> Self {
        Self {
            abs_tol: 1e-1 * ABS_TOL,
            check_minimum: true,
            max_steps: 100,
        }
    }
}

impl<H: Tensor, J: Tensor, X: Tensor> SecondOrder<H, J, X> for NewtonRaphson
where
    J: Div<H, Output = X>,
{
    fn minimize<const D: usize, const U: usize>(
        &self,
        jacobian: impl Fn(&X) -> Result<J, OptimizeError>,
        hessian: impl Fn(&X) -> Result<H, OptimizeError>,
        initial_guess: X,
        dirichlet: Option<Dirichlet<U>>,
        neumann: Option<Neumann>,
    ) -> Result<X, OptimizeError> where [(); D - U]: {
        let max_steps = self.max_steps * initial_guess.iter().count();
        let mut residual;
        let mut solution = initial_guess;
        let mut tangent;
        if let Some(ref bc) = dirichlet {
            bc.places
                .iter()
                .zip(bc.values.iter())
                .for_each(|(place, value)| *solution.get_at_mut(place) = *value);
            
            // need to translate [a][i] indices to [3*a + i] indices
            // but that is specialized!
            // this is getting out of hand?
            let test = solution.copy().eliminate(bc.places);
            // maybe it's finally time for Vec...
            // make a separate set of Tensors as like TensorRank2Vec and TensorRank2Vec2D etc.
            // ALWAYS USE LISTS WHEN POSSIBLE, and size used Vecs whenver possible, but can:
            // (1) get back off nightly
            // (2) do these eliminate-type-things without this giant mess
            // (3) reduce overall templating, maybe compile times,
            // (4) enable using the code without re-compiling (Block params D, E, ...; Yeoh)
            // (5) avoid stack overflow since you managed to get literally everything on the stack,
            // (6) Julia API can finally use fem
            //
            // leave this branch, make a new one off main
            // make some benchmarks! merge em in
            // make new --Vec types for anything that you couldn't use with like Julia
            // benchmark again! make sure nothing is lost!

        }
        for step in 0..max_steps {
            residual = jacobian(&solution)?;
            tangent = hessian(&solution)?;
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
            // need to pass in the types for reduced (J, H) as an Option
            // then how are you going to fill them? by skipping Dirichlet places?
            // or make a math method: fn eliminate(residual<D>, places<U>) -> residual<D - U>

            // how to get rid of the rigid body modes causing low-ass eigenvalues?
            // you might have to remove rows/colums of the Dirichlet BCs
            // but then you no longer know the size unless BC slices become sized
            // and you will further need the generic_const_exprs to subtract
            // at least you only need it for Dirichlet

            // println!("{:?}", (step, residual.norm()));
            // println!("{}", residual);
            // // println!("{}", tangent);
            // println!("{}", residual.copy() / tangent.copy());

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
            max_steps,
            format!("{:?}", &self),
        ))
    }
}
