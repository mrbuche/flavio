use crate::math::{Tensor, TensorRank0};

// do Vector, Matrix, and MatrixSym as the 1D and 2D Vec types
// MatrixSym so you can do Cholesky and store ~1/2 the data
// try to keep MatrixSym from impl less-efficient things that don't use things like Cholesky

// can also have ode solvers use tspan like matlab/etc. now
// might want to do that now since they use zero(), i.e. turn U into Type(Vec<Y>)

// remember: trying to get something more pliable (TensorVec, Matrix, etc.) to go back to Lagrange multipliers instead of constraints

/// ???
pub struct SquareMatrix(Vec<Vec<TensorRank0>>);

impl Tensor for SquareMatrix {
    type Item = Vec<TensorRank0>;
    fn copy(&self) -> Self {
        self.iter().map(|entry| entry.copy()).collect()
    }
    fn get_at(&self, indices: &[usize]) -> &TensorRank0 {
        &self[indices[0]][indices[1]]
    }
    fn get_at_mut(&mut self, indices: &[usize]) -> &mut TensorRank0 {
        &mut self[indices[0]][indices[1]]
    }
    fn is_positive_definite(&self) -> bool {
        todo!()
    }
    fn iter(&self) -> impl Iterator<Item = &Self::Item> {
        self.0.iter()
    }
    fn iter_mut(&mut self) -> impl Iterator<Item = &mut Self::Item> {
        self.0.iter_mut()
    }
}