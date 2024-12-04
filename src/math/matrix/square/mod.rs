use crate::math::{Tensor, TensorRank0, TensorVec, Vector};

// do Vector, Matrix, and MatrixSym as the 1D and 2D Vec types
// MatrixSym so you can do Cholesky and store ~1/2 the data
// try to keep MatrixSym from impl less-efficient things that don't use things like Cholesky

// can also have ode solvers use tspan like matlab/etc. now

// remember: trying to get something more pliable (TensorVec, Matrix, etc.) to go back to Lagrange multipliers instead of constraints
//           and might want to try to get rid of get_at methods if possible?

/// A square matrix.
pub struct SquareMatrix(Vec<Vector>);

impl FromIterator<Vector> for SquareMatrix {
    fn from_iter<Ii: IntoIterator<Item = Vector>>(into_iterator: Ii) -> Self {
        Self(Vec::from_iter(into_iterator))
    }
}

// impl Tensor for SquareMatrix {
//     type Item = Vec<TensorRank0>;
//     fn copy(&self) -> Self {
//         self.iter().map(|entry| entry.copy()).collect()
//     }
//     fn get_at(&self, indices: &[usize]) -> &TensorRank0 {
//         &self[indices[0]][indices[1]]
//     }
//     fn get_at_mut(&mut self, indices: &[usize]) -> &mut TensorRank0 {
//         &mut self[indices[0]][indices[1]]
//     }
//     fn iter(&self) -> impl Iterator<Item = &Self::Item> {
//         self.0.iter()
//     }
//     fn iter_mut(&mut self) -> impl Iterator<Item = &mut Self::Item> {
//         self.0.iter_mut()
//     }
// }

impl<'a> TensorVec<'a> for SquareMatrix {
    type Item = Vector;
    type Slice = &'a [&'a [TensorRank0]];
    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
    fn len(&self) -> usize {
        self.0.len()
    }
    fn new(slice: Self::Slice) -> Self {
        slice
            .iter()
            .map(|slice_entry| Self::Item::new(slice_entry))
            .collect()
    }
    fn zero(len: usize) -> Self {
        (0..len).map(|_| Self::Item::zero(len)).collect()
    }
}
