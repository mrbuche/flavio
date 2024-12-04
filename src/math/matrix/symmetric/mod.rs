#[cfg(test)]
use crate::math::test::ErrorTensor;

use crate::math::{Hessian, Rank2, Tensor, TensorRank0, TensorVec, Vector, write_tensor_rank_0, tensor::TensorError};
use std::{
    fmt,
    ops::{Add, AddAssign, Div, DivAssign, Index, IndexMut, Mul, MulAssign, Sub, SubAssign},
};

/// A symmetric matrix.
#[derive(Debug)]
pub struct SymmetricMatrix(Vec<Vector>);

//
// the packed (1D) storage is allegedly bad for vectorization
// so why not store a Vec of Vectors
// and put None at the end of each where they would go into the other half?
// or something like that
// but then you get into weird issues when trying to vectorize things like matrix multiplication, right?

//
// you can always benchmark against SquareMatrix later
// if you can't fill faster, and can't computation faster, then need other algorithms
// storage being the other factor
//

//
// maybe just worry about trying to match the speed of SquareMatrix for now
// while also testing that the storage size is correct (N*(N+1)/2 instead of N*N)

//
// use Cholesky for inverse and things?
// but being symmetric doesn't mean it's positive-semidefinite
// so inverse should have an option to do Cholesky or not
// can apply the same option to SquareMatrix too
//

//
// can also have ode solvers use tspan like matlab/etc. now
//

//
// remember: trying to get something more pliable (TensorVec, Matrix, etc.) to go back to Lagrange multipliers instead of constraints
//           and might want to try to get rid of get_at methods if possible?
//

#[cfg(test)]
impl ErrorTensor for SymmetricMatrix {
    fn error(
        &self,
        comparator: &Self,
        tol_abs: &TensorRank0,
        tol_rel: &TensorRank0,
    ) -> Option<usize> {
        todo!()
    }
    fn error_fd(&self, comparator: &Self, epsilon: &TensorRank0) -> Option<(bool, usize)> {
        todo!()
    }
}

impl fmt::Display for SymmetricMatrix {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        todo!()
    }
}

impl PartialEq for SymmetricMatrix {
    fn eq(&self, other: &Self) -> bool {
        todo!()
    }
}

impl FromIterator<Vector> for SymmetricMatrix {
    fn from_iter<Ii: IntoIterator<Item = Vector>>(into_iterator: Ii) -> Self {
        todo!();
        Self(Vec::from_iter(into_iterator))
    }
}

impl Index<usize> for SymmetricMatrix {
    type Output = Vector;
    fn index(&self, index: usize) -> &Self::Output {
        todo!();
        &self.0[index]
    }
}

impl IndexMut<usize> for SymmetricMatrix {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        todo!();
        &mut self.0[index]
    }
}

impl Hessian for SymmetricMatrix {
    fn is_positive_definite(&self) -> bool {
        self.cholesky_decomposition().is_ok()
    }
}

impl Rank2 for SymmetricMatrix {
    fn cholesky_decomposition(&self) -> Result<SymmetricMatrix, TensorError> {
        todo!()
    }
}

impl Tensor for SymmetricMatrix {
    type Item = Vector;
    fn copy(&self) -> Self {
        todo!()
    }
    fn get_at(&self, indices: &[usize]) -> &TensorRank0 {
        todo!()
    }
    fn get_at_mut(&mut self, indices: &[usize]) -> &mut TensorRank0 {
        todo!()
    }
    fn iter(&self) -> impl Iterator<Item = &Self::Item> {
        todo!();
        self.0.iter()
    }
    fn iter_mut(&mut self) -> impl Iterator<Item = &mut Self::Item> {
        todo!();
        self.0.iter_mut()
    }
}

impl<'a> TensorVec<'a> for SymmetricMatrix {
    type Item = Vector;
    type Slice = &'a [&'a [TensorRank0]];
    fn is_empty(&self) -> bool {
        todo!()
    }
    fn len(&self) -> usize {
        todo!()
    }
    fn new(slice: Self::Slice) -> Self {
        todo!()
    }
    fn zero(len: usize) -> Self {
        todo!()
    }
}

impl Div<TensorRank0> for SymmetricMatrix {
    type Output = Self;
    fn div(mut self, tensor_rank_0: TensorRank0) -> Self::Output {
        self /= &tensor_rank_0;
        self
    }
}

impl Div<&TensorRank0> for SymmetricMatrix {
    type Output = Self;
    fn div(mut self, tensor_rank_0: &TensorRank0) -> Self::Output {
        self /= tensor_rank_0;
        self
    }
}

impl DivAssign<TensorRank0> for SymmetricMatrix {
    fn div_assign(&mut self, tensor_rank_0: TensorRank0) {
        todo!()
    }
}

impl DivAssign<&TensorRank0> for SymmetricMatrix {
    fn div_assign(&mut self, tensor_rank_0: &TensorRank0) {
        todo!()
    }
}

impl Mul<TensorRank0> for SymmetricMatrix {
    type Output = Self;
    fn mul(mut self, tensor_rank_0: TensorRank0) -> Self::Output {
        self *= &tensor_rank_0;
        self
    }
}
impl Mul<&TensorRank0> for SymmetricMatrix {
    type Output = Self;
    fn mul(mut self, tensor_rank_0: &TensorRank0) -> Self::Output {
        self *= tensor_rank_0;
        self
    }
}

impl Mul<&TensorRank0> for &SymmetricMatrix {
    type Output = SymmetricMatrix;
    fn mul(self, tensor_rank_0: &TensorRank0) -> Self::Output {
        todo!()
    }
}

impl MulAssign<TensorRank0> for SymmetricMatrix {
    fn mul_assign(&mut self, tensor_rank_0: TensorRank0) {
        todo!()
    }
}

impl MulAssign<&TensorRank0> for SymmetricMatrix {
    fn mul_assign(&mut self, tensor_rank_0: &TensorRank0) {
        todo!()
    }
}

impl Add for SymmetricMatrix {
    type Output = Self;
    fn add(mut self, vector: Self) -> Self::Output {
        self += vector;
        self
    }
}

impl Add<&Self> for SymmetricMatrix {
    type Output = Self;
    fn add(mut self, vector: &Self) -> Self::Output {
        self += vector;
        self
    }
}

impl AddAssign for SymmetricMatrix {
    fn add_assign(&mut self, vector: Self) {
        todo!()
    }
}

impl AddAssign<&Self> for SymmetricMatrix {
    fn add_assign(&mut self, vector: &Self) {
        todo!()
    }
}

impl Sub for SymmetricMatrix {
    type Output = Self;
    fn sub(mut self, vector: Self) -> Self::Output {
        self -= vector;
        self
    }
}

impl Sub<&Self> for SymmetricMatrix {
    type Output = Self;
    fn sub(mut self, vector: &Self) -> Self::Output {
        self -= vector;
        self
    }
}

impl SubAssign for SymmetricMatrix {
    fn sub_assign(&mut self, vector: Self) {
        todo!()
    }
}

impl SubAssign<&Self> for SymmetricMatrix {
    fn sub_assign(&mut self, vector: &Self) {
        todo!()
    }
}

