#[cfg(test)]
mod test;

pub mod list;

use super::Tensor;

/// A tensor of rank 0 (a scalar).
pub type TensorRank0 = f64;

/// Implementation of [`Tensor`] for [`TensorRank0`].
impl Tensor for TensorRank0 {
    type Array = [Self; 1];
    fn as_array(&self) -> Self::Array {
        [self.copy()]
    }
    fn copy(&self) -> TensorRank0 {
        *self
    }
    fn new(array: Self::Array) -> Self {
        array[0]
    }
    fn norm(&self) -> TensorRank0 {
        self.abs()
    }
    fn zero() -> Self {
        0.0
    }
}
