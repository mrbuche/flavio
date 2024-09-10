#[cfg(test)]
mod test;

pub mod list;

use super::Tensor;

/// A tensor of rank 0 (a scalar).
pub type TensorRank0 = f64;

/// Implementation of [`Tensor`] for [`TensorRank0`].
impl Tensor for TensorRank0 {
    type Array = [Self; 1];
    type Item = TensorRank0;
    fn as_array(&self) -> Self::Array {
        [self.copy()]
    }
    fn copy(&self) -> TensorRank0 {
        *self
    }
    fn identity() -> Self {
        1.0
    }
    fn is_zero(&self) -> bool {
        self == &0.0
    }
    fn iter(&self) -> impl Iterator<Item = &Self::Item> {
        [0.0].iter()
    }
    fn iter_mut(&mut self) -> impl Iterator<Item = &mut Self::Item> {
        [self].into_iter()
    }
    fn new(array: Self::Array) -> Self {
        array[0]
    }
    fn norm_squared(&self) -> TensorRank0 {
        self.powi(2)
    }
    fn normalized(&self) -> Self {
        1.0
    }
    fn zero() -> Self {
        0.0
    }
}
