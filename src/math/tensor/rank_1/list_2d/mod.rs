#[cfg(test)]
mod test;

use super::{super::Tensors, list::TensorRank1List, TensorRank0};
use std::array::from_fn;

/// A 2D list of *d*-dimensional tensors of rank 1.
///
/// `D` is the dimension, `I` is the configuration, `W` and `X` are the list lengths.
pub struct TensorRank1List2D<const D: usize, const I: usize, const W: usize, const X: usize>(
    pub [TensorRank1List<D, I, W>; X],
);

/// Implementation of [`Tensors`] for [`TensorRank1List2D`].
impl<const D: usize, const I: usize, const W: usize, const X: usize> Tensors
    for TensorRank1List2D<D, I, W, X>
{
    type Array = [[[TensorRank0; D]; W]; X];
    type Item = TensorRank1List<D, I, W>;
    fn as_array(&self) -> Self::Array {
        let mut array = [[[0.0; D]; W]; X];
        array
            .iter_mut()
            .zip(self.iter())
            .for_each(|(entry, tensor_rank_1_list)| *entry = tensor_rank_1_list.as_array());
        array
    }
    fn iter(&self) -> impl Iterator<Item = &TensorRank1List<D, I, W>> {
        self.0.iter()
    }
    fn iter_mut(&mut self) -> impl Iterator<Item = &mut Self::Item> {
        self.0.iter_mut()
    }
    fn new(array: Self::Array) -> Self {
        array
            .iter()
            .map(|array_i| TensorRank1List::new(*array_i))
            .collect()
    }
    fn zero() -> Self {
        Self(from_fn(|_| TensorRank1List::zero()))
    }
}

impl<const D: usize, const I: usize, const W: usize, const X: usize>
    FromIterator<TensorRank1List<D, I, W>> for TensorRank1List2D<D, I, W, X>
{
    fn from_iter<Ii: IntoIterator<Item = TensorRank1List<D, I, W>>>(into_iterator: Ii) -> Self {
        let mut tensor_rank_1_list_2d = Self::zero();
        tensor_rank_1_list_2d
            .iter_mut()
            .zip(into_iterator)
            .for_each(|(tensor_rank_1_list, entry)| *tensor_rank_1_list = entry);
        tensor_rank_1_list_2d
    }
}
