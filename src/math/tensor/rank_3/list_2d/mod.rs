#[cfg(test)]
pub mod test;

use super::{super::Tensors, list::TensorRank3List, TensorRank0};
use std::array::from_fn;

/// A 2D list of *d*-dimensional tensors of rank 3.
///
/// `D` is the dimension, `I`, `J`, `K` are the configurations `W` and `X` are the list lengths.
pub struct TensorRank3List2D<
    const D: usize,
    const I: usize,
    const J: usize,
    const K: usize,
    const W: usize,
    const X: usize,
>([TensorRank3List<D, I, J, K, W>; X]);

/// Implementation of [`Tensors`] for [`TensorRank3List2D`].
impl<
        const D: usize,
        const I: usize,
        const J: usize,
        const K: usize,
        const W: usize,
        const X: usize,
    > Tensors for TensorRank3List2D<D, I, J, K, W, X>
{
    type Array = [[[[[TensorRank0; D]; D]; D]; W]; X];
    type Item = TensorRank3List<D, I, J, K, W>;
    fn as_array(&self) -> Self::Array {
        let mut array = [[[[[0.0; D]; D]; D]; W]; X];
        array
            .iter_mut()
            .zip(self.iter())
            .for_each(|(entry_rank_4_list, tensor_rank_4_list)| {
                *entry_rank_4_list = tensor_rank_4_list.as_array()
            });
        array
    }
    fn copy(&self) -> Self {
        self.iter().map(|entry| entry.copy()).collect()
    }
    fn iter(&self) -> impl Iterator<Item = &Self::Item> {
        self.0.iter()
    }
    fn iter_mut(&mut self) -> impl Iterator<Item = &mut Self::Item> {
        self.0.iter_mut()
    }
    fn new(array: Self::Array) -> Self {
        array
            .iter()
            .map(|array_i| TensorRank3List::new(*array_i))
            .collect()
    }
    fn zero() -> Self {
        Self(from_fn(|_| TensorRank3List::zero()))
    }
}

impl<
        const D: usize,
        const I: usize,
        const J: usize,
        const K: usize,
        const W: usize,
        const X: usize,
    > FromIterator<TensorRank3List<D, I, J, K, W>> for TensorRank3List2D<D, I, J, K, W, X>
{
    fn from_iter<Ii: IntoIterator<Item = TensorRank3List<D, I, J, K, W>>>(
        into_iterator: Ii,
    ) -> Self {
        let mut tensor_rank_3_list_2d = Self::zero();
        tensor_rank_3_list_2d
            .iter_mut()
            .zip(into_iterator)
            .for_each(|(tensor_rank_3_list, entry)| *tensor_rank_3_list = entry);
        tensor_rank_3_list_2d
    }
}
