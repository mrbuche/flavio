#[cfg(test)]
pub mod test;

use std::array::from_fn;
use super::{
    list_2d::TensorRank3List2D,
    TensorRank0, super::Tensors
};

type MakeClippyHappy<const D: usize> = [[[TensorRank0; D]; D]; D];

/// A 3D list of *d*-dimensional tensors of rank 3.
///
/// `D` is the dimension, `I`, `J`, `K` are the configurations `W`, `X`, and `Y` are the list lengths.
pub struct TensorRank3List3D<
    const D: usize,
    const I: usize,
    const J: usize,
    const K: usize,
    const W: usize,
    const X: usize,
    const Y: usize,
>([TensorRank3List2D<D, I, J, K, W, X>; Y]);

/// Implementation of [`Tensors`] for [`TensorRank3List3D`].
impl<const D: usize, const I: usize, const J: usize, const K: usize, const W: usize, const X: usize, const Y: usize> Tensors for TensorRank3List3D<D, I, J, K, W, X, Y> {
    type Array = [[[MakeClippyHappy<D>; W]; X]; Y];
    type Item = TensorRank3List2D<D, I, J, K, W, X>;
    fn as_array(&self) -> Self::Array {
        let mut array = [[[[[[0.0; D]; D]; D]; W]; X]; Y];
        array
            .iter_mut()
            .zip(self.iter())
            .for_each(|(entry_rank_4_list_2d, tensor_rank_4_list_2d)| *entry_rank_4_list_2d = tensor_rank_4_list_2d.as_array());
        array
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
            .map(|array_i| TensorRank3List2D::new(*array_i))
            .collect()
    }
    fn zero() -> Self {
        Self(from_fn(|_| TensorRank3List2D::zero()))
    }
}

impl<
        const D: usize,
        const I: usize,
        const J: usize,
        const K: usize,
        const W: usize,
        const X: usize,
        const Y: usize,
    > FromIterator<TensorRank3List2D<D, I, J, K, W, X>> for TensorRank3List3D<D, I, J, K, W, X, Y>
{
    fn from_iter<Ii: IntoIterator<Item = TensorRank3List2D<D, I, J, K, W, X>>>(
        into_iterator: Ii,
    ) -> Self {
        let mut tensor_rank_3_list_3d = Self::zero();
        tensor_rank_3_list_3d
            .iter_mut()
            .zip(into_iterator)
            .for_each(|(tensor_rank_3_list_2d, entry)| *tensor_rank_3_list_2d = entry);
        tensor_rank_3_list_3d
    }
}
