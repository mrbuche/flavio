#[cfg(test)]
mod test;

use std::{
    array::from_fn,
    ops::{Index, IndexMut},
};

use super::{super::Tensors, Tensor, TensorRank0, TensorRank4};

/// A list of *d*-dimensional tensor of rank 4.
///
/// `D` is the dimension, `I`, `J`, `K`, `L` are the configurations, `W` is the list length.
pub struct TensorRank4List<
    const D: usize,
    const I: usize,
    const J: usize,
    const K: usize,
    const L: usize,
    const W: usize,
>([TensorRank4<D, I, J, K, L>; W]);

/// Implementation of [`Tensors`] for [`TensorRank4List`].
impl<
        const D: usize,
        const I: usize,
        const J: usize,
        const K: usize,
        const L: usize,
        const W: usize,
    > Tensors for TensorRank4List<D, I, J, K, L, W>
{
    type Array = [[[[[TensorRank0; D]; D]; D]; D]; W];
    type Item = TensorRank4<D, I, J, K, L>;
    fn as_array(&self) -> Self::Array {
        let mut array = [[[[[0.0; D]; D]; D]; D]; W];
        array
            .iter_mut()
            .zip(self.iter())
            .for_each(|(entry_rank_4, tensor_rank_4)| *entry_rank_4 = tensor_rank_4.as_array());
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
            .map(|array_i| TensorRank4::new(*array_i))
            .collect()
    }
    fn zero() -> Self {
        Self(from_fn(|_| TensorRank4::zero()))
    }
}

impl<
        const D: usize,
        const I: usize,
        const J: usize,
        const K: usize,
        const L: usize,
        const W: usize,
    > FromIterator<TensorRank4<D, I, J, K, L>> for TensorRank4List<D, I, J, K, L, W>
{
    fn from_iter<Ii: IntoIterator<Item = TensorRank4<D, I, J, K, L>>>(into_iterator: Ii) -> Self {
        let mut tensor_rank_4_list = Self::zero();
        tensor_rank_4_list
            .iter_mut()
            .zip(into_iterator)
            .for_each(|(tensor_rank_4_list_entry, entry)| *tensor_rank_4_list_entry = entry);
        tensor_rank_4_list
    }
}

impl<
        const D: usize,
        const I: usize,
        const J: usize,
        const K: usize,
        const L: usize,
        const W: usize,
    > Index<usize> for TensorRank4List<D, I, J, K, L, W>
{
    type Output = TensorRank4<D, I, J, K, L>;
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl<
        const D: usize,
        const I: usize,
        const J: usize,
        const K: usize,
        const L: usize,
        const W: usize,
    > IndexMut<usize> for TensorRank4List<D, I, J, K, L, W>
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}
