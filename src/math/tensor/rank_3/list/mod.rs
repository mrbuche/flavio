#[cfg(test)]
pub mod test;

use std::{array::from_fn, ops::AddAssign};

use super::{TensorRank0, TensorRank3, Tensor, super::Tensors};

/// A list of *d*-dimensional tensors of rank 3.
///
/// `D` is the dimension, `I`, `J`, `K` are the configurations `W` is the list length.
pub struct TensorRank3List<
    const D: usize,
    const I: usize,
    const J: usize,
    const K: usize,
    const W: usize,
>([TensorRank3<D, I, J, K>; W]);

/// Implementation of [`Tensors`] for [`TensorRank3List`].
impl<const D: usize, const I: usize, const J: usize, const K: usize, const W: usize> Tensors for TensorRank3List<D, I, J, K, W> {
    type Array = [[[[TensorRank0; D]; D]; D]; W];
    type Item = TensorRank3<D, I, J, K>;
    fn as_array(&self) -> Self::Array {
        let mut array = [[[[0.0; D]; D]; D]; W];
        array
            .iter_mut()
            .zip(self.iter())
            .for_each(|(entry_rank_3, tensor_rank_3)| *entry_rank_3 = tensor_rank_3.as_array());
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
            .map(|array_i| TensorRank3::new(*array_i))
            .collect()
    }
    fn zero() -> Self {
        Self(from_fn(|_| TensorRank3::zero()))
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize, const W: usize>
    FromIterator<TensorRank3<D, I, J, K>> for TensorRank3List<D, I, J, K, W>
{
    fn from_iter<Ii: IntoIterator<Item = TensorRank3<D, I, J, K>>>(into_iterator: Ii) -> Self {
        let mut tensor_rank_3_list = Self::zero();
        tensor_rank_3_list
            .iter_mut()
            .zip(into_iterator)
            .for_each(|(tensor_rank_3_list_entry, entry)| *tensor_rank_3_list_entry = entry);
        tensor_rank_3_list
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize, const W: usize> std::iter::Sum
    for TensorRank3List<D, I, J, K, W>
{
    fn sum<Ii>(iter: Ii) -> Self
    where
        Ii: Iterator<Item = Self>,
    {
        let mut output = Self::zero();
        iter.for_each(|item| output += item);
        output
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize, const W: usize> AddAssign
    for TensorRank3List<D, I, J, K, W>
{
    fn add_assign(&mut self, tensor_rank_3_list: Self) {
        self.iter_mut()
            .zip(tensor_rank_3_list.iter())
            .for_each(|(self_i, tensor_rank_3)| *self_i += tensor_rank_3);
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize, const W: usize>
    AddAssign<&Self> for TensorRank3List<D, I, J, K, W>
{
    fn add_assign(&mut self, tensor_rank_3_list: &Self) {
        self.iter_mut()
            .zip(tensor_rank_3_list.iter())
            .for_each(|(self_i, tensor_rank_3)| *self_i += tensor_rank_3);
    }
}
