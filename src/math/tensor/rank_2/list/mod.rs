#[cfg(test)]
mod test;

use std::{array::from_fn, ops::{Add, AddAssign, Index, IndexMut}};

use super::{TensorRank0, TensorRank2, Tensor, super::Tensors};

/// A list of *d*-dimensional tensors of rank 2.
///
/// `D` is the dimension, `I`, `J` are the configurations `W` is the list length.
pub struct TensorRank2List<const D: usize, const I: usize, const J: usize, const W: usize>(
    pub [TensorRank2<D, I, J>; W],
);

/// Inherent implementation of [`TensorRank2List`].
impl<const D: usize, const I: usize, const J: usize, const W: usize> TensorRank2List<D, I, J, W> {
    /// Returns a list of rank-2 identity tensors.
    pub fn identity() -> Self {
        Self(from_fn(|_| TensorRank2::identity()))
    }
}

/// Implementation of [`Tensors`] for [`TensorRank2List`].
impl<const D: usize, const I: usize, const J: usize, const W: usize> Tensors for TensorRank2List<D, I, J, W> {
    type Array = [[[TensorRank0; D]; D]; W];
    type Item = TensorRank2<D, I, J>;
    fn as_array(&self) -> Self::Array {
        let mut array = [[[0.0; D]; D]; W];
        array
            .iter_mut()
            .zip(self.iter())
            .for_each(|(entry_rank_2, tensor_rank_2)| *entry_rank_2 = tensor_rank_2.as_array());
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
            .map(|array_i| TensorRank2::new(*array_i))
            .collect()
    }
    fn zero() -> Self {
        Self(from_fn(|_| TensorRank2::zero()))
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize>
    FromIterator<TensorRank2<D, I, J>> for TensorRank2List<D, I, J, W>
{
    fn from_iter<Ii: IntoIterator<Item = TensorRank2<D, I, J>>>(into_iterator: Ii) -> Self {
        let mut tensor_rank_2_list = Self::zero();
        tensor_rank_2_list
            .iter_mut()
            .zip(into_iterator)
            .for_each(|(tensor_rank_2_list_entry, entry)| *tensor_rank_2_list_entry = entry);
        tensor_rank_2_list
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> Index<usize>
    for TensorRank2List<D, I, J, W>
{
    type Output = TensorRank2<D, I, J>;
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> IndexMut<usize>
    for TensorRank2List<D, I, J, W>
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> Add
    for TensorRank2List<D, I, J, W>
{
    type Output = Self;
    fn add(mut self, tensor_rank_2_list: Self) -> Self::Output {
        self += tensor_rank_2_list;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> Add<&Self>
    for TensorRank2List<D, I, J, W>
{
    type Output = Self;
    fn add(mut self, tensor_rank_2_list: &Self) -> Self::Output {
        self += tensor_rank_2_list;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize>
    Add<TensorRank2List<D, I, J, W>> for &TensorRank2List<D, I, J, W>
{
    type Output = TensorRank2List<D, I, J, W>;
    fn add(self, mut tensor_rank_2_list: TensorRank2List<D, I, J, W>) -> Self::Output {
        tensor_rank_2_list += self;
        tensor_rank_2_list
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> AddAssign
    for TensorRank2List<D, I, J, W>
{
    fn add_assign(&mut self, tensor_rank_2_list: Self) {
        self.iter_mut()
            .zip(tensor_rank_2_list.iter())
            .for_each(|(self_entry, tensor_rank_2)| *self_entry += tensor_rank_2);
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> AddAssign<&Self>
    for TensorRank2List<D, I, J, W>
{
    fn add_assign(&mut self, tensor_rank_2_list: &Self) {
        self.iter_mut()
            .zip(tensor_rank_2_list.iter())
            .for_each(|(self_entry, tensor_rank_2)| *self_entry += tensor_rank_2);
    }
}
