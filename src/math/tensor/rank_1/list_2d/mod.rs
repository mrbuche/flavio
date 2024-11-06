#[cfg(test)]
mod test;

use super::{super::Tensors, list::TensorRank1List, TensorRank0};
use std::array::from_fn;
use std::{
    fmt::{Display, Formatter, Result},
    ops::{Add, AddAssign},
};

/// A 2D list of *d*-dimensional tensors of rank 1.
///
/// `D` is the dimension, `I` is the configuration, `W` and `X` are the list lengths.
#[derive(Debug)]
pub struct TensorRank1List2D<const D: usize, const I: usize, const W: usize, const X: usize>(
    pub [TensorRank1List<D, I, W>; X],
);

impl<const D: usize, const I: usize, const W: usize, const X: usize> Display
    for TensorRank1List2D<D, I, W, X>
{
    fn fmt(&self, _f: &mut Formatter) -> Result {
        Ok(())
    }
}

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
    fn copy(&self) -> Self {
        self.iter().map(|entry| entry.copy()).collect()
    }
    fn full_contraction(&self, tensor_rank_1_list_2d: &Self) -> TensorRank0 {
        self.iter()
            .zip(tensor_rank_1_list_2d.iter())
            .map(|(self_entry, tensor_rank_1_list)| self_entry.full_contraction(tensor_rank_1_list))
            .sum()
    }
    fn identity() -> Self {
        Self(from_fn(|_| Self::Item::identity()))
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

impl<const D: usize, const I: usize, const W: usize, const X: usize> Add
    for TensorRank1List2D<D, I, W, X>
{
    type Output = Self;
    fn add(mut self, tensor_rank_1_list_2d: Self) -> Self::Output {
        self += tensor_rank_1_list_2d;
        self
    }
}

impl<const D: usize, const I: usize, const W: usize, const X: usize> Add<&Self>
    for TensorRank1List2D<D, I, W, X>
{
    type Output = Self;
    fn add(mut self, tensor_rank_1_list_2d: &Self) -> Self::Output {
        self += tensor_rank_1_list_2d;
        self
    }
}

impl<const D: usize, const I: usize, const W: usize, const X: usize> AddAssign
    for TensorRank1List2D<D, I, W, X>
{
    fn add_assign(&mut self, tensor_rank_1_list_2d: Self) {
        self.iter_mut()
            .zip(tensor_rank_1_list_2d.iter())
            .for_each(|(self_entry, tensor_rank_1_list)| *self_entry += tensor_rank_1_list);
    }
}

impl<const D: usize, const I: usize, const W: usize, const X: usize> AddAssign<&Self>
    for TensorRank1List2D<D, I, W, X>
{
    fn add_assign(&mut self, tensor_rank_1_list_2d: &Self) {
        self.iter_mut()
            .zip(tensor_rank_1_list_2d.iter())
            .for_each(|(self_entry, tensor_rank_1_list)| *self_entry += tensor_rank_1_list);
    }
}
