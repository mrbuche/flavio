#[cfg(test)]
mod test;

pub mod list;
pub mod list_2d;
pub mod list_3d;

use std::{
    array::from_fn,
    fmt::{Display, Formatter, Result},
    ops::{Add, AddAssign, Div, DivAssign, Index, IndexMut, Mul, MulAssign, Sub, SubAssign},
};

use super::{rank_0::TensorRank0, rank_2::TensorRank2, Tensor};

/// Returns the rank-3 Levi-Civita symbol.
pub fn levi_civita<const I: usize, const J: usize, const K: usize>() -> TensorRank3<3, I, J, K> {
    TensorRank3::new([
        [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]],
        [[0.0, 0.0, -1.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
    ])
}

/// A *d*-dimensional tensor of rank 3.
///
/// `D` is the dimension, `I`, `J`, `K` are the configurations.
#[derive(Debug)]
pub struct TensorRank3<const D: usize, const I: usize, const J: usize, const K: usize>(
    pub [TensorRank2<D, J, K>; D],
);

impl<const D: usize, const I: usize, const J: usize, const K: usize> Display
    for TensorRank3<D, I, J, K>
{
    fn fmt(&self, f: &mut Formatter) -> Result {
        todo!()
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> PartialEq
    for TensorRank3<D, I, J, K>
{
    fn eq(&self, other: &Self) -> bool {
        let mut result = true;
        self.iter().zip(other.iter()).for_each(|(self_i, other_i)| {
            if self_i != other_i {
                result = false
            }
        });
        result
    }
}

/// Implementation of [`Tensor`] for [`TensorRank3`].
impl<const D: usize, const I: usize, const J: usize, const K: usize> Tensor
    for TensorRank3<D, I, J, K>
{
    type Array = [[[TensorRank0; D]; D]; D];
    type Item = TensorRank2<D, J, K>;
    fn as_array(&self) -> Self::Array {
        let mut array = [[[0.0; D]; D]; D];
        array
            .iter_mut()
            .zip(self.iter())
            .for_each(|(entry_rank_2, tensor_rank_2)| *entry_rank_2 = tensor_rank_2.as_array());
        array
    }
    fn copy(&self) -> Self {
        self.iter()
            .map(|entry_rank_2| entry_rank_2.copy())
            .collect()
    }
    fn identity() -> Self {
        panic!()
    }
    fn is_zero(&self) -> bool {
        self.iter()
            .map(|entry_rank_2| {
                entry_rank_2
                    .iter()
                    .map(|entry_rank_1| {
                        entry_rank_1
                            .iter()
                            .map(|entry_rank_0| (entry_rank_0 == &0.0) as u8)
                            .sum::<u8>()
                    })
                    .sum::<u8>()
            })
            .sum::<u8>()
            == ((D * D * D) as u8)
    }
    fn iter(&self) -> impl Iterator<Item = &Self::Item> {
        self.0.iter()
    }
    fn iter_mut(&mut self) -> impl Iterator<Item = &mut Self::Item> {
        self.0.iter_mut()
    }
    fn new(array: Self::Array) -> Self {
        array.iter().map(|entry| TensorRank2::new(*entry)).collect()
    }
    fn norm_squared(&self) -> TensorRank0 {
        panic!()
    }
    fn normalized(&self) -> Self {
        panic!()
    }
    fn zero() -> Self {
        Self(from_fn(|_| TensorRank2::zero()))
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize>
    FromIterator<TensorRank2<D, J, K>> for TensorRank3<D, I, J, K>
{
    fn from_iter<Ii: IntoIterator<Item = TensorRank2<D, J, K>>>(into_iterator: Ii) -> Self {
        let mut tensor_rank_3 = Self::zero();
        tensor_rank_3
            .iter_mut()
            .zip(into_iterator)
            .for_each(|(tensor_rank_3_i, value_i)| *tensor_rank_3_i = value_i);
        tensor_rank_3
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> Index<usize>
    for TensorRank3<D, I, J, K>
{
    type Output = TensorRank2<D, J, K>;
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> IndexMut<usize>
    for TensorRank3<D, I, J, K>
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> Div<TensorRank0>
    for TensorRank3<D, I, J, K>
{
    type Output = Self;
    fn div(mut self, tensor_rank_0: TensorRank0) -> Self::Output {
        self /= &tensor_rank_0;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> Div<&TensorRank0>
    for TensorRank3<D, I, J, K>
{
    type Output = Self;
    fn div(mut self, tensor_rank_0: &TensorRank0) -> Self::Output {
        self /= tensor_rank_0;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> DivAssign<TensorRank0>
    for TensorRank3<D, I, J, K>
{
    fn div_assign(&mut self, tensor_rank_0: TensorRank0) {
        self.iter_mut().for_each(|self_i| *self_i /= &tensor_rank_0);
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> DivAssign<&TensorRank0>
    for TensorRank3<D, I, J, K>
{
    fn div_assign(&mut self, tensor_rank_0: &TensorRank0) {
        self.iter_mut().for_each(|self_i| *self_i /= tensor_rank_0);
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> Mul<TensorRank0>
    for TensorRank3<D, I, J, K>
{
    type Output = Self;
    fn mul(mut self, tensor_rank_0: TensorRank0) -> Self::Output {
        self *= &tensor_rank_0;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> Mul<&TensorRank0>
    for TensorRank3<D, I, J, K>
{
    type Output = Self;
    fn mul(mut self, tensor_rank_0: &TensorRank0) -> Self::Output {
        self *= tensor_rank_0;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> MulAssign<TensorRank0>
    for TensorRank3<D, I, J, K>
{
    fn mul_assign(&mut self, tensor_rank_0: TensorRank0) {
        self.iter_mut().for_each(|self_i| *self_i *= &tensor_rank_0);
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> MulAssign<&TensorRank0>
    for TensorRank3<D, I, J, K>
{
    fn mul_assign(&mut self, tensor_rank_0: &TensorRank0) {
        self.iter_mut().for_each(|self_i| *self_i *= tensor_rank_0);
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> Add
    for TensorRank3<D, I, J, K>
{
    type Output = Self;
    fn add(mut self, tensor_rank_3: Self) -> Self::Output {
        self += tensor_rank_3;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> Add<&Self>
    for TensorRank3<D, I, J, K>
{
    type Output = Self;
    fn add(mut self, tensor_rank_3: &Self) -> Self::Output {
        self += tensor_rank_3;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> Add<TensorRank3<D, I, J, K>>
    for &TensorRank3<D, I, J, K>
{
    type Output = TensorRank3<D, I, J, K>;
    fn add(self, mut tensor_rank_3: TensorRank3<D, I, J, K>) -> Self::Output {
        tensor_rank_3 += self;
        tensor_rank_3
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> AddAssign
    for TensorRank3<D, I, J, K>
{
    fn add_assign(&mut self, tensor_rank_3: Self) {
        self.iter_mut()
            .zip(tensor_rank_3.iter())
            .for_each(|(self_i, tensor_rank_3_i)| *self_i += tensor_rank_3_i);
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> AddAssign<&Self>
    for TensorRank3<D, I, J, K>
{
    fn add_assign(&mut self, tensor_rank_3: &Self) {
        self.iter_mut()
            .zip(tensor_rank_3.iter())
            .for_each(|(self_i, tensor_rank_3_i)| *self_i += tensor_rank_3_i);
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> Sub
    for TensorRank3<D, I, J, K>
{
    type Output = Self;
    fn sub(mut self, tensor_rank_3: Self) -> Self::Output {
        self -= tensor_rank_3;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> Sub<&Self>
    for TensorRank3<D, I, J, K>
{
    type Output = Self;
    fn sub(mut self, tensor_rank_3: &Self) -> Self::Output {
        self -= tensor_rank_3;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> SubAssign
    for TensorRank3<D, I, J, K>
{
    fn sub_assign(&mut self, tensor_rank_3: Self) {
        self.iter_mut()
            .zip(tensor_rank_3.iter())
            .for_each(|(self_i, tensor_rank_3_i)| *self_i -= tensor_rank_3_i);
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> SubAssign<&Self>
    for TensorRank3<D, I, J, K>
{
    fn sub_assign(&mut self, tensor_rank_3: &Self) {
        self.iter_mut()
            .zip(tensor_rank_3.iter())
            .for_each(|(self_i, tensor_rank_3_i)| *self_i -= tensor_rank_3_i);
    }
}
