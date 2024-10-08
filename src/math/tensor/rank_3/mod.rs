#[cfg(test)]
mod test;

pub mod list;
pub mod list_2d;
pub mod list_3d;

use std::ops::{Add, AddAssign, Div, DivAssign, Index, IndexMut, Mul, MulAssign, Sub, SubAssign};

use super::{
    rank_0::TensorRank0,
    rank_2::{TensorRank2, TensorRank2Trait},
};

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
pub struct TensorRank3<const D: usize, const I: usize, const J: usize, const K: usize>(
    [TensorRank2<D, J, K>; D],
);

/// Inherent implementation of [`TensorRank3`].
impl<const D: usize, const I: usize, const J: usize, const K: usize> TensorRank3<D, I, J, K> {
    /// Returns an iterator.
    ///
    /// The iterator yields all items from start to end. [Read more](https://doc.rust-lang.org/std/iter/)
    pub fn iter(&self) -> impl Iterator<Item = &TensorRank2<D, J, K>> {
        self.0.iter()
    }
    /// Returns an iterator that allows modifying each value.
    ///
    /// The iterator yields all items from start to end. [Read more](https://doc.rust-lang.org/std/iter/)
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut TensorRank2<D, J, K>> {
        self.0.iter_mut()
    }
    /// Returns the rank-3 zero tensor.
    pub fn zero() -> Self {
        Self(std::array::from_fn(|_| TensorRank2::zero()))
    }
}

/// Required methods for rank-3 tensors.
pub trait TensorRank3Trait<const D: usize> {
    /// Returns the rank-3 tensor as an array.
    fn as_array(&self) -> [[[TensorRank0; D]; D]; D];
    /// Returns a rank-3 tensor given an array.
    fn new(array: [[[TensorRank0; D]; D]; D]) -> Self;
}

/// Implementation of [`TensorRank3Trait`] for [`TensorRank3`].
impl<const D: usize, const I: usize, const J: usize, const K: usize> TensorRank3Trait<D>
    for TensorRank3<D, I, J, K>
{
    fn as_array(&self) -> [[[TensorRank0; D]; D]; D] {
        let mut array = [[[0.0; D]; D]; D];
        array
            .iter_mut()
            .zip(self.iter())
            .for_each(|(entry_rank_2, tensor_rank_2)| *entry_rank_2 = tensor_rank_2.as_array());
        array
    }
    fn new(array: [[[TensorRank0; D]; D]; D]) -> Self {
        array
            .iter()
            .map(|array_i| TensorRank2::new(*array_i))
            .collect()
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
