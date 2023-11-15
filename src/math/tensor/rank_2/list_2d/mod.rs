#[cfg(test)]
mod test;

use std::ops::
{
    Index,
    IndexMut,
    Mul
};

use super::
{
    TensorRank0,
    TensorRank2,
    list::
    {
        TensorRank2List,
        TensorRank2ListTrait
    }
};

/// A 2D list of *d*-dimensional tensors of rank 2.
///
/// `D` is the dimension, `I`, `J` are the configurations `W` is the list length.
pub struct TensorRank2List2D<const D: usize, const I: usize, const J: usize, const W: usize>
(
    /// An array of rank-2 tensors.
    pub [TensorRank2List<D, I, J, W>; W]
);

/// Inherent implementation of [`TensorRank2List2D`].
impl<const D: usize, const I: usize, const J: usize, const W: usize> TensorRank2List2D<D, I, J, W>
{
    /// Returns an iterator.
    ///
    /// The iterator yields all items from start to end. [Read more](https://doc.rust-lang.org/std/iter/)
    pub fn iter(&self) -> impl Iterator<Item = &TensorRank2List<D, I, J, W>>
    {
        self.0.iter()
    }
    /// Returns an iterator that allows modifying each value.
    ///
    /// The iterator yields all items from start to end. [Read more](https://doc.rust-lang.org/std/iter/)
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut TensorRank2List<D, I, J, W>>
    {
        self.0.iter_mut()
    }
}

/// Required methods for 2D rank-2 tensor lists.
pub trait TensorRank2List2DTrait<const D: usize, const W: usize>
{
    /// Returns a list of rank-2 tensors given an array.
    fn new(array: [[[[TensorRank0; D]; D]; W]; W]) -> Self;
    /// Returns a list of rank-2 zero tensors.
    fn zero() -> Self;
}

/// Implementation of [`TensorRank2List2DTrait`] for [`TensorRank2List2D`].
impl<const D: usize, const I: usize, const J: usize, const W: usize> TensorRank2List2DTrait<D, W> for TensorRank2List2D<D, I, J, W>
{
    fn new(array: [[[[TensorRank0; D]; D]; W]; W]) -> Self
    {
        array.iter().map(|array_i|
            TensorRank2List::new(*array_i)
        ).collect()
    }
    fn zero() -> Self
    {
        Self(std::array::from_fn(|_| TensorRank2List::zero()))
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> FromIterator<TensorRank2List<D, I, J, W>> for TensorRank2List2D<D, I, J, W>
{
    fn from_iter<Ii: IntoIterator<Item=TensorRank2List<D, I, J, W>>>(into_iterator: Ii) -> Self
    {
        let mut tensor_rank_2_list_2d = Self::zero();
        tensor_rank_2_list_2d.iter_mut().zip(into_iterator).for_each(|(tensor_rank_2_list, entry)|
            *tensor_rank_2_list = entry
        );
        tensor_rank_2_list_2d
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> Index<usize> for TensorRank2List2D<D, I, J, W>
{
    type Output = TensorRank2List<D, I, J, W>;
    fn index(&self, index: usize) -> &Self::Output
    {
        &self.0[index]
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> IndexMut<usize> for TensorRank2List2D<D, I, J, W>
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output
    {
        &mut self.0[index]
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize, const W: usize> Mul<TensorRank2<D, J, K>> for TensorRank2List2D<D, I, J, W>
{
    type Output = TensorRank2List2D<D, I, K, W>;
    fn mul(self, tensor_rank_2: TensorRank2<D, J, K>) -> Self::Output
    {
        self.iter().map(|self_entry|
            self_entry.iter().map(|self_tensor_rank_2|
                self_tensor_rank_2 * &tensor_rank_2
            ).collect()
        ).collect()
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize, const W: usize> Mul<&TensorRank2<D, J, K>> for TensorRank2List2D<D, I, J, W>
{
    type Output = TensorRank2List2D<D, I, K, W>;
    fn mul(self, tensor_rank_2: &TensorRank2<D, J, K>) -> Self::Output
    {
        self.iter().map(|self_entry|
            self_entry.iter().map(|self_tensor_rank_2|
                self_tensor_rank_2 * tensor_rank_2
            ).collect()
        ).collect()
    }
}