#[cfg(test)]
mod test;

use std::ops::
{
    Index,
    IndexMut
};
use std::array::from_fn;

use super::
{
    TensorRank0,
    TensorRank4,
    TensorRank4Trait
};

/// A list of *d*-dimensional tensor of rank 4.
///
/// `D` is the dimension, `I`, `J`, `K`, `L` are the configurations, `W` is the list length.
pub struct TensorRank4List<const D: usize, const I: usize, const J: usize, const K: usize, const L: usize, const W: usize>
(
    [TensorRank4<D, I, J, K, L>; W]
);

/// Inherent implementation of [`TensorRank4List`].
impl<const D: usize, const I: usize, const J: usize, const K: usize, const L: usize, const W: usize> TensorRank4List<D, I, J, K, L, W>
{
    /// Returns an iterator.
    ///
    /// The iterator yields all items from start to end. [Read more](https://doc.rust-lang.org/std/iter/)
    pub fn iter(&self) -> impl Iterator<Item = &TensorRank4<D, I, J, K, L>>
    {
        self.0.iter()
    }
    /// Returns an iterator that allows modifying each value.
    ///
    /// The iterator yields all items from start to end. [Read more](https://doc.rust-lang.org/std/iter/)
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut TensorRank4<D, I, J, K, L>>
    {
        self.0.iter_mut()
    }
    /// Returns a list of rank-4 zero tensors.
    pub fn zero() -> Self
    {
        Self(std::array::from_fn(|_| TensorRank4::zero()))
    }
}

/// Required methods for rank-4 tensor lists.
pub trait TensorRank4ListTrait<const D: usize, const W: usize>
{
    /// Returns the rank-4 tensor list as an array.
    fn as_array(&self) -> [[[[[TensorRank0; D]; D]; D]; D]; W];
    /// Returns a list of rank-4 tensors given an array.
    fn new(array: [[[[[TensorRank0; D]; D]; D]; D]; W]) -> Self;
    /// Returns a list of rank-4 zero tensors.
    fn zero() -> Self;
}

/// Implementation of [`TensorRank4ListTrait`] for [`TensorRank4List`].
impl<const D: usize, const I: usize, const J: usize, const K: usize, const L: usize, const W: usize> TensorRank4ListTrait<D, W> for TensorRank4List<D, I, J, K, L, W>
{
    fn as_array(&self) -> [[[[[TensorRank0; D]; D]; D]; D]; W]
    {
        let mut array = [[[[[0.0; D]; D]; D]; D]; W];
        array.iter_mut()
        .zip(self.iter())
        .for_each(|(entry_rank_4, tensor_rank_4)|
            *entry_rank_4 = tensor_rank_4.as_array()
        );
        array
    }
    fn new(array: [[[[[TensorRank0; D]; D]; D]; D]; W]) -> Self
    {
        array.iter().map(|array_i|
            TensorRank4::new(*array_i)
        ).collect()
    }
    fn zero() -> Self
    {
        Self(from_fn(|_| TensorRank4::zero()))
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize, const L: usize, const W: usize> FromIterator<TensorRank4<D, I, J, K, L>> for TensorRank4List<D, I, J, K, L, W>
{
    fn from_iter<Ii: IntoIterator<Item=TensorRank4<D, I, J, K, L>>>(into_iterator: Ii) -> Self
    {
        let mut tensor_rank_4_list = Self::zero();
        tensor_rank_4_list.iter_mut()
        .zip(into_iterator)
        .for_each(|(tensor_rank_4_list_entry, entry)|
            *tensor_rank_4_list_entry = entry
        );
        tensor_rank_4_list
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize, const L: usize, const W: usize> Index<usize> for TensorRank4List<D, I, J, K, L, W>
{
    type Output = TensorRank4<D, I, J, K, L>;
    fn index(&self, index: usize) -> &Self::Output
    {
        &self.0[index]
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize, const L: usize, const W: usize> IndexMut<usize> for TensorRank4List<D, I, J, K, L, W>
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output
    {
        &mut self.0[index]
    }
}
