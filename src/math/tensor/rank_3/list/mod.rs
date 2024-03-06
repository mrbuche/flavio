#[cfg(test)]
pub mod test;

use super::
{
    TensorRank0,
    TensorRank3,
    TensorRank3Trait
};

/// A list of *d*-dimensional tensors of rank 3.
///
/// `D` is the dimension, `I`, `J`, `K` are the configurations `W` is the list length.
pub struct TensorRank3List<const D: usize, const I: usize, const J: usize, const K: usize, const W: usize>
(
    [TensorRank3<D, I, J, K>; W]
);

/// Inherent implementation of [`TensorRank3List`].
impl<const D: usize, const I: usize, const J: usize, const K: usize, const W: usize> TensorRank3List<D, I, J, K, W>
{
    /// Returns an iterator.
    ///
    /// The iterator yields all items from start to end. [Read more](https://doc.rust-lang.org/std/iter/)
    pub fn iter(&self) -> impl Iterator<Item = &TensorRank3<D, I, J, K>>
    {
        self.0.iter()
    }
    /// Returns an iterator that allows modifying each value.
    ///
    /// The iterator yields all items from start to end. [Read more](https://doc.rust-lang.org/std/iter/)
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut TensorRank3<D, I, J, K>>
    {
        self.0.iter_mut()
    }
}

/// Required methods for rank-3 tensor lists.
pub trait TensorRank3ListTrait<const D: usize, const W: usize>
{
    /// Returns the rank-3 tensor list as an array.
    fn as_array(&self) -> [[[[TensorRank0; D]; D]; D]; W];
    /// Returns a list of rank-3 tensors given an array.
    fn new(array: [[[[TensorRank0; D]; D]; D]; W]) -> Self;
    /// Returns a list of rank-3 zero tensors.
    fn zero() -> Self;
}

/// Implementation of [`TensorRank3ListTrait`] for [`TensorRank3List`].
impl<const D: usize, const I: usize, const J: usize, const K: usize, const W: usize> TensorRank3ListTrait<D, W> for TensorRank3List<D, I, J, K, W>
{
    fn as_array(&self) -> [[[[TensorRank0; D]; D]; D]; W]
    {
        let mut array = [[[[0.0; D]; D]; D]; W];
        array.iter_mut()
        .zip(self.iter())
        .for_each(|(entry_rank_3, tensor_rank_3)|
            *entry_rank_3 = tensor_rank_3.as_array()
        );
        array
    }
    fn new(array: [[[[TensorRank0; D]; D]; D]; W]) -> Self
    {
        array.iter().map(|array_i|
            TensorRank3::new(*array_i)
        ).collect()
    }
    fn zero() -> Self
    {
        Self(std::array::from_fn(|_| TensorRank3::zero()))
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize, const W: usize> FromIterator<TensorRank3<D, I, J, K>> for TensorRank3List<D, I, J, K, W>
{
    fn from_iter<Ii: IntoIterator<Item=TensorRank3<D, I, J, K>>>(into_iterator: Ii) -> Self
    {
        let mut tensor_rank_3_list = Self::zero();
        tensor_rank_3_list.iter_mut().zip(into_iterator).for_each(|(tensor_rank_3_list_entry, entry)|
            *tensor_rank_3_list_entry = entry
        );
        tensor_rank_3_list
    }
}