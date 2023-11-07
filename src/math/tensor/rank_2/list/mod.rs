#[cfg(test)]
mod test;

use std::ops::
{
    Index,
    IndexMut
};

use super::
{
    TensorRank0,
    TensorRank1,
    TensorRank1Trait,
    TensorRank2
};

/// A list of *d*-dimensional tensors of rank 2.
///
/// `D` is the dimension, `I`, `J` are the configurations `L` is the list length.
pub struct TensorRank2List<const D: usize, const I: usize, const J: usize, const L: usize>
(
    [TensorRank2<D, I, J>; L]
);

impl<const D: usize, const I: usize, const J: usize, const L: usize> TensorRank2List<D, I, J, L>
{
    /// Returns an iterator.
    ///
    /// The iterator yields all items from start to end. [Read more](https://doc.rust-lang.org/std/iter/)
    pub fn iter(&self) -> impl Iterator<Item = &TensorRank2<D, I, J>>
    {
        self.0.iter()
    }
    /// Returns an iterator that allows modifying each value.
    ///
    /// The iterator yields all items from start to end. [Read more](https://doc.rust-lang.org/std/iter/)
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut TensorRank2<D, I, J>>
    {
        self.0.iter_mut()
    }
}

/// Required methods for rank-2 tensor lists.
pub trait TensorRank2ListTrait<const D: usize, const L: usize>
{
    /// Returns a list of rank-2 tensors given an array.
    fn new(array: [[[TensorRank0; D]; D]; L]) -> Self;
    /// Returns a list of rank-2 zero tensors.
    fn zero() -> Self;
}

/// Implementation of [`TensorRank2ListTrait`] for [`TensorRank2List`].
impl<const D: usize, const I: usize, const J: usize, const L: usize> TensorRank2ListTrait<D, L> for TensorRank2List<D, I, J, L>
{
    fn new(array: [[[TensorRank0; D]; D]; L]) -> Self
    {
        array.iter().map(|array_i|
            array_i.iter().map(|array_ij|
                TensorRank1::new(*array_ij)
            ).collect()
        ).collect()
    }
    fn zero() -> Self
    {
        Self(std::array::from_fn(|_| TensorRank2::zero()))
    }
}

impl<const D: usize, const I: usize, const J: usize, const L: usize> FromIterator<TensorRank2<D, I, J>> for TensorRank2List<D, I, J, L>
{
    fn from_iter<Ii: IntoIterator<Item=TensorRank2<D, I, J>>>(into_iterator: Ii) -> Self
    {
        let mut tensor_rank_2_list = Self::zero();
        tensor_rank_2_list.iter_mut().zip(into_iterator).for_each(|(tensor_rank_2_list_entry, entry)|
            *tensor_rank_2_list_entry = entry
        );
        tensor_rank_2_list
    }
}

impl<const D: usize, const I: usize, const J: usize, const L: usize> Index<usize> for TensorRank2List<D, I, J, L>
{
    type Output = TensorRank2<D, I, J>;
    fn index(&self, index: usize) -> &Self::Output
    {
        &self.0[index]
    }
}

impl<const D: usize, const I: usize, const J: usize, const L: usize> IndexMut<usize> for TensorRank2List<D, I, J, L>
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output
    {
        &mut self.0[index]
    }
}
