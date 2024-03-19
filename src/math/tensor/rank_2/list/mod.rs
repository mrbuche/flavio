#[cfg(test)]
mod test;

use std::ops::
{
    Add,
    AddAssign,
    Div,
    DivAssign,
    Index,
    IndexMut
};

use super::
{
    TensorRank0,
    TensorRank2,
    TensorRank2Trait
};

/// A list of *d*-dimensional tensors of rank 2.
///
/// `D` is the dimension, `I`, `J` are the configurations `W` is the list length.
pub struct TensorRank2List<const D: usize, const I: usize, const J: usize, const W: usize>
(
    [TensorRank2<D, I, J>; W]
);

/// Inherent implementation of [`TensorRank2List`].
impl<const D: usize, const I: usize, const J: usize, const W: usize> TensorRank2List<D, I, J, W>
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
pub trait TensorRank2ListTrait<const D: usize, const W: usize>
{
    /// Returns the rank-2 tensor list as an array.
    fn as_array(&self) -> [[[TensorRank0; D]; D]; W];
    /// Returns a list of rank-2 tensors given an array.
    fn new(array: [[[TensorRank0; D]; D]; W]) -> Self;
    /// Returns a list of rank-2 zero tensors.
    fn zero() -> Self;
}

/// Implementation of [`TensorRank2ListTrait`] for [`TensorRank2List`].
impl<const D: usize, const I: usize, const J: usize, const W: usize> TensorRank2ListTrait<D, W> for TensorRank2List<D, I, J, W>
{
    fn as_array(&self) -> [[[TensorRank0; D]; D]; W]
    {
        let mut array = [[[0.0; D]; D]; W];
        array.iter_mut()
        .zip(self.iter())
        .for_each(|(entry_rank_2, tensor_rank_2)|
            *entry_rank_2 = tensor_rank_2.as_array()
        );
        array
    }
    fn new(array: [[[TensorRank0; D]; D]; W]) -> Self
    {
        array.iter().map(|array_i|
            TensorRank2::new(*array_i)
        ).collect()
    }
    fn zero() -> Self
    {
        Self(std::array::from_fn(|_| TensorRank2::zero()))
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> FromIterator<TensorRank2<D, I, J>> for TensorRank2List<D, I, J, W>
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

impl<const D: usize, const I: usize, const J: usize, const W: usize> Index<usize> for TensorRank2List<D, I, J, W>
{
    type Output = TensorRank2<D, I, J>;
    fn index(&self, index: usize) -> &Self::Output
    {
        &self.0[index]
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> IndexMut<usize> for TensorRank2List<D, I, J, W>
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output
    {
        &mut self.0[index]
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> std::iter::Sum for TensorRank2List<D, I, J, W>
{
    fn sum<Ii>(iter: Ii) -> Self
    where
        Ii: Iterator<Item = Self>
    {
        let mut output = TensorRank2List::zero();
        iter.for_each(|item|
            output += item
        );
        output
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> Div<TensorRank0> for TensorRank2List<D, I, J, W>
{
    type Output = Self;
    fn div(mut self, tensor_rank_0: TensorRank0) -> Self::Output
    {
        self /= &tensor_rank_0;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> Div<TensorRank0> for &TensorRank2List<D, I, J, W>
{
    type Output = TensorRank2List<D, I, J, W>;
    fn div(self, tensor_rank_0: TensorRank0) -> Self::Output
    {
        self.iter().map(|self_i|
            self_i / tensor_rank_0
        ).collect()
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> Div<&TensorRank0> for TensorRank2List<D, I, J, W>
{
    type Output = Self;
    fn div(mut self, tensor_rank_0: &TensorRank0) -> Self::Output
    {
        self /= tensor_rank_0;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> Div<&TensorRank0> for &TensorRank2List<D, I, J, W>
{
    type Output = TensorRank2List<D, I, J, W>;
    fn div(self, tensor_rank_0: &TensorRank0) -> Self::Output
    {
        self.iter().map(|self_i|
            self_i / tensor_rank_0
        ).collect()
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> DivAssign<TensorRank0> for TensorRank2List<D, I, J, W>
{
    fn div_assign(&mut self, tensor_rank_0: TensorRank0)
    {
        self.iter_mut().for_each(|entry|
            *entry /= &tensor_rank_0
        );
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> DivAssign<&TensorRank0> for TensorRank2List<D, I, J, W>
{
    fn div_assign(&mut self, tensor_rank_0: &TensorRank0)
    {
        self.iter_mut().for_each(|entry|
            *entry /= tensor_rank_0
        );
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> Add for TensorRank2List<D, I, J, W>
{
    type Output = Self;
    fn add(mut self, tensor_rank_2_list: Self) -> Self::Output
    {
        self += tensor_rank_2_list;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> Add<&Self> for TensorRank2List<D, I, J, W>
{
    type Output = Self;
    fn add(mut self, tensor_rank_2_list: &Self) -> Self::Output
    {
        self += tensor_rank_2_list;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> Add<TensorRank2List<D, I, J, W>> for &TensorRank2List<D, I, J, W>
{
    type Output = TensorRank2List<D, I, J, W>;
    fn add(self, mut tensor_rank_2_list: TensorRank2List<D, I, J, W>) -> Self::Output
    {
        tensor_rank_2_list += self;
        tensor_rank_2_list
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> AddAssign for TensorRank2List<D, I, J, W>
{
    fn add_assign(&mut self, tensor_rank_2_list: Self)
    {
        self.iter_mut().zip(tensor_rank_2_list.iter())
        .for_each(|(self_entry, tensor_rank_2)|
            *self_entry += tensor_rank_2
        );
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> AddAssign<&Self> for TensorRank2List<D, I, J, W>
{
    fn add_assign(&mut self, tensor_rank_2_list: &Self)
    {
        self.iter_mut().zip(tensor_rank_2_list.iter())
        .for_each(|(self_entry, tensor_rank_2)|
            *self_entry += tensor_rank_2
        );
    }
}