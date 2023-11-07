#[cfg(test)]
mod test;

use std::ops::
{
    Add,
    AddAssign,
    Div,
    DivAssign,
    Index,
    IndexMut,
    Mul,
    MulAssign,
    Sub,
    SubAssign
};
use std::array::from_fn;

use crate::math::
{
    TensorRank2,
    TensorRank2Trait
};

use super::
{
    TensorRank0,
    TensorRank1,
    TensorRank1Trait
};

/// A list of *d*-dimensional tensors of rank 1.
///
/// `D` is the dimension, `I` is the configuration, `L` is the list length.
pub struct TensorRank1List<const D: usize, const I: usize, const W: usize>
(
    [TensorRank1<D, I>; W]
);

impl<const D: usize, const I: usize, const W: usize> TensorRank1List<D, I, W>
{
    /// Returns an iterator.
    ///
    /// The iterator yields all items from start to end. [Read more](https://doc.rust-lang.org/std/iter/)
    pub fn iter(&self) -> impl Iterator<Item = &TensorRank1<D, I>>
    {
        self.0.iter()
    }
    /// Returns an iterator that allows modifying each value.
    ///
    /// The iterator yields all items from start to end. [Read more](https://doc.rust-lang.org/std/iter/)
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut TensorRank1<D, I>>
    {
        self.0.iter_mut()
    }
}

/// Required methods for rank-1 tensor lists.
pub trait TensorRank1ListTrait<const D: usize, const W: usize>
{
    /// Returns a list of rank-1 tensors given an array.
    fn new(array: [[TensorRank0; D]; W]) -> Self;
    /// Returns a list of rank-1 zero tensors.
    fn zero() -> Self;
}

/// Implementation of [`TensorRank1ListTrait`] for [`TensorRank1List`].
impl<const D: usize, const I: usize, const W: usize> TensorRank1ListTrait<D, W> for TensorRank1List<D, I, W>
{
    fn new(array: [[TensorRank0; D]; W]) -> Self
    {
        array.iter().map(|array_i|
            TensorRank1::new(*array_i)
        ).collect()
    }
    fn zero() -> Self
    {
        Self(from_fn(|_| TensorRank1::zero()))
    }
}

impl<const D: usize, const I: usize, const W: usize> FromIterator<TensorRank1<D, I>> for TensorRank1List<D, I, W>
{
    fn from_iter<Ii: IntoIterator<Item=TensorRank1<D, I>>>(into_iterator: Ii) -> Self
    {
        let mut tensor_rank_1_list = Self::zero();
        tensor_rank_1_list.iter_mut().zip(into_iterator).for_each(|(tensor_rank_1_list_entry, entry)|
            *tensor_rank_1_list_entry = entry
        );
        tensor_rank_1_list
    }
}

impl<const D: usize, const I: usize, const W: usize> Index<usize> for TensorRank1List<D, I, W>
{
    type Output = TensorRank1<D, I>;
    fn index(&self, index: usize) -> &Self::Output
    {
        &self.0[index]
    }
}

impl<const D: usize, const I: usize, const W: usize> IndexMut<usize> for TensorRank1List<D, I, W>
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output
    {
        &mut self.0[index]
    }
}

impl<const D: usize, const I: usize, const W: usize> Div<TensorRank0> for TensorRank1List<D, I, W>
{
    type Output = Self;
    fn div(mut self, tensor_rank_0: TensorRank0) -> Self::Output
    {
        self /= &tensor_rank_0;
        self
    }
}

impl<const D: usize, const I: usize, const W: usize> Div<TensorRank0> for &TensorRank1List<D, I, W>
{
    type Output = TensorRank1List<D, I, W>;
    fn div(self, tensor_rank_0: TensorRank0) -> Self::Output
    {
        self.iter().map(|self_i|
            self_i / tensor_rank_0
        ).collect()
    }
}

impl<const D: usize, const I: usize, const W: usize> Div<&TensorRank0> for TensorRank1List<D, I, W>
{
    type Output = Self;
    fn div(mut self, tensor_rank_0: &TensorRank0) -> Self::Output
    {
        self /= tensor_rank_0;
        self
    }
}

impl<const D: usize, const I: usize, const W: usize> Div<&TensorRank0> for &TensorRank1List<D, I, W>
{
    type Output = TensorRank1List<D, I, W>;
    fn div(self, tensor_rank_0: &TensorRank0) -> Self::Output
    {
        self.iter().map(|self_i|
            self_i / tensor_rank_0
        ).collect()
    }
}

impl<const D: usize, const I: usize, const W: usize> DivAssign<TensorRank0> for TensorRank1List<D, I, W>
{
    fn div_assign(&mut self, tensor_rank_0: TensorRank0)
    {
        self.iter_mut().for_each(|entry|
            *entry /= &tensor_rank_0
        );
    }
}

impl<const D: usize, const I: usize, const W: usize> DivAssign<&TensorRank0> for TensorRank1List<D, I, W>
{
    fn div_assign(&mut self, tensor_rank_0: &TensorRank0)
    {
        self.iter_mut().for_each(|entry|
            *entry /= tensor_rank_0
        );
    }
}

impl<const D: usize, const I: usize, const W: usize> Mul<TensorRank0> for TensorRank1List<D, I, W>
{
    type Output = Self;
    fn mul(mut self, tensor_rank_0: TensorRank0) -> Self::Output
    {
        self *= &tensor_rank_0;
        self
    }
}

impl<const D: usize, const I: usize, const W: usize> Mul<TensorRank0> for &TensorRank1List<D, I, W>
{
    type Output = TensorRank1List<D, I, W>;
    fn mul(self, tensor_rank_0: TensorRank0) -> Self::Output
    {
        self.iter().map(|self_i|
            self_i * tensor_rank_0
        ).collect()
    }
}

impl<const D: usize, const I: usize, const W: usize> Mul<&TensorRank0> for TensorRank1List<D, I, W>
{
    type Output = Self;
    fn mul(mut self, tensor_rank_0: &TensorRank0) -> Self::Output
    {
        self *= tensor_rank_0;
        self
    }
}

impl<const D: usize, const I: usize, const W: usize> Mul<&TensorRank0> for &TensorRank1List<D, I, W>
{
    type Output = TensorRank1List<D, I, W>;
    fn mul(self, tensor_rank_0: &TensorRank0) -> Self::Output
    {
        self.iter().map(|self_i|
            self_i * tensor_rank_0
        ).collect()
    }
}

impl<const D: usize, const I: usize, const W: usize> MulAssign<TensorRank0> for TensorRank1List<D, I, W>
{
    fn mul_assign(&mut self, tensor_rank_0: TensorRank0)
    {
        self.iter_mut().for_each(|entry|
            *entry *= &tensor_rank_0
        );
    }
}

impl<const D: usize, const I: usize, const W: usize> MulAssign<&TensorRank0> for TensorRank1List<D, I, W>
{
    fn mul_assign(&mut self, tensor_rank_0: &TensorRank0)
    {
        self.iter_mut().for_each(|entry|
            *entry *= tensor_rank_0
        );
    }
}

impl<const D: usize, const I: usize, const W: usize> Add for TensorRank1List<D, I, W>
{
    type Output = Self;
    fn add(mut self, tensor_rank_1_list: Self) -> Self::Output
    {
        self += tensor_rank_1_list;
        self
    }
}

impl<const D: usize, const I: usize, const W: usize> Add<&Self> for TensorRank1List<D, I, W>
{
    type Output = Self;
    fn add(mut self, tensor_rank_1_list: &Self) -> Self::Output
    {
        self += tensor_rank_1_list;
        self
    }
}

impl<const D: usize, const I: usize, const W: usize> Add<TensorRank1List<D, I, W>> for &TensorRank1List<D, I, W>
{
    type Output = TensorRank1List<D, I, W>;
    fn add(self, mut tensor_rank_1_list: TensorRank1List<D, I, W>) -> Self::Output
    {
        tensor_rank_1_list += self;
        tensor_rank_1_list
    }
}

impl<const D: usize, const I: usize, const W: usize> AddAssign for TensorRank1List<D, I, W>
{
    fn add_assign(&mut self, tensor_rank_1_list: Self)
    {
        self.iter_mut().zip(tensor_rank_1_list.iter()).for_each(|(self_entry, tensor_rank_1_list_entry)|
            *self_entry += tensor_rank_1_list_entry
        );
    }
}

impl<const D: usize, const I: usize, const W: usize> AddAssign<&Self> for TensorRank1List<D, I, W>
{
    fn add_assign(&mut self, tensor_rank_1_list: &Self)
    {
        self.iter_mut().zip(tensor_rank_1_list.iter()).for_each(|(self_entry, tensor_rank_1_list_entry)|
            *self_entry += tensor_rank_1_list_entry
        );
    }
}

impl<const I: usize, const J: usize, const W: usize> Mul<TensorRank1List<3, J, W>> for TensorRank1List<3, I, W>
{
    type Output = TensorRank2<3, I, J>;
    fn mul(self, tensor_rank_1_list: TensorRank1List<3, J, W>) -> Self::Output
    {
        self.iter().zip(tensor_rank_1_list.iter()).map(|(self_entry, tensor_rank_1_list_entry)|
            TensorRank2::dyad(self_entry, tensor_rank_1_list_entry)
        ).sum()
    }
}

impl<const I: usize, const J: usize, const W: usize> Mul<&TensorRank1List<3, J, W>> for TensorRank1List<3, I, W>
{
    type Output = TensorRank2<3, I, J>;
    fn mul(self, tensor_rank_1_list: &TensorRank1List<3, J, W>) -> Self::Output
    {
        self.iter().zip(tensor_rank_1_list.iter()).map(|(self_entry, tensor_rank_1_list_entry)|
            TensorRank2::dyad(self_entry, tensor_rank_1_list_entry)
        ).sum()
    }
}

impl<const I: usize, const J: usize, const W: usize> Mul<TensorRank1List<3, J, W>> for &TensorRank1List<3, I, W>
{
    type Output = TensorRank2<3, I, J>;
    fn mul(self, tensor_rank_1_list: TensorRank1List<3, J, W>) -> Self::Output
    {
        self.iter().zip(tensor_rank_1_list.iter()).map(|(self_entry, tensor_rank_1_list_entry)|
            TensorRank2::dyad(self_entry, tensor_rank_1_list_entry)
        ).sum()
    }
}

impl<const I: usize, const J: usize, const W: usize> Mul<&TensorRank1List<3, J, W>> for &TensorRank1List<3, I, W>
{
    type Output = TensorRank2<3, I, J>;
    fn mul(self, tensor_rank_1_list: &TensorRank1List<3, J, W>) -> Self::Output
    {
        self.iter().zip(tensor_rank_1_list.iter()).map(|(self_entry, tensor_rank_1_list_entry)|
            TensorRank2::dyad(self_entry, tensor_rank_1_list_entry)
        ).sum()
    }
}

impl<const D: usize, const I: usize, const W: usize> Sub for TensorRank1List<D, I, W>
{
    type Output = Self;
    fn sub(mut self, tensor_rank_1_list: Self) -> Self::Output
    {
        self -= tensor_rank_1_list;
        self
    }
}

impl<const D: usize, const I: usize, const W: usize> Sub<&Self> for TensorRank1List<D, I, W>
{
    type Output = Self;
    fn sub(mut self, tensor_rank_1_list: &Self) -> Self::Output
    {
        self -= tensor_rank_1_list;
        self
    }
}

impl<const D: usize, const I: usize, const W: usize> SubAssign for TensorRank1List<D, I, W>
{
    fn sub_assign(&mut self, tensor_rank_1_list: Self)
    {
        self.iter_mut().zip(tensor_rank_1_list.iter()).for_each(|(self_entry, tensor_rank_1_list_entry)|
            *self_entry -= tensor_rank_1_list_entry
        );
    }
}

impl<const D: usize, const I: usize, const W: usize> SubAssign<&Self> for TensorRank1List<D, I, W>
{
    fn sub_assign(&mut self, tensor_rank_1_list: &Self)
    {
        self.iter_mut().zip(tensor_rank_1_list.iter()).for_each(|(self_entry, tensor_rank_1_list_entry)|
            *self_entry -= tensor_rank_1_list_entry
        );
    }
}