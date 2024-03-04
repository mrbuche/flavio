#[cfg(test)]
mod test;

pub mod list;

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

use super::
{
    Convert,
    rank_0::TensorRank0
};

/// A *d*-dimensional tensor of rank 1.
///
/// `D` is the dimension, `I` is the configuration.
pub struct TensorRank1<const D: usize, const I: usize>
(
    [TensorRank0; D]
);

/// Inherent implementation of [`TensorRank1`].
impl<const D: usize, const I: usize> TensorRank1<D, I>
{
    /// Returns an iterator.
    ///
    /// The iterator yields all items from start to end. [Read more](https://doc.rust-lang.org/std/iter/)
    pub fn iter(&self) -> impl Iterator<Item = &TensorRank0>
    {
        self.0.iter()
    }
    /// Returns an iterator that allows modifying each value.
    ///
    /// The iterator yields all items from start to end. [Read more](https://doc.rust-lang.org/std/iter/)
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut TensorRank0>
    {
        self.0.iter_mut()
    }
}

/// Required methods for rank-1 tensors.
pub trait TensorRank1Trait<const D: usize, const I: usize>
{
    /// Returns the rank-1 tensor as an array.
    fn as_array(&self) -> [TensorRank0; D];
    /// Returns the cross product with another rank-1 tensor.
    fn cross(&self, tensor_rank_1: &Self) -> Self;
    /// Returns a rank-1 tensor given an array.
    fn new(array: [TensorRank0; D]) -> Self;
    /// Returns the rank-1 tensor norm.
    fn norm(&self) -> TensorRank0;
    /// Returns the rank-1 tensor normalized.
    fn normalized(&self) -> Self;
    /// Returns the rank-1 tensor in the current configuration.
    fn to_current_configuration(self) -> TensorRank1<D, 1>;
    /// Returns the rank-1 tensor in the intermediate configuration.
    fn to_intermediate_configuration(self) -> TensorRank1<D, 2>;
    /// Returns the rank-1 tensor in the reference configuration.
    fn to_reference_configuration(self) -> TensorRank1<D, 0>;
    /// Returns the rank-1 zero tensor.
    fn zero() -> Self;
}

/// Implementation of [`TensorRank1Trait`] for [`TensorRank1`].
impl<const D: usize, const I: usize> TensorRank1Trait<D, I> for TensorRank1<D, I>
{
    fn as_array(&self) -> [TensorRank0; D]
    {
        self.0
    }
    fn cross(&self, tensor_rank_1: &Self) -> Self
    {
        if D == 3
        {
            let mut output = TensorRank1::<D, I>::zero();
            output[0] = self[1]*tensor_rank_1[2] - self[2]*tensor_rank_1[1];
            output[1] = self[2]*tensor_rank_1[0] - self[0]*tensor_rank_1[2];
            output[2] = self[0]*tensor_rank_1[1] - self[1]*tensor_rank_1[0];
            output
        }
        else
        {
            panic!()
        }
    }
    fn new(array: [TensorRank0; D]) -> Self
    {
        array.into_iter().collect()
    }
    fn norm(&self) -> TensorRank0
    {
        (self * self).sqrt()
    }
    fn normalized(&self) -> Self
    {
        self / self.norm()
    }
    fn to_current_configuration(self) -> TensorRank1<D, 1>
    {
        TensorRank1(self.0)
    }
    fn to_intermediate_configuration(self) -> TensorRank1<D, 2>
    {
        TensorRank1(self.0)
    }
    fn to_reference_configuration(self) -> TensorRank1<D, 0>
    {
        TensorRank1(self.0)
    }
    fn zero() -> Self
    {
        Self([0.0; D])
    }
}

impl<const D: usize, const I: usize, const J: usize> Convert<TensorRank1<D, J>> for TensorRank1<D, I>
{
    fn convert(&self) -> TensorRank1<D, J>
    {
        self.iter().map(|self_i|
            *self_i
        ).collect()
    }
}

impl<const D: usize, const I: usize> FromIterator<TensorRank0> for TensorRank1<D, I>
{
    fn from_iter<Ii: IntoIterator<Item=TensorRank0>>(into_iterator: Ii) -> Self
    {
        let mut tensor_rank_1 = Self::zero();
        tensor_rank_1.iter_mut().zip(into_iterator).for_each(|(tensor_rank_1_i, value_i)|
            *tensor_rank_1_i = value_i
        );
        tensor_rank_1
    }
}

impl<const D: usize, const I: usize> Index<usize> for TensorRank1<D, I>
{
    type Output = TensorRank0;
    fn index(&self, index: usize) -> &Self::Output
    {
        &self.0[index]
    }
}

impl<const D: usize, const I: usize> IndexMut<usize> for TensorRank1<D, I>
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output
    {
        &mut self.0[index]
    }
}

impl<const D: usize, const I: usize> std::iter::Sum for TensorRank1<D, I>
{
    fn sum<Ii>(iter: Ii) -> Self
    where
        Ii: Iterator<Item = Self>
    {
        let mut output = TensorRank1::zero();
        iter.for_each(|item|
            output += item
        );
        output
    }
}

impl<const D: usize, const I: usize> Div<TensorRank0> for TensorRank1<D, I>
{
    type Output = Self;
    fn div(mut self, tensor_rank_0: TensorRank0) -> Self::Output
    {
        self /= tensor_rank_0;
        self
    }
}

impl<const D: usize, const I: usize> Div<TensorRank0> for &TensorRank1<D, I>
{
    type Output = TensorRank1<D, I>;
    fn div(self, tensor_rank_0: TensorRank0) -> Self::Output
    {
        self.iter().map(|self_i|
            self_i / tensor_rank_0
        ).collect()
    }
}

impl<const D: usize, const I: usize> Div<&TensorRank0> for TensorRank1<D, I>
{
    type Output = Self;
    fn div(mut self, tensor_rank_0: &TensorRank0) -> Self::Output
    {
        self /= tensor_rank_0;
        self
    }
}

impl<const D: usize, const I: usize> Div<&TensorRank0> for &TensorRank1<D, I>
{
    type Output = TensorRank1<D, I>;
    fn div(self, tensor_rank_0: &TensorRank0) -> Self::Output
    {
        self.iter().map(|self_i|
            self_i / tensor_rank_0
        ).collect()
    }
}

impl<const D: usize, const I: usize> DivAssign<TensorRank0> for TensorRank1<D, I>
{
    fn div_assign(&mut self, tensor_rank_0: TensorRank0)
    {
        self.iter_mut().for_each(|self_i|
            *self_i /= &tensor_rank_0
        );
    }
}

impl<const D: usize, const I: usize> DivAssign<&TensorRank0> for TensorRank1<D, I>
{
    fn div_assign(&mut self, tensor_rank_0: &TensorRank0)
    {
        self.iter_mut().for_each(|self_i|
            *self_i /= tensor_rank_0
        );
    }
}

impl<const D: usize, const I: usize> Mul<TensorRank0> for TensorRank1<D, I>
{
    type Output = Self;
    fn mul(mut self, tensor_rank_0: TensorRank0) -> Self::Output
    {
        self *= tensor_rank_0;
        self
    }
}

impl<const D: usize, const I: usize> Mul<TensorRank0> for &TensorRank1<D, I>
{
    type Output = TensorRank1<D, I>;
    fn mul(self, tensor_rank_0: TensorRank0) -> Self::Output
    {
        self.iter().map(|self_i|
            self_i * tensor_rank_0
        ).collect()
    }
}

impl<const D: usize, const I: usize> Mul<&TensorRank0> for TensorRank1<D, I>
{
    type Output = Self;
    fn mul(mut self, tensor_rank_0: &TensorRank0) -> Self::Output
    {
        self *= tensor_rank_0;
        self
    }
}

impl<const D: usize, const I: usize> Mul<&TensorRank0> for &TensorRank1<D, I>
{
    type Output = TensorRank1<D, I>;
    fn mul(self, tensor_rank_0: &TensorRank0) -> Self::Output
    {
        self.iter().map(|self_i|
            self_i * tensor_rank_0
        ).collect()
    }
}

impl<const D: usize, const I: usize> MulAssign<TensorRank0> for TensorRank1<D, I>
{
    fn mul_assign(&mut self, tensor_rank_0: TensorRank0)
    {
        self.iter_mut().for_each(|self_i|
            *self_i *= &tensor_rank_0
        );
    }
}

impl<const D: usize, const I: usize> MulAssign<&TensorRank0> for TensorRank1<D, I>
{
    fn mul_assign(&mut self, tensor_rank_0: &TensorRank0)
    {
        self.iter_mut().for_each(|self_i|
            *self_i *= tensor_rank_0
        );
    }
}

impl<const D: usize, const I: usize> Add for TensorRank1<D, I>
{
    type Output = Self;
    fn add(mut self, tensor_rank_1: Self) -> Self::Output
    {
        self += tensor_rank_1;
        self
    }
}

impl<const D: usize, const I: usize> Add<&Self> for TensorRank1<D, I>
{
    type Output = Self;
    fn add(mut self, tensor_rank_1: &Self) -> Self::Output
    {
        self += tensor_rank_1;
        self
    }
}

impl<const D: usize, const I: usize> Add<TensorRank1<D, I>> for &TensorRank1<D, I>
{
    type Output = TensorRank1<D, I>;
    fn add(self, mut tensor_rank_1: TensorRank1<D, I>) -> Self::Output
    {
        tensor_rank_1 += self;
        tensor_rank_1
    }
}

impl<const D: usize, const I: usize> AddAssign for TensorRank1<D, I>
{
    fn add_assign(&mut self, tensor_rank_1: Self)
    {
        self.iter_mut().zip(tensor_rank_1.iter()).for_each(|(self_i, tensor_rank_1_i)|
            *self_i += tensor_rank_1_i
        );
    }
}

impl<const D: usize, const I: usize> AddAssign<&Self> for TensorRank1<D, I>
{
    fn add_assign(&mut self, tensor_rank_1: &Self)
    {
        self.iter_mut().zip(tensor_rank_1.iter()).for_each(|(self_i, tensor_rank_1_i)|
            *self_i += tensor_rank_1_i
        );
    }
}

impl<const D: usize, const I: usize> Sub for TensorRank1<D, I>
{
    type Output = Self;
    fn sub(mut self, tensor_rank_1: Self) -> Self::Output
    {
        self -= tensor_rank_1;
        self
    }
}

impl<const D: usize, const I: usize> Sub<&Self> for TensorRank1<D, I>
{
    type Output = Self;
    fn sub(mut self, tensor_rank_1: &Self) -> Self::Output
    {
        self -= tensor_rank_1;
        self
    }
}

impl<const D: usize, const I: usize> Sub<TensorRank1<D, I>> for &TensorRank1<D, I>
{
    type Output = TensorRank1<D, I>;
    fn sub(self, mut tensor_rank_1: TensorRank1<D, I>) -> Self::Output
    {
        tensor_rank_1.iter_mut().zip(self.iter()).for_each(|(tensor_rank_1_i, self_i)|
            *tensor_rank_1_i = self_i - *tensor_rank_1_i
        );
        tensor_rank_1
    }
}

impl<const D: usize, const I: usize> SubAssign for TensorRank1<D, I>
{
    fn sub_assign(&mut self, tensor_rank_1: Self)
    {
        self.iter_mut().zip(tensor_rank_1.iter()).for_each(|(self_i, tensor_rank_1_i)|
            *self_i -= tensor_rank_1_i
        );
    }
}

impl<const D: usize, const I: usize> SubAssign<&Self> for TensorRank1<D, I>
{
    fn sub_assign(&mut self, tensor_rank_1: &Self)
    {
        self.iter_mut().zip(tensor_rank_1.iter()).for_each(|(self_i, tensor_rank_1_i)|
            *self_i -= tensor_rank_1_i
        );
    }
}

impl<const D: usize, const I: usize> Mul for TensorRank1<D, I>
{
    type Output = TensorRank0;
    fn mul(self, tensor_rank_1: Self) -> Self::Output
    {
        self.iter().zip(tensor_rank_1.iter()).map(|(self_i, tensor_rank_1_i)|
            self_i * tensor_rank_1_i
        ).sum()
    }
}

impl<const D: usize, const I: usize> Mul<&Self> for TensorRank1<D, I>
{
    type Output = TensorRank0;
    fn mul(self, tensor_rank_1: &Self) -> Self::Output
    {
        self.iter().zip(tensor_rank_1.iter()).map(|(self_i, tensor_rank_1_i)|
            self_i * tensor_rank_1_i
        ).sum()
    }
}

impl<const D: usize, const I: usize> Mul<TensorRank1<D, I>> for &TensorRank1<D, I>
{
    type Output = TensorRank0;
    fn mul(self, tensor_rank_1: TensorRank1<D, I>) -> Self::Output
    {
        self.iter().zip(tensor_rank_1.iter()).map(|(self_i, tensor_rank_1_i)|
            self_i * tensor_rank_1_i
        ).sum()
    }
}

impl<const D: usize, const I: usize> Mul for &TensorRank1<D, I>
{
    type Output = TensorRank0;
    fn mul(self, tensor_rank_1: Self) -> Self::Output
    {
        self.iter().zip(tensor_rank_1.iter()).map(|(self_i, tensor_rank_1_i)|
            self_i * tensor_rank_1_i
        ).sum()
    }
}
