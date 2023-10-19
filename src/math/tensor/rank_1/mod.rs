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

use super::
{
    rank_0::TensorRank0,
    rank_2::TensorRank2
};

pub struct TensorRank1<const D: usize>
(
    pub [TensorRank0; D]
);

impl<const D: usize> TensorRank1<D>
{
    pub fn iter(&self) -> impl Iterator<Item = &TensorRank0>
    {
        self.0.iter()
    }
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut TensorRank0>
    {
        self.0.iter_mut()
    }
}

pub trait TensorRank1Traits<'a, const D: usize>
where
    Self: 'a
        + FromIterator<TensorRank0>
        + Mul<Output = TensorRank0>
        + Sized,
    &'a Self: Mul<&'a Self, Output = TensorRank0>
{
    fn new(array: [TensorRank0; D]) -> Self
    {
        array.into_iter().collect()
    }
    fn norm(&'a self) -> TensorRank0
    {
        (self * self).sqrt()
    }
    fn zero() -> Self
    {
        Self::new([0.0; D])
    }
}

impl<'a, const D: usize> TensorRank1Traits<'a, D> for TensorRank1<D> {}

impl<const D: usize> FromIterator<TensorRank0> for TensorRank1<D>
{
    fn from_iter<I: IntoIterator<Item=TensorRank0>>(into_iterator: I) -> Self
    {
        let mut tensor_rank_1 = Self([0.0; D]);
        tensor_rank_1.iter_mut().zip(into_iterator).for_each(|(tensor_rank_1_i, value_i)|
            *tensor_rank_1_i = value_i
        );
        tensor_rank_1
    }
}

impl<const D: usize> Index<usize> for TensorRank1<D>
{
    type Output = TensorRank0;
    fn index(&self, index: usize) -> &Self::Output
    {
        &self.0[index]
    }
}

impl<const D: usize> IndexMut<usize> for TensorRank1<D>
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output
    {
        &mut self.0[index]
    }
}

impl<const D: usize> std::iter::Sum for TensorRank1<D>
{
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = Self>
    {
        let mut output = TensorRank1::zero();
        iter.for_each(|item|
            output += item
        );
        output
    }
}

impl<const D: usize> Div<TensorRank0> for TensorRank1<D>
{
    type Output = Self;
    fn div(mut self, tensor_rank_0: TensorRank0) -> Self::Output
    {
        self /= tensor_rank_0;
        self
    }
}

impl<const D: usize> Div<&TensorRank0> for TensorRank1<D>
{
    type Output = Self;
    fn div(mut self, tensor_rank_0: &TensorRank0) -> Self::Output
    {
        self /= tensor_rank_0;
        self
    }
}

impl<const D: usize> DivAssign<TensorRank0> for TensorRank1<D>
{
    fn div_assign(&mut self, tensor_rank_0: TensorRank0)
    {
        self.iter_mut().for_each(|self_i|
            *self_i /= &tensor_rank_0
        );
    }
}

impl<const D: usize> DivAssign<&TensorRank0> for TensorRank1<D>
{
    fn div_assign(&mut self, tensor_rank_0: &TensorRank0)
    {
        self.iter_mut().for_each(|self_i|
            *self_i /= tensor_rank_0
        );
    }
}

impl<const D: usize> Mul<TensorRank0> for TensorRank1<D>
{
    type Output = Self;
    fn mul(mut self, tensor_rank_0: TensorRank0) -> Self::Output
    {
        self *= tensor_rank_0;
        self
    }
}

impl<const D: usize> Mul<TensorRank0> for &TensorRank1<D>
{
    type Output = TensorRank1<D>;
    fn mul(self, tensor_rank_0: TensorRank0) -> Self::Output
    {
        self.iter().map(|self_i|
            self_i * tensor_rank_0
        ).collect()
    }
}

impl<const D: usize> Mul<&TensorRank0> for TensorRank1<D>
{
    type Output = Self;
    fn mul(mut self, tensor_rank_0: &TensorRank0) -> Self::Output
    {
        self *= tensor_rank_0;
        self
    }
}

impl<const D: usize> Mul<&TensorRank0> for &TensorRank1<D>
{
    type Output = TensorRank1<D>;
    fn mul(self, tensor_rank_0: &TensorRank0) -> Self::Output
    {
        self.iter().map(|self_i|
            self_i * tensor_rank_0
        ).collect()
    }
}

impl<const D: usize> MulAssign<TensorRank0> for TensorRank1<D>
{
    fn mul_assign(&mut self, tensor_rank_0: TensorRank0)
    {
        self.iter_mut().for_each(|self_i|
            *self_i *= &tensor_rank_0
        );
    }
}

impl<const D: usize> MulAssign<&TensorRank0> for TensorRank1<D>
{
    fn mul_assign(&mut self, tensor_rank_0: &TensorRank0)
    {
        self.iter_mut().for_each(|self_i|
            *self_i *= tensor_rank_0
        );
    }
}

impl<const D: usize> Add for TensorRank1<D>
{
    type Output = Self;
    fn add(mut self, tensor_rank_1: Self) -> Self::Output
    {
        self += tensor_rank_1;
        self
    }
}

impl<const D: usize> Add<&Self> for TensorRank1<D>
{
    type Output = Self;
    fn add(mut self, tensor_rank_1: &Self) -> Self::Output
    {
        self += tensor_rank_1;
        self
    }
}

impl<const D: usize> Add<TensorRank1<D>> for &TensorRank1<D>
{
    type Output = TensorRank1<D>;
    fn add(self, mut tensor_rank_1: TensorRank1<D>) -> Self::Output
    {
        tensor_rank_1 += self;
        tensor_rank_1
    }
}

impl<const D: usize> AddAssign for TensorRank1<D>
{
    fn add_assign(&mut self, tensor_rank_1: Self)
    {
        self.iter_mut().zip(tensor_rank_1.iter()).for_each(|(self_i, tensor_rank_1_i)|
            *self_i += tensor_rank_1_i
        );
    }
}

impl<const D: usize> AddAssign<&Self> for TensorRank1<D>
{
    fn add_assign(&mut self, tensor_rank_1: &Self)
    {
        self.iter_mut().zip(tensor_rank_1.iter()).for_each(|(self_i, tensor_rank_1_i)|
            *self_i += tensor_rank_1_i
        );
    }
}

impl<const D: usize> Sub for TensorRank1<D>
{
    type Output = Self;
    fn sub(mut self, tensor_rank_1: Self) -> Self::Output
    {
        self -= tensor_rank_1;
        self
    }
}

impl<const D: usize> Sub<&Self> for TensorRank1<D>
{
    type Output = Self;
    fn sub(mut self, tensor_rank_1: &Self) -> Self::Output
    {
        self -= tensor_rank_1;
        self
    }
}

impl<const D: usize> Sub<TensorRank1<D>> for &TensorRank1<D>
{
    type Output = TensorRank1<D>;
    fn sub(self, mut tensor_rank_1: TensorRank1<D>) -> Self::Output
    {
        tensor_rank_1.iter_mut().zip(self.iter()).for_each(|(tensor_rank_1_i, self_i)|
            *tensor_rank_1_i = self_i - *tensor_rank_1_i
        );
        tensor_rank_1
    }
}

impl<const D: usize> SubAssign for TensorRank1<D>
{
    fn sub_assign(&mut self, tensor_rank_1: Self)
    {
        self.iter_mut().zip(tensor_rank_1.iter()).for_each(|(self_i, tensor_rank_1_i)|
            *self_i -= tensor_rank_1_i
        );
    }
}

impl<const D: usize> SubAssign<&Self> for TensorRank1<D>
{
    fn sub_assign(&mut self, tensor_rank_1: &Self)
    {
        self.iter_mut().zip(tensor_rank_1.iter()).for_each(|(self_i, tensor_rank_1_i)|
            *self_i -= tensor_rank_1_i
        );
    }
}

impl<const D: usize> Mul for TensorRank1<D>
{
    type Output = TensorRank0;
    fn mul(self, tensor_rank_1: Self) -> Self::Output
    {
        self.iter().zip(tensor_rank_1.iter()).map(|(self_i, tensor_rank_1_i)|
            self_i * tensor_rank_1_i
        ).sum()
    }
}

impl<const D: usize> Mul<&Self> for TensorRank1<D>
{
    type Output = TensorRank0;
    fn mul(self, tensor_rank_1: &Self) -> Self::Output
    {
        self.iter().zip(tensor_rank_1.iter()).map(|(self_i, tensor_rank_1_i)|
            self_i * tensor_rank_1_i
        ).sum()
    }
}

impl<const D: usize> Mul<TensorRank1<D>> for &TensorRank1<D>
{
    type Output = TensorRank0;
    fn mul(self, tensor_rank_1: TensorRank1<D>) -> Self::Output
    {
        self.iter().zip(tensor_rank_1.iter()).map(|(self_i, tensor_rank_1_i)|
            self_i * tensor_rank_1_i
        ).sum()
    }
}

impl<const D: usize> Mul for &TensorRank1<D>
{
    type Output = TensorRank0;
    fn mul(self, tensor_rank_1: Self) -> Self::Output
    {
        self.iter().zip(tensor_rank_1.iter()).map(|(self_i, tensor_rank_1_i)|
            self_i * tensor_rank_1_i
        ).sum()
    }
}
