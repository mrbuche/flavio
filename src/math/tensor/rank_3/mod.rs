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
    rank_1::
    {
        TensorRank1,
        TensorRank1Trait
    },
    rank_2::TensorRank2
};

pub struct TensorRank3<const D: usize>
(
    [TensorRank2<D>; D]
);

impl<const D: usize> TensorRank3<D>
{
    pub fn iter(&self) -> impl Iterator<Item=&TensorRank2<D>>
    {
        self.0.iter()
    }
    pub fn iter_mut(&mut self) -> impl Iterator<Item=&mut TensorRank2<D>>
    {
        self.0.iter_mut()
    }
    pub fn zero() -> Self
    {
        Self(std::array::from_fn(|_| TensorRank2::zero()))
    }
}

pub trait TensorRank3Trait<const D: usize>
where
    Self: FromIterator<TensorRank2<D>>
        + Index<usize, Output = TensorRank2<D>>
        + IndexMut<usize, Output = TensorRank2<D>>
        + Sized
{
    fn new(array: [[[TensorRank0; D]; D]; D]) -> Self
    {
        array.iter().map(|array_i|
            array_i.iter().map(|array_ij|
                TensorRank1::new(*array_ij)
            ).collect()
        ).collect()
    }
}

impl<const D: usize> TensorRank3Trait<D> for TensorRank3<D> {}

impl<const D: usize> FromIterator<TensorRank2<D>> for TensorRank3<D>
{
    fn from_iter<I: IntoIterator<Item=TensorRank2<D>>>(into_iterator: I) -> Self
    {
        let mut tensor_rank_3 = Self::zero();
        tensor_rank_3.iter_mut().zip(into_iterator).for_each(|(tensor_rank_3_i, value_i)|
            *tensor_rank_3_i = value_i
        );
        tensor_rank_3
    }
}

impl<const D: usize> Index<usize> for TensorRank3<D>
{
    type Output = TensorRank2<D>;
    fn index(&self, index: usize) -> &Self::Output
    {
        &self.0[index]
    }
}

impl<const D: usize> IndexMut<usize> for TensorRank3<D>
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output
    {
        &mut self.0[index]
    }
}

impl<const D: usize> Div<TensorRank0> for TensorRank3<D>
{
    type Output = Self;
    fn div(mut self, tensor_rank_0: TensorRank0) -> Self::Output
    {
        self /= &tensor_rank_0;
        self
    }
}

impl<const D: usize> Div<&TensorRank0> for TensorRank3<D>
{
    type Output = Self;
    fn div(mut self, tensor_rank_0: &TensorRank0) -> Self::Output
    {
        self /= tensor_rank_0;
        self
    }
}

impl<const D: usize> DivAssign<TensorRank0> for TensorRank3<D>
{
    fn div_assign(&mut self, tensor_rank_0: TensorRank0)
    {
        self.iter_mut().for_each(|self_i|
            *self_i /= &tensor_rank_0
        );
    }
}

impl<const D: usize> DivAssign<&TensorRank0> for TensorRank3<D>
{
    fn div_assign(&mut self, tensor_rank_0: &TensorRank0)
    {
        self.iter_mut().for_each(|self_i|
            *self_i /= tensor_rank_0
        );
    }
}

impl<const D: usize> Mul<TensorRank0> for TensorRank3<D>
{
    type Output = Self;
    fn mul(mut self, tensor_rank_0: TensorRank0) -> Self::Output
    {
        self *= &tensor_rank_0;
        self
    }
}

impl<const D: usize> Mul<&TensorRank0> for TensorRank3<D>
{
    type Output = Self;
    fn mul(mut self, tensor_rank_0: &TensorRank0) -> Self::Output
    {
        self *= tensor_rank_0;
        self
    }
}

impl<const D: usize> MulAssign<TensorRank0> for TensorRank3<D>
{
    fn mul_assign(&mut self, tensor_rank_0: TensorRank0)
    {
        self.iter_mut().for_each(|self_i|
            *self_i *= &tensor_rank_0
        );
    }
}

impl<const D: usize> MulAssign<&TensorRank0> for TensorRank3<D>
{
    fn mul_assign(&mut self, tensor_rank_0: &TensorRank0)
    {
        self.iter_mut().for_each(|self_i|
            *self_i *= tensor_rank_0
        );
    }
}

impl<const D: usize> Add for TensorRank3<D>
{
    type Output = Self;
    fn add(mut self, tensor_rank_3: Self) -> Self::Output
    {
        self += tensor_rank_3;
        self
    }
}

impl<const D: usize> Add<&Self> for TensorRank3<D>
{
    type Output = Self;
    fn add(mut self, tensor_rank_3: &Self) -> Self::Output
    {
        self += tensor_rank_3;
        self
    }
}

impl<const D: usize> Add<TensorRank3<D>> for &TensorRank3<D>
{
    type Output = TensorRank3<D>;
    fn add(self, mut tensor_rank_3: TensorRank3<D>) -> Self::Output
    {
        tensor_rank_3 += self;
        tensor_rank_3
    }
}

impl<const D: usize> AddAssign for TensorRank3<D>
{
    fn add_assign(&mut self, tensor_rank_3: Self)
    {
        self.iter_mut().zip(tensor_rank_3.iter()).for_each(|(self_i, tensor_rank_3_i)|
            *self_i += tensor_rank_3_i
        );
    }
}

impl<const D: usize> AddAssign<&Self> for TensorRank3<D>
{
    fn add_assign(&mut self, tensor_rank_3: &Self)
    {
        self.iter_mut().zip(tensor_rank_3.iter()).for_each(|(self_i, tensor_rank_3_i)|
            *self_i += tensor_rank_3_i
        );
    }
}

impl<const D: usize> Sub for TensorRank3<D>
{
    type Output = Self;
    fn sub(mut self, tensor_rank_3: Self) -> Self::Output
    {
        self -= tensor_rank_3;
        self
    }
}

impl<const D: usize> Sub<&Self> for TensorRank3<D>
{
    type Output = Self;
    fn sub(mut self, tensor_rank_3: &Self) -> Self::Output
    {
        self -= tensor_rank_3;
        self
    }
}

impl<const D: usize> SubAssign for TensorRank3<D>
{
    fn sub_assign(&mut self, tensor_rank_3: Self)
    {
        self.iter_mut().zip(tensor_rank_3.iter()).for_each(|(self_i, tensor_rank_3_i)|
            *self_i -= tensor_rank_3_i
        );
    }
}

impl<const D: usize> SubAssign<&Self> for TensorRank3<D>
{
    fn sub_assign(&mut self, tensor_rank_3: &Self)
    {
        self.iter_mut().zip(tensor_rank_3.iter()).for_each(|(self_i, tensor_rank_3_i)|
            *self_i -= tensor_rank_3_i
        );
    }
}