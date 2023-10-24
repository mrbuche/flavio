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

use crate::math::
{
    TensorRank2,
    TensorRank2Traits
};

use super::
{
    TensorRank0,
    TensorRank1,
    TensorRank1Traits
};

pub struct TensorRank1List<const D: usize, const L: usize>
(
    pub [TensorRank1<D>; L]
);

impl<const D: usize, const L: usize> TensorRank1List<D, L>
{
    pub fn iter(&self) -> impl Iterator<Item = &TensorRank1<D>>
    {
        self.0.iter()
    }
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut TensorRank1<D>>
    {
        self.0.iter_mut()
    }
}

pub trait TensorRank1ListTraits<const D: usize, const L: usize>
where
    Self: FromIterator<TensorRank1<D>>
        + Index<usize, Output = TensorRank1<D>>
        + IndexMut<usize, Output = TensorRank1<D>>
        + Sized
{
    fn new(array: [[TensorRank0; D]; L]) -> Self
    {
        array.iter().map(|array_i|
            TensorRank1(*array_i)
        ).collect()
    }
    fn zero() -> Self;
}

impl<const D: usize, const L: usize> TensorRank1ListTraits<D, L> for TensorRank1List<D, L>
{
    fn zero() -> Self
    {
        Self(std::array::from_fn(|_| TensorRank1::zero()))
    }
}

impl<const D: usize, const L: usize> FromIterator<TensorRank1<D>> for TensorRank1List<D, L>
{
    fn from_iter<I: IntoIterator<Item=TensorRank1<D>>>(into_iterator: I) -> Self
    {
        let mut tensor_rank_1_list = Self::zero();
        tensor_rank_1_list.iter_mut().zip(into_iterator).for_each(|(tensor_rank_1_list_entry, entry)|
            *tensor_rank_1_list_entry = entry
        );
        tensor_rank_1_list
    }
}

impl<const D: usize, const L: usize> Index<usize> for TensorRank1List<D, L>
{
    type Output = TensorRank1<D>;
    fn index(&self, index: usize) -> &Self::Output
    {
        &self.0[index]
    }
}

impl<const D: usize, const L: usize> IndexMut<usize> for TensorRank1List<D, L>
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output
    {
        &mut self.0[index]
    }
}

impl<const D: usize, const L: usize> Div<TensorRank0> for TensorRank1List<D, L>
{
    type Output = Self;
    fn div(mut self, tensor_rank_0: TensorRank0) -> Self::Output
    {
        self /= &tensor_rank_0;
        self
    }
}

impl<const D: usize, const L: usize> Div<TensorRank0> for &TensorRank1List<D, L>
{
    type Output = TensorRank1List<D, L>;
    fn div(self, tensor_rank_0: TensorRank0) -> Self::Output
    {
        self.iter().map(|self_i|
            self_i / tensor_rank_0
        ).collect()
    }
}

impl<const D: usize, const L: usize> Div<&TensorRank0> for TensorRank1List<D, L>
{
    type Output = Self;
    fn div(mut self, tensor_rank_0: &TensorRank0) -> Self::Output
    {
        self /= tensor_rank_0;
        self
    }
}

impl<const D: usize, const L: usize> Div<&TensorRank0> for &TensorRank1List<D, L>
{
    type Output = TensorRank1List<D, L>;
    fn div(self, tensor_rank_0: &TensorRank0) -> Self::Output
    {
        self.iter().map(|self_i|
            self_i / tensor_rank_0
        ).collect()
    }
}

impl<const D: usize, const L: usize> DivAssign<TensorRank0> for TensorRank1List<D, L>
{
    fn div_assign(&mut self, tensor_rank_0: TensorRank0)
    {
        self.iter_mut().for_each(|entry|
            *entry /= &tensor_rank_0
        );
    }
}

impl<const D: usize, const L: usize> DivAssign<&TensorRank0> for TensorRank1List<D, L>
{
    fn div_assign(&mut self, tensor_rank_0: &TensorRank0)
    {
        self.iter_mut().for_each(|entry|
            *entry /= tensor_rank_0
        );
    }
}

impl<const D: usize, const L: usize> Mul<TensorRank0> for TensorRank1List<D, L>
{
    type Output = Self;
    fn mul(mut self, tensor_rank_0: TensorRank0) -> Self::Output
    {
        self *= &tensor_rank_0;
        self
    }
}

impl<const D: usize, const L: usize> Mul<TensorRank0> for &TensorRank1List<D, L>
{
    type Output = TensorRank1List<D, L>;
    fn mul(self, tensor_rank_0: TensorRank0) -> Self::Output
    {
        self.iter().map(|self_i|
            self_i * tensor_rank_0
        ).collect()
    }
}

impl<const D: usize, const L: usize> Mul<&TensorRank0> for TensorRank1List<D, L>
{
    type Output = Self;
    fn mul(mut self, tensor_rank_0: &TensorRank0) -> Self::Output
    {
        self *= tensor_rank_0;
        self
    }
}

impl<const D: usize, const L: usize> Mul<&TensorRank0> for &TensorRank1List<D, L>
{
    type Output = TensorRank1List<D, L>;
    fn mul(self, tensor_rank_0: &TensorRank0) -> Self::Output
    {
        self.iter().map(|self_i|
            self_i * tensor_rank_0
        ).collect()
    }
}

impl<const D: usize, const L: usize> MulAssign<TensorRank0> for TensorRank1List<D, L>
{
    fn mul_assign(&mut self, tensor_rank_0: TensorRank0)
    {
        self.iter_mut().for_each(|entry|
            *entry *= &tensor_rank_0
        );
    }
}

impl<const D: usize, const L: usize> MulAssign<&TensorRank0> for TensorRank1List<D, L>
{
    fn mul_assign(&mut self, tensor_rank_0: &TensorRank0)
    {
        self.iter_mut().for_each(|entry|
            *entry *= tensor_rank_0
        );
    }
}

impl<const L: usize> Mul for TensorRank1List<3, L>
{
    type Output = TensorRank2<3>;
    fn mul(self, tensor_rank_1_list: Self) -> Self::Output
    {
        self.iter().zip(tensor_rank_1_list.iter()).map(|(self_entry, tensor_rank_1_list_entry)|
            TensorRank2::dyad(self_entry, tensor_rank_1_list_entry)
        ).sum()
    }
}

impl<const L: usize> Mul<&Self> for TensorRank1List<3, L>
{
    type Output = TensorRank2<3>;
    fn mul(self, tensor_rank_1_list: &Self) -> Self::Output
    {
        self.iter().zip(tensor_rank_1_list.iter()).map(|(self_entry, tensor_rank_1_list_entry)|
            TensorRank2::dyad(self_entry, tensor_rank_1_list_entry)
        ).sum()
    }
}

impl<const L: usize> Mul<TensorRank1List<3, L>> for &TensorRank1List<3, L>
{
    type Output = TensorRank2<3>;
    fn mul(self, tensor_rank_1_list: TensorRank1List<3, L>) -> Self::Output
    {
        self.iter().zip(tensor_rank_1_list.iter()).map(|(self_entry, tensor_rank_1_list_entry)|
            TensorRank2::dyad(self_entry, tensor_rank_1_list_entry)
        ).sum()
    }
}

impl<const L: usize> Mul<Self> for &TensorRank1List<3, L>
{
    type Output = TensorRank2<3>;
    fn mul(self, tensor_rank_1_list: Self) -> Self::Output
    {
        self.iter().zip(tensor_rank_1_list.iter()).map(|(self_entry, tensor_rank_1_list_entry)|
            TensorRank2::dyad(self_entry, tensor_rank_1_list_entry)
        ).sum()
    }
}

impl<const D: usize, const L: usize> Add for TensorRank1List<D, L>
{
    type Output = Self;
    fn add(mut self, tensor_rank_1_list: Self) -> Self::Output
    {
        self += tensor_rank_1_list;
        self
    }
}

impl<const D: usize, const L: usize> Add<&Self> for TensorRank1List<D, L>
{
    type Output = Self;
    fn add(mut self, tensor_rank_1_list: &Self) -> Self::Output
    {
        self += tensor_rank_1_list;
        self
    }
}

impl<const D: usize, const L: usize> Add<TensorRank1List<D, L>> for &TensorRank1List<D, L>
{
    type Output = TensorRank1List<D, L>;
    fn add(self, mut tensor_rank_1_list: TensorRank1List<D, L>) -> Self::Output
    {
        tensor_rank_1_list += self;
        tensor_rank_1_list
    }
}

impl<const D: usize, const L: usize> AddAssign for TensorRank1List<D, L>
{
    fn add_assign(&mut self, tensor_rank_1_list: Self)
    {
        self.iter_mut().zip(tensor_rank_1_list.iter()).for_each(|(self_entry, tensor_rank_1_list_entry)|
            *self_entry += tensor_rank_1_list_entry
        );
    }
}

impl<const D: usize, const L: usize> AddAssign<&Self> for TensorRank1List<D, L>
{
    fn add_assign(&mut self, tensor_rank_1_list: &Self)
    {
        self.iter_mut().zip(tensor_rank_1_list.iter()).for_each(|(self_entry, tensor_rank_1_list_entry)|
            *self_entry += tensor_rank_1_list_entry
        );
    }
}
impl<const D: usize, const L: usize> Sub for TensorRank1List<D, L>
{
    type Output = Self;
    fn sub(mut self, tensor_rank_1_list: Self) -> Self::Output
    {
        self -= tensor_rank_1_list;
        self
    }
}

impl<const D: usize, const L: usize> Sub<&Self> for TensorRank1List<D, L>
{
    type Output = Self;
    fn sub(mut self, tensor_rank_1_list: &Self) -> Self::Output
    {
        self -= tensor_rank_1_list;
        self
    }
}

impl<const D: usize, const L: usize> SubAssign for TensorRank1List<D, L>
{
    fn sub_assign(&mut self, tensor_rank_1_list: Self)
    {
        self.iter_mut().zip(tensor_rank_1_list.iter()).for_each(|(self_entry, tensor_rank_1_list_entry)|
            *self_entry -= tensor_rank_1_list_entry
        );
    }
}

impl<const D: usize, const L: usize> SubAssign<&Self> for TensorRank1List<D, L>
{
    fn sub_assign(&mut self, tensor_rank_1_list: &Self)
    {
        self.iter_mut().zip(tensor_rank_1_list.iter()).for_each(|(self_entry, tensor_rank_1_list_entry)|
            *self_entry -= tensor_rank_1_list_entry
        );
    }
}