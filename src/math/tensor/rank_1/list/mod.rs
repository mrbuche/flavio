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

pub struct TensorRank1List<const D: usize, const I: usize, const L: usize>
(
    [TensorRank1<D, I>; L]
);

impl<const D: usize, const I: usize, const L: usize> TensorRank1List<D, I, L>
{
    pub fn iter(&self) -> impl Iterator<Item = &TensorRank1<D, I>>
    {
        self.0.iter()
    }
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut TensorRank1<D, I>>
    {
        self.0.iter_mut()
    }
}

pub trait TensorRank1ListTrait<const D: usize, const L: usize>
{
    fn new(array: [[TensorRank0; D]; L]) -> Self;
    fn to_current_configuration(self) -> TensorRank1List<D, 1, L>;
    fn to_intermediate_configuration(self) -> TensorRank1List<D, 2, L>;
    fn to_reference_configuration(self) -> TensorRank1List<D, 0, L>;
    fn zero() -> Self;
}

impl<const D: usize, const I: usize, const L: usize> TensorRank1ListTrait<D, L> for TensorRank1List<D, I, L>
{
    fn new(array: [[TensorRank0; D]; L]) -> Self
    {
        array.iter().map(|array_i|
            TensorRank1::new(*array_i)
        ).collect()
    }
    fn to_current_configuration(self) -> TensorRank1List<D, 1, L>
    {
        TensorRank1List(
            from_fn(|entry|
                self.0[entry].to_current_configuration()
            )
        )
    }
    fn to_intermediate_configuration(self) -> TensorRank1List<D, 2, L>
    {
        TensorRank1List(
            from_fn(|entry|
                self.0[entry].to_intermediate_configuration()
            )
        )
    }
    fn to_reference_configuration(self) -> TensorRank1List<D, 0, L>
    {
        TensorRank1List(
            from_fn(|entry|
                self.0[entry].to_reference_configuration()
            )
        )
    }
    fn zero() -> Self
    {
        Self(from_fn(|_| TensorRank1::zero()))
    }
}

impl<const D: usize, const I: usize, const L: usize> FromIterator<TensorRank1<D, I>> for TensorRank1List<D, I, L>
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

impl<const D: usize, const I: usize, const L: usize> Index<usize> for TensorRank1List<D, I, L>
{
    type Output = TensorRank1<D, I>;
    fn index(&self, index: usize) -> &Self::Output
    {
        &self.0[index]
    }
}

impl<const D: usize, const I: usize, const L: usize> IndexMut<usize> for TensorRank1List<D, I, L>
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output
    {
        &mut self.0[index]
    }
}

impl<const D: usize, const I: usize, const L: usize> Div<TensorRank0> for TensorRank1List<D, I, L>
{
    type Output = Self;
    fn div(mut self, tensor_rank_0: TensorRank0) -> Self::Output
    {
        self /= &tensor_rank_0;
        self
    }
}

impl<const D: usize, const I: usize, const L: usize> Div<TensorRank0> for &TensorRank1List<D, I, L>
{
    type Output = TensorRank1List<D, I, L>;
    fn div(self, tensor_rank_0: TensorRank0) -> Self::Output
    {
        self.iter().map(|self_i|
            self_i / tensor_rank_0
        ).collect()
    }
}

impl<const D: usize, const I: usize, const L: usize> Div<&TensorRank0> for TensorRank1List<D, I, L>
{
    type Output = Self;
    fn div(mut self, tensor_rank_0: &TensorRank0) -> Self::Output
    {
        self /= tensor_rank_0;
        self
    }
}

impl<const D: usize, const I: usize, const L: usize> Div<&TensorRank0> for &TensorRank1List<D, I, L>
{
    type Output = TensorRank1List<D, I, L>;
    fn div(self, tensor_rank_0: &TensorRank0) -> Self::Output
    {
        self.iter().map(|self_i|
            self_i / tensor_rank_0
        ).collect()
    }
}

impl<const D: usize, const I: usize, const L: usize> DivAssign<TensorRank0> for TensorRank1List<D, I, L>
{
    fn div_assign(&mut self, tensor_rank_0: TensorRank0)
    {
        self.iter_mut().for_each(|entry|
            *entry /= &tensor_rank_0
        );
    }
}

impl<const D: usize, const I: usize, const L: usize> DivAssign<&TensorRank0> for TensorRank1List<D, I, L>
{
    fn div_assign(&mut self, tensor_rank_0: &TensorRank0)
    {
        self.iter_mut().for_each(|entry|
            *entry /= tensor_rank_0
        );
    }
}

impl<const D: usize, const I: usize, const L: usize> Mul<TensorRank0> for TensorRank1List<D, I, L>
{
    type Output = Self;
    fn mul(mut self, tensor_rank_0: TensorRank0) -> Self::Output
    {
        self *= &tensor_rank_0;
        self
    }
}

impl<const D: usize, const I: usize, const L: usize> Mul<TensorRank0> for &TensorRank1List<D, I, L>
{
    type Output = TensorRank1List<D, I, L>;
    fn mul(self, tensor_rank_0: TensorRank0) -> Self::Output
    {
        self.iter().map(|self_i|
            self_i * tensor_rank_0
        ).collect()
    }
}

impl<const D: usize, const I: usize, const L: usize> Mul<&TensorRank0> for TensorRank1List<D, I, L>
{
    type Output = Self;
    fn mul(mut self, tensor_rank_0: &TensorRank0) -> Self::Output
    {
        self *= tensor_rank_0;
        self
    }
}

impl<const D: usize, const I: usize, const L: usize> Mul<&TensorRank0> for &TensorRank1List<D, I, L>
{
    type Output = TensorRank1List<D, I, L>;
    fn mul(self, tensor_rank_0: &TensorRank0) -> Self::Output
    {
        self.iter().map(|self_i|
            self_i * tensor_rank_0
        ).collect()
    }
}

impl<const D: usize, const I: usize, const L: usize> MulAssign<TensorRank0> for TensorRank1List<D, I, L>
{
    fn mul_assign(&mut self, tensor_rank_0: TensorRank0)
    {
        self.iter_mut().for_each(|entry|
            *entry *= &tensor_rank_0
        );
    }
}

impl<const D: usize, const I: usize, const L: usize> MulAssign<&TensorRank0> for TensorRank1List<D, I, L>
{
    fn mul_assign(&mut self, tensor_rank_0: &TensorRank0)
    {
        self.iter_mut().for_each(|entry|
            *entry *= tensor_rank_0
        );
    }
}

impl<const D: usize, const I: usize, const L: usize> Add for TensorRank1List<D, I, L>
{
    type Output = Self;
    fn add(mut self, tensor_rank_1_list: Self) -> Self::Output
    {
        self += tensor_rank_1_list;
        self
    }
}

impl<const D: usize, const I: usize, const L: usize> Add<&Self> for TensorRank1List<D, I, L>
{
    type Output = Self;
    fn add(mut self, tensor_rank_1_list: &Self) -> Self::Output
    {
        self += tensor_rank_1_list;
        self
    }
}

impl<const D: usize, const I: usize, const L: usize> Add<TensorRank1List<D, I, L>> for &TensorRank1List<D, I, L>
{
    type Output = TensorRank1List<D, I, L>;
    fn add(self, mut tensor_rank_1_list: TensorRank1List<D, I, L>) -> Self::Output
    {
        tensor_rank_1_list += self;
        tensor_rank_1_list
    }
}

impl<const D: usize, const I: usize, const L: usize> AddAssign for TensorRank1List<D, I, L>
{
    fn add_assign(&mut self, tensor_rank_1_list: Self)
    {
        self.iter_mut().zip(tensor_rank_1_list.iter()).for_each(|(self_entry, tensor_rank_1_list_entry)|
            *self_entry += tensor_rank_1_list_entry
        );
    }
}

impl<const D: usize, const I: usize, const L: usize> AddAssign<&Self> for TensorRank1List<D, I, L>
{
    fn add_assign(&mut self, tensor_rank_1_list: &Self)
    {
        self.iter_mut().zip(tensor_rank_1_list.iter()).for_each(|(self_entry, tensor_rank_1_list_entry)|
            *self_entry += tensor_rank_1_list_entry
        );
    }
}

impl<const I: usize, const J: usize, const L: usize> Mul<TensorRank1List<3, J, L>> for TensorRank1List<3, I, L>
{
    type Output = TensorRank2<3, I, J>;
    fn mul(self, tensor_rank_1_list: TensorRank1List<3, J, L>) -> Self::Output
    {
        self.iter().zip(tensor_rank_1_list.iter()).map(|(self_entry, tensor_rank_1_list_entry)|
            TensorRank2::dyad(self_entry, tensor_rank_1_list_entry)
        ).sum()
    }
}

impl<const I: usize, const J: usize, const L: usize> Mul<&TensorRank1List<3, J, L>> for TensorRank1List<3, I, L>
{
    type Output = TensorRank2<3, I, J>;
    fn mul(self, tensor_rank_1_list: &TensorRank1List<3, J, L>) -> Self::Output
    {
        self.iter().zip(tensor_rank_1_list.iter()).map(|(self_entry, tensor_rank_1_list_entry)|
            TensorRank2::dyad(self_entry, tensor_rank_1_list_entry)
        ).sum()
    }
}

impl<const I: usize, const J: usize, const L: usize> Mul<TensorRank1List<3, J, L>> for &TensorRank1List<3, I, L>
{
    type Output = TensorRank2<3, I, J>;
    fn mul(self, tensor_rank_1_list: TensorRank1List<3, J, L>) -> Self::Output
    {
        self.iter().zip(tensor_rank_1_list.iter()).map(|(self_entry, tensor_rank_1_list_entry)|
            TensorRank2::dyad(self_entry, tensor_rank_1_list_entry)
        ).sum()
    }
}

impl<const I: usize, const J: usize, const L: usize> Mul<&TensorRank1List<3, J, L>> for &TensorRank1List<3, I, L>
{
    type Output = TensorRank2<3, I, J>;
    fn mul(self, tensor_rank_1_list: &TensorRank1List<3, J, L>) -> Self::Output
    {
        self.iter().zip(tensor_rank_1_list.iter()).map(|(self_entry, tensor_rank_1_list_entry)|
            TensorRank2::dyad(self_entry, tensor_rank_1_list_entry)
        ).sum()
    }
}

impl<const D: usize, const I: usize, const L: usize> Sub for TensorRank1List<D, I, L>
{
    type Output = Self;
    fn sub(mut self, tensor_rank_1_list: Self) -> Self::Output
    {
        self -= tensor_rank_1_list;
        self
    }
}

impl<const D: usize, const I: usize, const L: usize> Sub<&Self> for TensorRank1List<D, I, L>
{
    type Output = Self;
    fn sub(mut self, tensor_rank_1_list: &Self) -> Self::Output
    {
        self -= tensor_rank_1_list;
        self
    }
}

impl<const D: usize, const I: usize, const L: usize> SubAssign for TensorRank1List<D, I, L>
{
    fn sub_assign(&mut self, tensor_rank_1_list: Self)
    {
        self.iter_mut().zip(tensor_rank_1_list.iter()).for_each(|(self_entry, tensor_rank_1_list_entry)|
            *self_entry -= tensor_rank_1_list_entry
        );
    }
}

impl<const D: usize, const I: usize, const L: usize> SubAssign<&Self> for TensorRank1List<D, I, L>
{
    fn sub_assign(&mut self, tensor_rank_1_list: &Self)
    {
        self.iter_mut().zip(tensor_rank_1_list.iter()).for_each(|(self_entry, tensor_rank_1_list_entry)|
            *self_entry -= tensor_rank_1_list_entry
        );
    }
}