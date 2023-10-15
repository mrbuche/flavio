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
        TensorRank1Traits
    }
};

// eliminate in order to identify explicit copying at some point
#[derive(Clone, Copy)]
pub struct TensorRank2<const D: usize>
(
    pub [TensorRank1<D>; D]
);

// move into TensorRank2Traits if ever becomes possible
impl<const D: usize> TensorRank2<D>
{
    fn iter(&self) -> impl Iterator<Item=&TensorRank1<D>>
    {
        self.0.iter()
    }
    fn iter_mut(&mut self) -> impl Iterator<Item=&mut TensorRank1<D>>
    {
        self.0.iter_mut()
    }
}

pub trait TensorRank2Traits<const D: usize>
where
    Self: Index<usize, Output = TensorRank1<D>>
{
    type Inverse;
    fn inverse(&self) -> Self::Inverse
    {
        todo!("Can use the LU decomposition as the default implementation and override for certain dimensions.");
    }
    fn trace(&self) -> TensorRank0
    {
        (0..D).map(|i| self[i][i]).sum()
    }
    fn transpose(&self) -> Self::Inverse
    {
        todo!();
    }
    fn zero() -> Self;
}

impl<const D: usize> TensorRank2Traits<D> for TensorRank2<D>
{
    type Inverse = Self;
    fn zero() -> Self
    {
        Self(std::array::from_fn(|_| TensorRank1::zero()))
    }
}

impl<const D: usize> Index<usize> for TensorRank2<D>
{
    type Output = TensorRank1<D>;
    fn index(&self, index: usize) -> &Self::Output
    {
        &self.0[index]
    }
}

impl<const D: usize> IndexMut<usize> for TensorRank2<D>
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output
    {
        &mut self.0[index]
    }
}