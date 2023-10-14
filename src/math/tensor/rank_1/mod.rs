#[cfg(test)]
mod test;

use super::
{
    rank_0::TensorRank0,
    // rank_2::TensorRank2
};

// eliminate in order to identify explicit copying at some point
#[derive(Clone, Copy)]
pub struct TensorRank1<const D: usize>
(
    pub [TensorRank0; D]
);


// move into TensorRank1Traits if ever becomes possible
impl<const D: usize> TensorRank1<D>
{
    fn iter(&self) -> impl Iterator<Item=&TensorRank0>
    {
        self.0.iter()
    }
    fn iter_mut(&mut self) -> impl Iterator<Item=&mut TensorRank0>
    {
        self.0.iter_mut()
    }
}

pub trait TensorRank1Traits
{
    fn zero() -> Self;
}

impl<const D: usize> TensorRank1Traits for TensorRank1<D>
{
    fn zero() -> Self
    {
        TensorRank1([0.0; D])
    }
}
