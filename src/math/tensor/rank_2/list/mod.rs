#[cfg(test)]
mod test;

use super::
{
    TensorRank0,
    TensorRank2,
    TensorRank2Traits
};

pub struct TensorRank2List<const D: usize, const L: usize>
(
    pub [TensorRank2<D>; L]
);

impl<const D: usize, const L: usize> TensorRank2List<D, L>
{
    pub fn iter(&self) -> impl Iterator<Item = &TensorRank2<D>>
    {
        self.0.iter()
    }
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut TensorRank2<D>>
    {
        self.0.iter_mut()
    }
}

pub trait TensorRank2ListTraits<const D: usize, const L: usize>
{}