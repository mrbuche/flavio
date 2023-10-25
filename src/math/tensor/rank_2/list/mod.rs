#[cfg(test)]
mod test;

use std::ops::
{
    Index,
    IndexMut
};

use super::
{
    TensorRank0,
    TensorRank1,
    TensorRank2
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
where
    Self: FromIterator<TensorRank2<D>>
        + Index<usize, Output = TensorRank2<D>>
        + IndexMut<usize, Output = TensorRank2<D>>
        + Sized
{
    fn new(array: [[[TensorRank0; D]; D]; L]) -> Self
    {
        array.iter().map(|array_i|
            array_i.iter().map(|array_ij|
                TensorRank1(*array_ij)
            ).collect()
        ).collect()
    }
    fn zero() -> Self;
}

impl<const D: usize, const L: usize> TensorRank2ListTraits<D, L> for TensorRank2List<D, L>
{
    fn zero() -> Self
    {
        Self(std::array::from_fn(|_| TensorRank2::zero()))
    }
}

impl<const D: usize, const L: usize> FromIterator<TensorRank2<D>> for TensorRank2List<D, L>
{
    fn from_iter<I: IntoIterator<Item=TensorRank2<D>>>(into_iterator: I) -> Self
    {
        let mut tensor_rank_2_list = Self::zero();
        tensor_rank_2_list.iter_mut().zip(into_iterator).for_each(|(tensor_rank_2_list_entry, entry)|
            *tensor_rank_2_list_entry = entry
        );
        tensor_rank_2_list
    }
}

impl<const D: usize, const L: usize> Index<usize> for TensorRank2List<D, L>
{
    type Output = TensorRank2<D>;
    fn index(&self, index: usize) -> &Self::Output
    {
        &self.0[index]
    }
}

impl<const D: usize, const L: usize> IndexMut<usize> for TensorRank2List<D, L>
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output
    {
        &mut self.0[index]
    }
}
