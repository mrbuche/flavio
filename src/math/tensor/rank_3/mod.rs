#[cfg(test)]
mod test;

use super::
{
    rank_2::
    {
        TensorRank2,
        TensorRank2Traits
    }
};

pub struct TensorRank3<const D: usize>
(
    pub [TensorRank2<D>; D]
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
