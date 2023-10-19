#[cfg(test)]
mod test;

use std::ops::
{
    Index,
    IndexMut
};

use super::
{
    rank_0::TensorRank0,
    rank_2::
    {
        TensorRank2,
        TensorRank2Traits
    },
    rank_3::TensorRank3
};

pub struct TensorRank4<const D: usize>
(
    pub [TensorRank3<D>; D]
);

impl<const D: usize> TensorRank4<D>
{
    pub fn iter(&self) -> impl Iterator<Item=&TensorRank3<D>>
    {
        self.0.iter()
    }
    pub fn iter_mut(&mut self) -> impl Iterator<Item=&mut TensorRank3<D>>
    {
        self.0.iter_mut()
    }
    pub fn zero() -> Self
    {
        Self(std::array::from_fn(|_| TensorRank3::zero()))
    }
}

pub trait TensorRank4Traits<const D: usize>
where
    Self: FromIterator<TensorRank3<D>>
        + Index<usize, Output = TensorRank3<D>>
        + IndexMut<usize, Output = TensorRank3<D>>
        + Sized
{
    fn determinant(&self) -> TensorRank0
    {
        panic!()
    }
    fn dyad_ij_kl(tensor_rank_2_a: &TensorRank2<D>, tensor_rank_2_b: &TensorRank2<D>) -> Self
    {
        tensor_rank_2_a.iter().map(|tensor_rank_2_a_i|
            tensor_rank_2_a_i.iter().map(|tensor_rank_2_a_ij|
                tensor_rank_2_b.iter().map(|tensor_rank_2_b_k|
                    tensor_rank_2_b_k.iter().map(|tensor_rank_2_b_kl|
                        tensor_rank_2_a_ij * tensor_rank_2_b_kl
                    ).collect()
                ).collect()
            ).collect()
        ).collect()
    }
    fn dyad_ik_jl(tensor_rank_2_a: &TensorRank2<D>, tensor_rank_2_b: &TensorRank2<D>) -> Self
    {
        tensor_rank_2_a.iter().map(|tensor_rank_2_a_i|
            tensor_rank_2_b.iter().map(|tensor_rank_2_b_j|
                tensor_rank_2_a_i.iter().map(|tensor_rank_2_a_ik|
                    tensor_rank_2_b_j.iter().map(|tensor_rank_2_b_jl|
                        tensor_rank_2_a_ik * tensor_rank_2_b_jl
                    ).collect()
                ).collect()
            ).collect()
        ).collect()
    }
    fn dyad_il_jk(tensor_rank_2_a: &TensorRank2<D>, tensor_rank_2_b: &TensorRank2<D>) -> Self
    {
        tensor_rank_2_a.iter().map(|tensor_rank_2_a_i|
            tensor_rank_2_b.iter().map(|tensor_rank_2_b_j|
                tensor_rank_2_b_j.iter().map(|tensor_rank_2_b_jk|
                    tensor_rank_2_a_i.iter().map(|tensor_rank_2_a_il|
                        tensor_rank_2_a_il * tensor_rank_2_b_jk
                    ).collect()
                ).collect()
            ).collect()
        ).collect()
    }
    fn dyad_il_kj(tensor_rank_2_a: &TensorRank2<D>, tensor_rank_2_b: &TensorRank2<D>) -> Self;
    fn zero() -> Self;
}

impl TensorRank4Traits<9> for TensorRank4<9>
{
    fn dyad_il_kj(tensor_rank_2_a: &TensorRank2<9>, tensor_rank_2_b: &TensorRank2<9>) -> Self
    {
        Self::dyad_il_jk(tensor_rank_2_a, &(tensor_rank_2_b.transpose()))
    }
    fn zero() -> Self
    {
        Self(std::array::from_fn(|_| TensorRank3::zero()))
    }
}

impl<const D: usize> FromIterator<TensorRank3<D>> for TensorRank4<D>
{
    fn from_iter<I: IntoIterator<Item=TensorRank3<D>>>(into_iterator: I) -> Self
    {
        let mut tensor_rank_4 = Self::zero();
        tensor_rank_4.iter_mut().zip(into_iterator).for_each(|(tensor_rank_4_i, value_i)|
            *tensor_rank_4_i = value_i
        );
        tensor_rank_4
    }
}

impl<const D: usize> Index<usize> for TensorRank4<D>
{
    type Output = TensorRank3<D>;
    fn index(&self, index: usize) -> &Self::Output
    {
        &self.0[index]
    }
}

impl<const D: usize> IndexMut<usize> for TensorRank4<D>
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output
    {
        &mut self.0[index]
    }
}
