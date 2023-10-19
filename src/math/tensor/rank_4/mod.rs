#[cfg(test)]
mod test;

use std::ops::
{
    Div,
    DivAssign,
    Index,
    IndexMut,
    Mul,
    MulAssign
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

pub trait TensorRank4Traits<const D: usize, T2A, T2B>
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
    fn dyad_ij_kl(tensor_rank_2_a: T2A, tensor_rank_2_b: T2B) -> Self;
    fn dyad_ik_jl(tensor_rank_2_a: T2A, tensor_rank_2_b: T2B) -> Self;
    fn dyad_il_jk(tensor_rank_2_a: T2A, tensor_rank_2_b: T2B) -> Self;
    fn dyad_il_kj(tensor_rank_2_a: T2A, tensor_rank_2_b: T2B) -> Self;
    fn zero() -> Self;
}

impl<const D: usize> TensorRank4Traits<D, &TensorRank2<D>, &TensorRank2<D>> for TensorRank4<D>
{
    fn dyad_ij_kl(tensor_rank_2_a: &TensorRank2<D>, tensor_rank_2_b: &TensorRank2<D>) -> Self
    {
        tensor_rank_2_a.iter().map(|tensor_rank_2_a_i|
            tensor_rank_2_a_i.iter().map(|tensor_rank_2_a_ij|
                tensor_rank_2_b * tensor_rank_2_a_ij
            ).collect()
        ).collect()
    }
    fn dyad_ik_jl(tensor_rank_2_a: &TensorRank2<D>, tensor_rank_2_b: &TensorRank2<D>) -> Self
    {
        tensor_rank_2_a.iter().map(|tensor_rank_2_a_i|
            tensor_rank_2_b.iter().map(|tensor_rank_2_b_j|
                tensor_rank_2_a_i.iter().map(|tensor_rank_2_a_ik|
                    tensor_rank_2_b_j * tensor_rank_2_a_ik
                ).collect()
            ).collect()
        ).collect()
    }
    fn dyad_il_jk(tensor_rank_2_a: &TensorRank2<D>, tensor_rank_2_b: &TensorRank2<D>) -> Self
    {
        tensor_rank_2_a.iter().map(|tensor_rank_2_a_i|
            tensor_rank_2_b.iter().map(|tensor_rank_2_b_j|
                tensor_rank_2_b_j.iter().map(|tensor_rank_2_b_jk|
                    tensor_rank_2_a_i * tensor_rank_2_b_jk
                ).collect()
            ).collect()
        ).collect()
    }
    fn dyad_il_kj(tensor_rank_2_a: &TensorRank2<D>, tensor_rank_2_b: &TensorRank2<D>) -> Self
    {
        todo!();
    }
    fn zero() -> Self
    {
        Self(std::array::from_fn(|_| TensorRank3::zero()))
    }
}

// try to make use of TensorRank2::sum() for some of the contractions !!!!!!!!!!!!!!!!!!!!!!!!!

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

impl<const D: usize> Div<TensorRank0> for TensorRank4<D>
{
    type Output = Self;
    fn div(mut self, tensor_rank_0: TensorRank0) -> Self::Output
    {
        self /= &tensor_rank_0;
        self
    }
}

impl<const D: usize> Div<&TensorRank0> for TensorRank4<D>
{
    type Output = Self;
    fn div(mut self, tensor_rank_0: &TensorRank0) -> Self::Output
    {
        self /= tensor_rank_0;
        self
    }
}

impl<const D: usize> DivAssign<TensorRank0> for TensorRank4<D>
{
    fn div_assign(&mut self, tensor_rank_0: TensorRank0)
    {
        self.iter_mut().for_each(|self_i|
            *self_i /= &tensor_rank_0
        );
    }
}

impl<const D: usize> DivAssign<&TensorRank0> for TensorRank4<D>
{
    fn div_assign(&mut self, tensor_rank_0: &TensorRank0)
    {
        self.iter_mut().for_each(|self_i|
            *self_i /= tensor_rank_0
        );
    }
}

impl<const D: usize> Mul<TensorRank0> for TensorRank4<D>
{
    type Output = Self;
    fn mul(mut self, tensor_rank_0: TensorRank0) -> Self::Output
    {
        self *= &tensor_rank_0;
        self
    }
}

impl<const D: usize> Mul<&TensorRank0> for TensorRank4<D>
{
    type Output = Self;
    fn mul(mut self, tensor_rank_0: &TensorRank0) -> Self::Output
    {
        self *= tensor_rank_0;
        self
    }
}

impl<const D: usize> MulAssign<TensorRank0> for TensorRank4<D>
{
    fn mul_assign(&mut self, tensor_rank_0: TensorRank0)
    {
        self.iter_mut().for_each(|self_i|
            *self_i *= &tensor_rank_0
        );
    }
}

impl<const D: usize> MulAssign<&TensorRank0> for TensorRank4<D>
{
    fn mul_assign(&mut self, tensor_rank_0: &TensorRank0)
    {
        self.iter_mut().for_each(|self_i|
            *self_i *= tensor_rank_0
        );
    }
}
