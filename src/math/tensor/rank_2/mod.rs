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

// not as big of a deal to have all the default trait implementations
// if just going to redo them in a struct-tuple deal as Tensor<D, I, J>
// only helps with things like the inverse, which want to be overridden for some specific <D> implementations here
pub trait TensorRank2Traits<const D: usize>
where
    Self: Index<usize, Output = TensorRank1<D>>
        + IndexMut<usize, Output = TensorRank1<D>>
        + Mul<Output = Self>
        + Sized
{
    fn construct_from_dyad_i_j(vector_a: TensorRank1<D>, vector_b: TensorRank1<D>) -> Self;
    fn full_contraction_with(&self, tensor_rank_2: &Self) -> TensorRank0;
    fn frobenius_norm(&self) -> TensorRank0;
    fn inverse(&self) -> Self
    {
        let (tensor_l, tensor_u) = self.lu_decomposition();
        tensor_u.inverse_upper_triangular() * tensor_l.inverse_lower_triangular()
    }
    fn inverse_lower_triangular(mut self) -> Self
    {
        let mut sum;
        for i in 0..D
        {
            self[i][i] = 1.0/self[i][i];
            for j in 0..i
            {
                sum = 0.0;
                for k in j..i
                {
                    sum += self[i][k] * self[k][j];
                }
                self[i][j] = -sum * self[i][i];
            }
        }
        self
    }
    fn inverse_upper_triangular(mut self) -> Self
    {
        let mut sum;
        for i in 0..D
        {
            self[i][i] = 1.0/self[i][i];
            for j in 0..i
            {
                sum = 0.0;
                for k in j..i
                {
                    sum += self[j][k] * self[k][i];
                }
                self[j][i] = -sum * self[i][i];
            }
        }
        self
    }
    fn lu_decomposition(&self) -> (Self, Self);
    fn second_invariant(self) -> TensorRank0
    {
        0.5*(self.trace().powi(2) - self.squared().trace())
    }
    fn squared(self) -> Self;
    fn trace(&self) -> TensorRank0
    {
        (0..D).map(|i| self[i][i]).sum()
    }
    fn transpose(&self) -> Self;
    fn zero() -> Self;
}

impl<const D: usize> TensorRank2Traits<D> for TensorRank2<D>
{
    fn construct_from_dyad_i_j(vector_a: TensorRank1<D>, vector_b: TensorRank1<D>) -> Self
    {
        let mut output = TensorRank2::<D>::zero();
        output.iter_mut().zip(vector_a.iter()).for_each(|(output_i, vector_a_i)|
            output_i.iter_mut().zip(vector_b.iter()).for_each(|(output_ij, vector_b_j)|
                *output_ij = vector_a_i * vector_b_j
            )
        );
        output
    }
    fn full_contraction_with(&self, tensor_rank_2: &Self) -> TensorRank0
    {
        self.iter().zip(tensor_rank_2.iter()).map(|(self_i, tensor_rank_2_i)|
            self_i.iter().zip(tensor_rank_2_i.iter()).map(|(self_ij, tensor_rank_2_ij)|
                self_ij * tensor_rank_2_ij
            ).sum::<TensorRank0>()
        ).sum()
    }
    fn frobenius_norm(&self) -> TensorRank0
    {
        (self.transpose()*self).trace().sqrt()
    }
    fn lu_decomposition(&self) -> (Self, Self)
    {
        let mut tensor_l = Self::zero();
        let mut tensor_u = Self::zero();
        self.iter().enumerate().for_each(|(i, self_i)|
        {
            self.iter().enumerate().zip(tensor_l.iter_mut()).for_each(|((j, self_j), tensor_l_j)|
            {
                if j >= i
                {
                    tensor_l_j[i] = self_j[i];
                    for k in 0..i
                    {
                        tensor_l_j[i] -= tensor_l_j[k] * tensor_u[k][i];
                    }
                }
            });
            self_i.iter().enumerate().for_each(|(j, self_ij)|
            {
                if j == i
                {
                    tensor_u[i][j] = 1.0;
                }
                else if j > i
                {
                    tensor_u[i][j] = self_ij / tensor_l[i][i];
                    for k in 0..i
                    {
                        tensor_u[i][j] -= (tensor_l[i][k] * tensor_u[k][j]) / tensor_l[i][i];
                    }
                }
            });
        });
        // for i in 0..D
        // {
        //     for j in 0..D
        //     {
        //         if j >= i
        //         {
        //             tensor_l[j][i] = self[j][i];
        //             for k in 0..i
        //             {
        //                 tensor_l[j][i] -=  tensor_l[j][k] * tensor_u[k][i];
        //             }
        //         }
        //     }
        //     for j in 0..D
        //     {
        //         if j == i
        //         {
        //             tensor_u[i][j] = 1.0;
        //         }
        //         else if j > i
        //         {
        //             tensor_u[i][j] = self[i][j] / tensor_l[i][i];
        //             for k in 0..i
        //             {
        //                 tensor_u[i][j] -= (tensor_l[i][k] * tensor_u[k][j]) / tensor_l[i][i];
        //             }
        //         }
        //     }
        // }
        (tensor_l, tensor_u)
    }
    fn squared(self) -> Self
    {
        self * &self
    }
    fn transpose(&self) -> Self
    {
        let mut transpose = Self::zero();
        transpose.iter_mut().enumerate().for_each(|(i, transpose_i)|
            transpose_i.iter_mut().enumerate().for_each(|(j, transpose_ij)|
                *transpose_ij = self[j][i]
            )
        );
        transpose
    }
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

impl<const D: usize> Mul<TensorRank0> for TensorRank2<D>
{
    type Output = Self;
    fn mul(mut self, tensor_rank_0: TensorRank0) -> Self::Output
    {
        self.iter_mut().for_each(|self_i|
            self_i.iter_mut().for_each(|self_ij|
                *self_ij *= &tensor_rank_0
            )
        );
        self
    }
}

impl<const D: usize> Mul<&TensorRank0> for TensorRank2<D>
{
    type Output = Self;
    fn mul(mut self, tensor_rank_0: &TensorRank0) -> Self::Output
    {
        self.iter_mut().for_each(|self_i|
            self_i.iter_mut().for_each(|self_ij|
                *self_ij *= tensor_rank_0
            )
        );
        self
    }
}

impl<const D: usize> Div<TensorRank0> for TensorRank2<D>
{
    type Output = Self;
    fn div(mut self, tensor_rank_0: TensorRank0) -> Self::Output
    {
        self.iter_mut().for_each(|self_i|
            self_i.iter_mut().for_each(|self_ij|
                *self_ij /= &tensor_rank_0
            )
        );
        self
    }
}

impl<const D: usize> Div<&TensorRank0> for TensorRank2<D>
{
    type Output = Self;
    fn div(mut self, tensor_rank_0: &TensorRank0) -> Self::Output
    {
        self.iter_mut().for_each(|self_i|
            self_i.iter_mut().for_each(|self_ij|
                *self_ij /= tensor_rank_0
            )
        );
        self
    }
}

impl<const D: usize> Mul<TensorRank1<D>> for TensorRank2<D>
{
    type Output = TensorRank1<D>;
    fn mul(self, tensor_rank_1: TensorRank1<D>) -> Self::Output
    {
        self.iter().map(|self_i|
            self_i * &tensor_rank_1
        ).collect()
    }
}

impl<const D: usize> Mul<&TensorRank1<D>> for TensorRank2<D>
{
    type Output = TensorRank1<D>;
    fn mul(self, tensor_rank_1: &TensorRank1<D>) -> Self::Output
    {
        self.iter().map(|self_i|
            self_i * tensor_rank_1
        ).collect()
    }
}

// look for places to use map(), sum(), and collect()

impl<const D: usize> Mul for TensorRank2<D>
{
    type Output = Self;
    fn mul(self, tensor_rank_2: Self) -> Self::Output
    {
        let mut output = TensorRank2::zero();
        let tensor_rank_2_transpose = tensor_rank_2.transpose();
        output.iter_mut().zip(self.iter()).for_each(|(output_i, self_i)|
            output_i.iter_mut().zip(tensor_rank_2_transpose.iter()).for_each(|(output_ij, tensor_rank_2_transpose_j)|
                *output_ij = *self_i * tensor_rank_2_transpose_j
            )
        );
        output
    }
}
impl<const D: usize> Mul<&Self> for TensorRank2<D>
{
    type Output = Self;
    fn mul(self, tensor_rank_2: &Self) -> Self::Output
    {
        let mut output = TensorRank2::zero();
        let tensor_rank_2_transpose = tensor_rank_2.transpose();
        output.iter_mut().zip(self.iter()).for_each(|(output_i, self_i)|
            output_i.iter_mut().zip(tensor_rank_2_transpose.iter()).for_each(|(output_ij, tensor_rank_2_transpose_j)|
                *output_ij = *self_i * tensor_rank_2_transpose_j
            )
        );
        output
    }
}