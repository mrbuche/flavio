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
// #[derive(Clone, Copy)]
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
    Self: FromIterator<TensorRank1<D>>
        + Index<usize, Output = TensorRank1<D>>
        + IndexMut<usize, Output = TensorRank1<D>>
        + Mul<TensorRank0, Output = Self>
        + Mul<Output = Self>
        + Sized
        + Sub<Output = Self>
{
    fn determinant(&self) -> TensorRank0
    {
        panic!("Determinant only implemented for D = 2, 3, or 4.");
    }
    fn deviatoric(&self) -> Self;
    fn dyad(vector_a: &TensorRank1<D>, vector_b: &TensorRank1<D>) -> Self;
    fn full_contraction(&self, tensor_rank_2: &Self) -> TensorRank0;
    fn identity() -> Self
    {
        (0..D).map(|i|
            (0..D).map(|j|
                ((i == j) as u8) as TensorRank0
            ).collect()
        ).collect()
    }
    fn inverse(&self) -> Self
    {
        let (tensor_l, tensor_u) = self.lu_decomposition();
        tensor_u.inverse_upper_triangular() * tensor_l.inverse_lower_triangular()
    }
    fn inverse_transpose(&self) -> Self
    {
        self.inverse().transpose()
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
    fn new(array: [[TensorRank0; D]; D]) -> Self
    {
        array.iter().map(|array_i|
            TensorRank1(*array_i)
        ).collect()
    }
    fn norm(&self) -> TensorRank0;
    fn second_invariant(self) -> TensorRank0
    {
        0.5*(self.trace().powi(2) - self.squared().trace())
    }
    fn squared(self) -> Self;
    fn trace(&self) -> TensorRank0
    {
        (0..D).map(|i|
            self[i][i]
        ).sum()
    }
    fn transpose(&self) -> Self
    {
        (0..D).map(|i|
            (0..D).map(|j|
                self[j][i]
            ).collect()
        ).collect()
    }
    fn zero() -> Self;
}

impl<const D: usize> TensorRank2Traits<D> for TensorRank2<D>
{
    fn deviatoric(&self) -> Self
    {
        let factor = -self.trace() / (D as TensorRank0);
        Self::identity() * factor + self
    }
    fn dyad(vector_a: &TensorRank1<D>, vector_b: &TensorRank1<D>) -> Self
    {
        vector_a.iter().map(|vector_a_i|
            vector_b.iter().map(|vector_b_j|
                vector_a_i * vector_b_j
            ).collect()
        ).collect()
    }
    fn full_contraction(&self, tensor_rank_2: &Self) -> TensorRank0
    {
        self.iter().zip(tensor_rank_2.iter()).map(|(self_i, tensor_rank_2_i)|
            self_i * tensor_rank_2_i
        ).sum()
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
    fn norm(&self) -> TensorRank0
    {
        ((self.transpose() * self).trace()/2.0).sqrt()
    }
    fn squared(self) -> Self
    {
        // self * &self
        todo!();
    }
    fn zero() -> Self
    {
        Self(std::array::from_fn(|_| TensorRank1::zero()))
    }
}

// even if get all the traits defaulted
// cannot implement them for all D and then for specific D
// would need to implement only for specific D
// but good to still have all the defaults
// otherwise need to re-implement for each D

impl TensorRank2<2>
{
    fn determinant(&self) -> TensorRank0
    {
        self[0][0] * self[1][1] - self[0][1] * self[1][0]
    }
    fn inverse(&self) -> Self
    {
        Self::new([
            [ self[1][1], -self[0][1]],
            [-self[1][0],  self[0][0]]
        ])/self.determinant()
    }
    fn inverse_transpose(&self) -> Self
    {
        Self::new([
            [ self[1][1], -self[1][0]],
            [-self[0][1],  self[0][0]]
        ])/self.determinant()
    }
}

impl TensorRank2<3>
{
    fn determinant(&self) -> TensorRank0
    {
          self[0][0] * (self[1][1] * self[2][2] - self[1][2] * self[2][1])
        + self[0][1] * (self[1][2] * self[2][0] - self[1][0] * self[2][2])
        + self[0][2] * (self[1][0] * self[2][1] - self[1][1] * self[2][0])
    }
    fn inverse(&self) -> Self
    {
        Self([
            TensorRank1([
                self[1][1] * self[2][2] - self[1][2] * self[2][1],
                self[0][2] * self[2][1] - self[0][1] * self[2][2],
                self[0][1] * self[1][2] - self[0][2] * self[1][1],
            ]),
            TensorRank1([
                self[1][2] * self[2][0] - self[1][0] * self[2][2],
                self[0][0] * self[2][2] - self[0][2] * self[2][0],
                self[0][2] * self[1][0] - self[0][0] * self[1][2],
            ]),
            TensorRank1([
                self[1][0] * self[2][1] - self[1][1] * self[2][0],
                self[0][1] * self[2][0] - self[0][0] * self[2][1],
                self[0][0] * self[1][1] - self[0][1] * self[1][0],
            ])
        ])/self.determinant()
    }
    fn inverse_transpose(&self) -> Self
    {
        TensorRank2([
            TensorRank1([
                self[1][1] * self[2][2] - self[1][2] * self[2][1],
                self[1][2] * self[2][0] - self[1][0] * self[2][2],
                self[1][0] * self[2][1] - self[1][1] * self[2][0],
            ]),
            TensorRank1([
                self[0][2] * self[2][1] - self[0][1] * self[2][2],
                self[0][0] * self[2][2] - self[0][2] * self[2][0],
                self[0][1] * self[2][0] - self[0][0] * self[2][1],
            ]),
            TensorRank1([
                self[0][1] * self[1][2] - self[0][2] * self[1][1],
                self[0][2] * self[1][0] - self[0][0] * self[1][2],
                self[0][0] * self[1][1] - self[0][1] * self[1][0],
            ])
        ])/self.determinant()
    }
}

impl TensorRank2<4>
{
    fn determinant(&self) -> TensorRank0
    {
        todo!();
    }
    fn inverse(&self) -> Self
    {
        todo!();
    }
    fn inverse_transpose(&self) -> Self
    {
        todo!();
    }
}

impl<const D: usize> FromIterator<TensorRank1<D>> for TensorRank2<D>
{
    fn from_iter<I: IntoIterator<Item=TensorRank1<D>>>(into_iterator: I) -> Self
    {
        let mut tensor_rank_2 = Self::zero();
        tensor_rank_2.iter_mut().zip(into_iterator).for_each(|(tensor_rank_2_i, value_i)|
            *tensor_rank_2_i = value_i
        );
        tensor_rank_2
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

impl<const D: usize> Div<TensorRank0> for TensorRank2<D>
{
    type Output = Self;
    fn div(mut self, tensor_rank_0: TensorRank0) -> Self::Output
    {
        self /= &tensor_rank_0;
        self
    }
}

impl<const D: usize> Div<&TensorRank0> for TensorRank2<D>
{
    type Output = Self;
    fn div(mut self, tensor_rank_0: &TensorRank0) -> Self::Output
    {
        self /= tensor_rank_0;
        self
    }
}

impl<const D: usize> DivAssign<TensorRank0> for TensorRank2<D>
{
    fn div_assign(&mut self, tensor_rank_0: TensorRank0)
    {
        self.iter_mut().for_each(|self_i|
            *self_i /= &tensor_rank_0
        );
    }
}

impl<const D: usize> DivAssign<&TensorRank0> for TensorRank2<D>
{
    fn div_assign(&mut self, tensor_rank_0: &TensorRank0)
    {
        self.iter_mut().for_each(|self_i|
            *self_i /= tensor_rank_0
        );
    }
}

impl<const D: usize> Mul<TensorRank0> for TensorRank2<D>
{
    type Output = Self;
    fn mul(mut self, tensor_rank_0: TensorRank0) -> Self::Output
    {
        self *= &tensor_rank_0;
        self
    }
}

impl<const D: usize> Mul<&TensorRank0> for TensorRank2<D>
{
    type Output = Self;
    fn mul(mut self, tensor_rank_0: &TensorRank0) -> Self::Output
    {
        self *= tensor_rank_0;
        self
    }
}

impl<const D: usize> MulAssign<TensorRank0> for TensorRank2<D>
{
    fn mul_assign(&mut self, tensor_rank_0: TensorRank0)
    {
        self.iter_mut().for_each(|self_i|
            *self_i *= &tensor_rank_0
        );
    }
}

impl<const D: usize> MulAssign<&TensorRank0> for TensorRank2<D>
{
    fn mul_assign(&mut self, tensor_rank_0: &TensorRank0)
    {
        self.iter_mut().for_each(|self_i|
            *self_i *= tensor_rank_0
        );
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

impl<const D: usize> Add for TensorRank2<D>
{
    type Output = Self;
    fn add(mut self, tensor_rank_2: Self) -> Self::Output
    {
        self += tensor_rank_2;
        self
    }
}

impl<const D: usize> Add<&Self> for TensorRank2<D>
{
    type Output = Self;
    fn add(mut self, tensor_rank_2: &Self) -> Self::Output
    {
        self += tensor_rank_2;
        self
    }
}

impl<const D: usize> Add<TensorRank2<D>> for &TensorRank2<D>
{
    type Output = TensorRank2<D>;
    fn add(self, mut tensor_rank_2: TensorRank2<D>) -> Self::Output
    {
        tensor_rank_2 += self;
        tensor_rank_2
    }
}

impl<const D: usize> AddAssign for TensorRank2<D>
{
    fn add_assign(&mut self, tensor_rank_2: Self)
    {
        self.iter_mut().zip(tensor_rank_2.iter()).for_each(|(self_i, tensor_rank_2_i)|
            *self_i += tensor_rank_2_i
        );
    }
}

impl<const D: usize> AddAssign<&Self> for TensorRank2<D>
{
    fn add_assign(&mut self, tensor_rank_2: &Self)
    {
        self.iter_mut().zip(tensor_rank_2.iter()).for_each(|(self_i, tensor_rank_2_i)|
            *self_i += tensor_rank_2_i
        );
    }
}

impl<const D: usize> Mul for TensorRank2<D>
{
    type Output = Self;
    fn mul(self, tensor_rank_2: Self) -> Self::Output
    {
        tensor_rank_2.transpose().iter().map(|tensor_rank_2_transpose_j|
            self.iter().map(|self_i|
                self_i * tensor_rank_2_transpose_j
            ).collect()
        ).collect()
    }
}

impl<const D: usize> Mul<&Self> for TensorRank2<D>
{
    type Output = Self;
    fn mul(self, tensor_rank_2: &Self) -> Self::Output
    {
        tensor_rank_2.transpose().iter().map(|tensor_rank_2_transpose_j|
            self.iter().map(|self_i|
                self_i * tensor_rank_2_transpose_j
            ).collect()
        ).collect()
    }
}

impl<const D: usize> Mul<TensorRank2<D>> for &TensorRank2<D>
{
    type Output = TensorRank2<D>;
    fn mul(self, tensor_rank_2: TensorRank2<D>) -> Self::Output
    {
        tensor_rank_2.transpose().iter().map(|tensor_rank_2_transpose_j|
            self.iter().map(|self_i|
                self_i * tensor_rank_2_transpose_j
            ).collect()
        ).collect()
    }
}

impl<const D: usize> Sub for TensorRank2<D>
{
    type Output = Self;
    fn sub(mut self, tensor_rank_2: Self) -> Self::Output
    {
        self -= tensor_rank_2;
        self
    }
}

impl<const D: usize> Sub<&Self> for TensorRank2<D>
{
    type Output = Self;
    fn sub(mut self, tensor_rank_2: &Self) -> Self::Output
    {
        self -= tensor_rank_2;
        self
    }
}

// impl<const D: usize> Sub<TensorRank2<D>> for &TensorRank2<D>
// {
//     type Output = TensorRank2<D>;
//     fn sub(self, mut tensor_rank_2: TensorRank2<D>) -> Self::Output
//     {
//         tensor_rank_2.iter_mut().zip(self.iter()).for_each(|(tensor_rank_2_i, self_i)|
//             *tensor_rank_2_i = self_i - *tensor_rank_2_i
//         );
//         tensor_rank_2
//     }
// }

impl<const D: usize> SubAssign for TensorRank2<D>
{
    fn sub_assign(&mut self, tensor_rank_2: Self)
    {
        self.iter_mut().zip(tensor_rank_2.iter()).for_each(|(self_i, tensor_rank_2_i)|
            *self_i -= tensor_rank_2_i
        );
    }
}

impl<const D: usize> SubAssign<&Self> for TensorRank2<D>
{
    fn sub_assign(&mut self, tensor_rank_2: &Self)
    {
        self.iter_mut().zip(tensor_rank_2.iter()).for_each(|(self_i, tensor_rank_2_i)|
            *self_i -= tensor_rank_2_i
        );
    }
}