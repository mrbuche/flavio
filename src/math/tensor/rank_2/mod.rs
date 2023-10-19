#[cfg(test)]
mod test;

use std::cmp::Ordering;
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

pub struct TensorRank2<const D: usize>
(
    pub [TensorRank1<D>; D]
);

impl<const D: usize> TensorRank2<D>
{
    pub fn iter(&self) -> impl Iterator<Item = &TensorRank1<D>>
    {
        self.0.iter()
    }
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut TensorRank1<D>>
    {
        self.0.iter_mut()
    }
    pub fn zero() -> Self
    {
        Self(std::array::from_fn(|_| TensorRank1::zero()))
    }
}

pub trait TensorRank2Traits<const D: usize>
where
    Self: FromIterator<TensorRank1<D>>
        + Index<usize, Output = TensorRank1<D>>
        + IndexMut<usize, Output = TensorRank1<D>>
        + Sized
{
    fn determinant(&self) -> TensorRank0;
    fn deviatoric(&self) -> Self;
    fn deviatoric_and_trace(&self) -> (Self, TensorRank0);
    fn dyad(vector_a: &TensorRank1<D>, vector_b: &TensorRank1<D>) -> Self
    {
        vector_a.iter().map(|vector_a_i|
            vector_b.iter().map(|vector_b_j|
                vector_a_i * vector_b_j
            ).collect()
        ).collect()
    }
    fn full_contraction(&self, tensor_rank_2: &Self) -> TensorRank0;
    fn identity() -> Self
    {
        (0..D).map(|i|
            (0..D).map(|j|
                ((i == j) as u8) as TensorRank0
            ).collect()
        ).collect()
    }
    fn inverse(&self) -> Self;
    fn inverse_and_determinant(&self) -> (Self, TensorRank0);
    fn inverse_transpose(&self) -> Self;
    fn inverse_transpose_and_determinant(&self) -> (Self, TensorRank0);
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
    fn lu_decomposition(&self) -> (Self, Self)
    {
        let mut tensor_l = Self::zero();
        let mut tensor_u = Self::zero();
        for i in 0..D
        {
            for j in 0..D
            {
                if j >= i
                {
                    tensor_l[j][i] = self[j][i];
                    for k in 0..i
                    {
                        tensor_l[j][i] -=  tensor_l[j][k] * tensor_u[k][i];
                    }
                }
            }
            for j in 0..D
            {
                match j.cmp(&i) {
                    Ordering::Equal =>
                    {
                        tensor_u[i][j] = 1.0;
                    }
                    Ordering::Greater =>
                    {
                        tensor_u[i][j] = self[i][j] / tensor_l[i][i];
                        for k in 0..i
                        {
                            tensor_u[i][j] -= (tensor_l[i][k] * tensor_u[k][j]) / tensor_l[i][i];
                        }
                    }
                    Ordering::Less => ()
                }
            }
        }
        (tensor_l, tensor_u)
    }
    fn new(array: [[TensorRank0; D]; D]) -> Self
    {
        array.iter().map(|array_i|
            TensorRank1(*array_i)
        ).collect()
    }
    fn norm(&self) -> TensorRank0;
    fn second_invariant(&self) -> TensorRank0
    {
        0.5*(self.trace().powi(2) - self.squared_trace())
    }
    fn squared_trace(&self) -> TensorRank0;
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

impl TensorRank2Traits<2> for TensorRank2<2>
{
    fn determinant(&self) -> TensorRank0
    {
        self[0][0] * self[1][1] - self[0][1] * self[1][0]
    }
    fn deviatoric(&self) -> Self
    {
        Self::identity() * (self.trace() / -2.0) + self
    }
    fn deviatoric_and_trace(&self) -> (Self, TensorRank0)
    {
        let trace = self.trace();
        (
            Self::identity() * (trace / -2.0) + self,
            trace
        )
    }
    fn full_contraction(&self, tensor_rank_2: &Self) -> TensorRank0
    {
        self.iter().zip(tensor_rank_2.iter()).map(|(self_i, tensor_rank_2_i)|
            self_i.iter().zip(tensor_rank_2_i.iter()).map(|(self_ij, tensor_rank_2_ij)|
                self_ij * tensor_rank_2_ij
            ).sum::<TensorRank0>()
        ).sum()
    }
    fn inverse(&self) -> Self
    {
        Self::new([
            [ self[1][1], -self[0][1]],
            [-self[1][0],  self[0][0]]
        ]) / self.determinant()
    }
    fn inverse_and_determinant(&self) -> (Self, TensorRank0)
    {
        let determinant = self.determinant();
        (
            Self::new([
                [ self[1][1], -self[0][1]],
                [-self[1][0],  self[0][0]]
            ]) / determinant,
            determinant
        )
    }
    fn inverse_transpose(&self) -> Self
    {
        Self::new([
            [ self[1][1], -self[1][0]],
            [-self[0][1],  self[0][0]]
        ]) / self.determinant()
    }
    fn inverse_transpose_and_determinant(&self) -> (Self, TensorRank0)
    {
        let determinant = self.determinant();
        (
            Self::new([
                [ self[1][1], -self[1][0]],
                [-self[0][1],  self[0][0]]
            ]) / determinant,
            determinant
        )
    }
    fn norm(&self) -> TensorRank0
    {
        (self.iter().map(|self_i|
            self_i.iter().map(|self_ij|
                self_ij.powi(2)
            ).sum::<TensorRank0>()
        ).sum::<TensorRank0>()/2.0).sqrt()
    }
    fn squared_trace(&self) -> TensorRank0
    {
        self.iter().enumerate().map(|(i, self_i)|
            self_i.iter().zip(self.iter()).map(|(self_ij, self_j)|
                self_ij * self_j[i]
            ).sum::<TensorRank0>()
        ).sum()
    }
    fn zero() -> Self
    {
        Self(std::array::from_fn(|_| TensorRank1::zero()))
    }
}

impl TensorRank2Traits<3> for TensorRank2<3>
{
    fn determinant(&self) -> TensorRank0
    {
        let c_00 = self[1][1] * self[2][2] - self[1][2] * self[2][1];
        let c_10 = self[1][2] * self[2][0] - self[1][0] * self[2][2];
        let c_20 = self[1][0] * self[2][1] - self[1][1] * self[2][0];
        self[0][0] * c_00 + self[0][1] * c_10 + self[0][2] * c_20
    }
    fn deviatoric(&self) -> Self
    {
        Self::identity() * (self.trace() / -3.0) + self
    }
    fn deviatoric_and_trace(&self) -> (Self, TensorRank0)
    {
        let trace = self.trace();
        (
            Self::identity() * (trace / -3.0) + self,
            trace
        )
    }
    fn full_contraction(&self, tensor_rank_2: &Self) -> TensorRank0
    {
        self.iter().zip(tensor_rank_2.iter()).map(|(self_i, tensor_rank_2_i)|
            self_i.iter().zip(tensor_rank_2_i.iter()).map(|(self_ij, tensor_rank_2_ij)|
                self_ij * tensor_rank_2_ij
            ).sum::<TensorRank0>()
        ).sum()
    }
    fn inverse(&self) -> Self
    {
        let c_00 = self[1][1] * self[2][2] - self[1][2] * self[2][1];
        let c_10 = self[1][2] * self[2][0] - self[1][0] * self[2][2];
        let c_20 = self[1][0] * self[2][1] - self[1][1] * self[2][0];
        Self([
            TensorRank1([
                c_00,
                self[0][2] * self[2][1] - self[0][1] * self[2][2],
                self[0][1] * self[1][2] - self[0][2] * self[1][1],
            ]),
            TensorRank1([
                c_10,
                self[0][0] * self[2][2] - self[0][2] * self[2][0],
                self[0][2] * self[1][0] - self[0][0] * self[1][2],
            ]),
            TensorRank1([
                c_20,
                self[0][1] * self[2][0] - self[0][0] * self[2][1],
                self[0][0] * self[1][1] - self[0][1] * self[1][0],
            ])
        ]) / (self[0][0] * c_00 + self[0][1] * c_10 + self[0][2] * c_20)
    }
    fn inverse_and_determinant(&self) -> (Self, TensorRank0)
    {
        let c_00 = self[1][1] * self[2][2] - self[1][2] * self[2][1];
        let c_10 = self[1][2] * self[2][0] - self[1][0] * self[2][2];
        let c_20 = self[1][0] * self[2][1] - self[1][1] * self[2][0];
        let determinant = (self[0][0] * c_00 + self[0][1] * c_10 + self[0][2] * c_20);
        (
            Self([
                TensorRank1([
                    c_00,
                    self[0][2] * self[2][1] - self[0][1] * self[2][2],
                    self[0][1] * self[1][2] - self[0][2] * self[1][1],
                ]),
                TensorRank1([
                    c_10,
                    self[0][0] * self[2][2] - self[0][2] * self[2][0],
                    self[0][2] * self[1][0] - self[0][0] * self[1][2],
                ]),
                TensorRank1([
                    c_20,
                    self[0][1] * self[2][0] - self[0][0] * self[2][1],
                    self[0][0] * self[1][1] - self[0][1] * self[1][0],
                ])
            ]) / determinant,
            determinant
        )
    }
    fn inverse_transpose(&self) -> Self
    {
        let c_00 = self[1][1] * self[2][2] - self[1][2] * self[2][1];
        let c_10 = self[1][2] * self[2][0] - self[1][0] * self[2][2];
        let c_20 = self[1][0] * self[2][1] - self[1][1] * self[2][0];
        Self([
            TensorRank1([
                c_00,
                c_10,
                c_20,
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
        ]) / (self[0][0] * c_00 + self[0][1] * c_10 + self[0][2] * c_20)
    }
    fn inverse_transpose_and_determinant(&self) -> (Self, TensorRank0)
    {
        let c_00 = self[1][1] * self[2][2] - self[1][2] * self[2][1];
        let c_10 = self[1][2] * self[2][0] - self[1][0] * self[2][2];
        let c_20 = self[1][0] * self[2][1] - self[1][1] * self[2][0];
        let determinant = (self[0][0] * c_00 + self[0][1] * c_10 + self[0][2] * c_20);
        (
            Self([
                TensorRank1([
                    c_00,
                    c_10,
                    c_20,
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
            ]) / determinant,
            determinant
        )
    }
    fn norm(&self) -> TensorRank0
    {
        (self.iter().map(|self_i|
            self_i.iter().map(|self_ij|
                self_ij.powi(2)
            ).sum::<TensorRank0>()
        ).sum::<TensorRank0>()/2.0).sqrt()
    }
    fn squared_trace(&self) -> TensorRank0
    {
        self.iter().enumerate().map(|(i, self_i)|
            self_i.iter().zip(self.iter()).map(|(self_ij, self_j)|
                self_ij * self_j[i]
            ).sum::<TensorRank0>()
        ).sum()
    }
    fn zero() -> Self
    {
        Self(std::array::from_fn(|_| TensorRank1::zero()))
    }
}

impl TensorRank2Traits<4> for TensorRank2<4>
{
    fn determinant(&self) -> TensorRank0
    {
        let s0 = self[0][0] * self[1][1] - self[0][1] * self[1][0];
        let s1 = self[0][0] * self[1][2] - self[0][2] * self[1][0];
        let s2 = self[0][0] * self[1][3] - self[0][3] * self[1][0];
        let s3 = self[0][1] * self[1][2] - self[0][2] * self[1][1];
        let s4 = self[0][1] * self[1][3] - self[0][3] * self[1][1];
        let s5 = self[0][2] * self[1][3] - self[0][3] * self[1][2];
        let c5 = self[2][2] * self[3][3] - self[2][3] * self[3][2];
        let c4 = self[2][1] * self[3][3] - self[2][3] * self[3][1];
        let c3 = self[2][1] * self[3][2] - self[2][2] * self[3][1];
        let c2 = self[2][0] * self[3][3] - self[2][3] * self[3][0];
        let c1 = self[2][0] * self[3][2] - self[2][2] * self[3][0];
        let c0 = self[2][0] * self[3][1] - self[2][1] * self[3][0];
        s0 * c5 - s1 * c4 + s2 * c3 + s3 * c2 - s4 * c1 + s5 * c0
    }
    fn deviatoric(&self) -> Self
    {
        Self::identity() * (self.trace() / -4.0) + self
    }
    fn deviatoric_and_trace(&self) -> (Self, TensorRank0)
    {
        let trace = self.trace();
        (
            Self::identity() * (trace / -4.0) + self,
            trace
        )
    }
    fn full_contraction(&self, tensor_rank_2: &Self) -> TensorRank0
    {
        self.iter().zip(tensor_rank_2.iter()).map(|(self_i, tensor_rank_2_i)|
            self_i.iter().zip(tensor_rank_2_i.iter()).map(|(self_ij, tensor_rank_2_ij)|
                self_ij * tensor_rank_2_ij
            ).sum::<TensorRank0>()
        ).sum()
    }
    fn inverse(&self) -> Self
    {
        let s0 = self[0][0] * self[1][1] - self[0][1] * self[1][0];
        let s1 = self[0][0] * self[1][2] - self[0][2] * self[1][0];
        let s2 = self[0][0] * self[1][3] - self[0][3] * self[1][0];
        let s3 = self[0][1] * self[1][2] - self[0][2] * self[1][1];
        let s4 = self[0][1] * self[1][3] - self[0][3] * self[1][1];
        let s5 = self[0][2] * self[1][3] - self[0][3] * self[1][2];
        let c5 = self[2][2] * self[3][3] - self[2][3] * self[3][2];
        let c4 = self[2][1] * self[3][3] - self[2][3] * self[3][1];
        let c3 = self[2][1] * self[3][2] - self[2][2] * self[3][1];
        let c2 = self[2][0] * self[3][3] - self[2][3] * self[3][0];
        let c1 = self[2][0] * self[3][2] - self[2][2] * self[3][0];
        let c0 = self[2][0] * self[3][1] - self[2][1] * self[3][0];
        TensorRank2([
            TensorRank1([
                self[1][1] * c5 - self[1][2] * c4 + self[1][3] * c3,
                self[0][2] * c4 - self[0][1] * c5 - self[0][3] * c3,
                self[3][1] * s5 - self[3][2] * s4 + self[3][3] * s3,
                self[2][2] * s4 - self[2][1] * s5 - self[2][3] * s3
            ]),
            TensorRank1([
                self[1][2] * c2 - self[1][0] * c5 - self[1][3] * c1,
                self[0][0] * c5 - self[0][2] * c2 + self[0][3] * c1,
                self[3][2] * s2 - self[3][0] * s5 - self[3][3] * s1,
                self[2][0] * s5 - self[2][2] * s2 + self[2][3] * s1
            ]),
            TensorRank1([
                self[1][0] * c4 - self[1][1] * c2 + self[1][3] * c0,
                self[0][1] * c2 - self[0][0] * c4 - self[0][3] * c0,
                self[3][0] * s4 - self[3][1] * s2 + self[3][3] * s0,
                self[2][1] * s2 - self[2][0] * s4 - self[2][3] * s0
            ]),
            TensorRank1([
                self[1][1] * c1 - self[1][0] * c3 - self[1][2] * c0,
                self[0][0] * c3 - self[0][1] * c1 + self[0][2] * c0,
                self[3][1] * s1 - self[3][0] * s3 - self[3][2] * s0,
                self[2][0] * s3 - self[2][1] * s1 + self[2][2] * s0
            ])
        ]) / (s0 * c5 - s1 * c4 + s2 * c3 + s3 * c2 - s4 * c1 + s5 * c0)
    }
    fn inverse_and_determinant(&self) -> (Self, TensorRank0)
    {
        let s0 = self[0][0] * self[1][1] - self[0][1] * self[1][0];
        let s1 = self[0][0] * self[1][2] - self[0][2] * self[1][0];
        let s2 = self[0][0] * self[1][3] - self[0][3] * self[1][0];
        let s3 = self[0][1] * self[1][2] - self[0][2] * self[1][1];
        let s4 = self[0][1] * self[1][3] - self[0][3] * self[1][1];
        let s5 = self[0][2] * self[1][3] - self[0][3] * self[1][2];
        let c5 = self[2][2] * self[3][3] - self[2][3] * self[3][2];
        let c4 = self[2][1] * self[3][3] - self[2][3] * self[3][1];
        let c3 = self[2][1] * self[3][2] - self[2][2] * self[3][1];
        let c2 = self[2][0] * self[3][3] - self[2][3] * self[3][0];
        let c1 = self[2][0] * self[3][2] - self[2][2] * self[3][0];
        let c0 = self[2][0] * self[3][1] - self[2][1] * self[3][0];
        let determinant = s0 * c5 - s1 * c4 + s2 * c3 + s3 * c2 - s4 * c1 + s5 * c0;
        (
            TensorRank2([
                TensorRank1([
                    self[1][1] * c5 - self[1][2] * c4 + self[1][3] * c3,
                    self[0][2] * c4 - self[0][1] * c5 - self[0][3] * c3,
                    self[3][1] * s5 - self[3][2] * s4 + self[3][3] * s3,
                    self[2][2] * s4 - self[2][1] * s5 - self[2][3] * s3
                ]),
                TensorRank1([
                    self[1][2] * c2 - self[1][0] * c5 - self[1][3] * c1,
                    self[0][0] * c5 - self[0][2] * c2 + self[0][3] * c1,
                    self[3][2] * s2 - self[3][0] * s5 - self[3][3] * s1,
                    self[2][0] * s5 - self[2][2] * s2 + self[2][3] * s1
                ]),
                TensorRank1([
                    self[1][0] * c4 - self[1][1] * c2 + self[1][3] * c0,
                    self[0][1] * c2 - self[0][0] * c4 - self[0][3] * c0,
                    self[3][0] * s4 - self[3][1] * s2 + self[3][3] * s0,
                    self[2][1] * s2 - self[2][0] * s4 - self[2][3] * s0
                ]),
                TensorRank1([
                    self[1][1] * c1 - self[1][0] * c3 - self[1][2] * c0,
                    self[0][0] * c3 - self[0][1] * c1 + self[0][2] * c0,
                    self[3][1] * s1 - self[3][0] * s3 - self[3][2] * s0,
                    self[2][0] * s3 - self[2][1] * s1 + self[2][2] * s0
                ])
            ]) / determinant,
            determinant
        )
    }
    fn inverse_transpose(&self) -> Self
    {
        let s0 = self[0][0] * self[1][1] - self[0][1] * self[1][0];
        let s1 = self[0][0] * self[1][2] - self[0][2] * self[1][0];
        let s2 = self[0][0] * self[1][3] - self[0][3] * self[1][0];
        let s3 = self[0][1] * self[1][2] - self[0][2] * self[1][1];
        let s4 = self[0][1] * self[1][3] - self[0][3] * self[1][1];
        let s5 = self[0][2] * self[1][3] - self[0][3] * self[1][2];
        let c5 = self[2][2] * self[3][3] - self[2][3] * self[3][2];
        let c4 = self[2][1] * self[3][3] - self[2][3] * self[3][1];
        let c3 = self[2][1] * self[3][2] - self[2][2] * self[3][1];
        let c2 = self[2][0] * self[3][3] - self[2][3] * self[3][0];
        let c1 = self[2][0] * self[3][2] - self[2][2] * self[3][0];
        let c0 = self[2][0] * self[3][1] - self[2][1] * self[3][0];
        Self([
            TensorRank1([
                self[1][1] * c5 - self[1][2] * c4 + self[1][3] * c3,
                self[1][2] * c2 - self[1][0] * c5 - self[1][3] * c1,
                self[1][0] * c4 - self[1][1] * c2 + self[1][3] * c0,
                self[1][1] * c1 - self[1][0] * c3 - self[1][2] * c0,
            ]),
            TensorRank1([
                self[0][2] * c4 - self[0][1] * c5 - self[0][3] * c3,
                self[0][0] * c5 - self[0][2] * c2 + self[0][3] * c1,
                self[0][1] * c2 - self[0][0] * c4 - self[0][3] * c0,
                self[0][0] * c3 - self[0][1] * c1 + self[0][2] * c0
            ]),
            TensorRank1([
                self[3][1] * s5 - self[3][2] * s4 + self[3][3] * s3,
                self[3][2] * s2 - self[3][0] * s5 - self[3][3] * s1,
                self[3][0] * s4 - self[3][1] * s2 + self[3][3] * s0,
                self[3][1] * s1 - self[3][0] * s3 - self[3][2] * s0
            ]),
            TensorRank1([
                self[2][2] * s4 - self[2][1] * s5 - self[2][3] * s3,
                self[2][0] * s5 - self[2][2] * s2 + self[2][3] * s1,
                self[2][1] * s2 - self[2][0] * s4 - self[2][3] * s0,
                self[2][0] * s3 - self[2][1] * s1 + self[2][2] * s0
            ])
        ]) / (s0 * c5 - s1 * c4 + s2 * c3 + s3 * c2 - s4 * c1 + s5 * c0)
    }
    fn inverse_transpose_and_determinant(&self) -> (Self, TensorRank0)
    {
        let s0 = self[0][0] * self[1][1] - self[0][1] * self[1][0];
        let s1 = self[0][0] * self[1][2] - self[0][2] * self[1][0];
        let s2 = self[0][0] * self[1][3] - self[0][3] * self[1][0];
        let s3 = self[0][1] * self[1][2] - self[0][2] * self[1][1];
        let s4 = self[0][1] * self[1][3] - self[0][3] * self[1][1];
        let s5 = self[0][2] * self[1][3] - self[0][3] * self[1][2];
        let c5 = self[2][2] * self[3][3] - self[2][3] * self[3][2];
        let c4 = self[2][1] * self[3][3] - self[2][3] * self[3][1];
        let c3 = self[2][1] * self[3][2] - self[2][2] * self[3][1];
        let c2 = self[2][0] * self[3][3] - self[2][3] * self[3][0];
        let c1 = self[2][0] * self[3][2] - self[2][2] * self[3][0];
        let c0 = self[2][0] * self[3][1] - self[2][1] * self[3][0];
        let determinant = s0 * c5 - s1 * c4 + s2 * c3 + s3 * c2 - s4 * c1 + s5 * c0;
        (
            Self([
                TensorRank1([
                    self[1][1] * c5 - self[1][2] * c4 + self[1][3] * c3,
                    self[1][2] * c2 - self[1][0] * c5 - self[1][3] * c1,
                    self[1][0] * c4 - self[1][1] * c2 + self[1][3] * c0,
                    self[1][1] * c1 - self[1][0] * c3 - self[1][2] * c0,
                ]),
                TensorRank1([
                    self[0][2] * c4 - self[0][1] * c5 - self[0][3] * c3,
                    self[0][0] * c5 - self[0][2] * c2 + self[0][3] * c1,
                    self[0][1] * c2 - self[0][0] * c4 - self[0][3] * c0,
                    self[0][0] * c3 - self[0][1] * c1 + self[0][2] * c0
                ]),
                TensorRank1([
                    self[3][1] * s5 - self[3][2] * s4 + self[3][3] * s3,
                    self[3][2] * s2 - self[3][0] * s5 - self[3][3] * s1,
                    self[3][0] * s4 - self[3][1] * s2 + self[3][3] * s0,
                    self[3][1] * s1 - self[3][0] * s3 - self[3][2] * s0
                ]),
                TensorRank1([
                    self[2][2] * s4 - self[2][1] * s5 - self[2][3] * s3,
                    self[2][0] * s5 - self[2][2] * s2 + self[2][3] * s1,
                    self[2][1] * s2 - self[2][0] * s4 - self[2][3] * s0,
                    self[2][0] * s3 - self[2][1] * s1 + self[2][2] * s0
                ])
            ]) / determinant,
            determinant
        )
    }
    fn norm(&self) -> TensorRank0
    {
        (self.iter().map(|self_i|
            self_i.iter().map(|self_ij|
                self_ij.powi(2)
            ).sum::<TensorRank0>()
        ).sum::<TensorRank0>()/2.0).sqrt()
    }
    fn squared_trace(&self) -> TensorRank0
    {
        self.iter().enumerate().map(|(i, self_i)|
            self_i.iter().zip(self.iter()).map(|(self_ij, self_j)|
                self_ij * self_j[i]
            ).sum::<TensorRank0>()
        ).sum()
    }
    fn zero() -> Self
    {
        Self(std::array::from_fn(|_| TensorRank1::zero()))
    }
}

impl TensorRank2Traits<9> for TensorRank2<9>
{
    fn determinant(&self) -> TensorRank0
    {
        let (tensor_l, tensor_u) = self.lu_decomposition();
        tensor_l.iter().enumerate().zip(tensor_u.iter()).map(|((i, tensor_l_i), tensor_u_i)|
            tensor_l_i[i] * tensor_u_i[i]
        ).product()
    }
    fn deviatoric(&self) -> Self
    {
        Self::identity() * (self.trace() / -9.0) + self
    }
    fn deviatoric_and_trace(&self) -> (Self, TensorRank0)
    {
        let trace = self.trace();
        (
            Self::identity() * (trace / -9.0) + self,
            trace
        )
    }
    fn full_contraction(&self, tensor_rank_2: &Self) -> TensorRank0
    {
        self.iter().zip(tensor_rank_2.iter()).map(|(self_i, tensor_rank_2_i)|
            self_i.iter().zip(tensor_rank_2_i.iter()).map(|(self_ij, tensor_rank_2_ij)|
                self_ij * tensor_rank_2_ij
            ).sum::<TensorRank0>()
        ).sum()
    }
    fn inverse(&self) -> Self
    {
        let (tensor_l, tensor_u) = self.lu_decomposition();
        tensor_u.inverse_upper_triangular() * tensor_l.inverse_lower_triangular()
    }
    fn inverse_and_determinant(&self) -> (Self, TensorRank0)
    {
        let (tensor_l, tensor_u) = self.lu_decomposition();
        let determinant = 
            tensor_l.iter().enumerate().zip(tensor_u.iter()).map(|((i, tensor_l_i), tensor_u_i)|
                tensor_l_i[i] * tensor_u_i[i]
            ).product();
        (
            tensor_u.inverse_upper_triangular() * tensor_l.inverse_lower_triangular(),
            determinant
        )
    }
    fn inverse_transpose(&self) -> Self
    {
        self.inverse().transpose()
    }
    fn inverse_transpose_and_determinant(&self) -> (Self, TensorRank0)
    {
        let (inverse, determinant) = self.inverse_and_determinant();
        (inverse.transpose(), determinant)
    }
    fn norm(&self) -> TensorRank0
    {
        (self.iter().map(|self_i|
            self_i.iter().map(|self_ij|
                self_ij.powi(2)
            ).sum::<TensorRank0>()
        ).sum::<TensorRank0>()/2.0).sqrt()
    }
    fn squared_trace(&self) -> TensorRank0
    {
        self.iter().enumerate().map(|(i, self_i)|
            self_i.iter().zip(self.iter()).map(|(self_ij, self_j)|
                self_ij * self_j[i]
            ).sum::<TensorRank0>()
        ).sum()
    }
    fn zero() -> Self
    {
        Self(std::array::from_fn(|_| TensorRank1::zero()))
    }
}

// consider implementing things like transpose and inverse for &TensorRank2 since they are copy anyway
//     same with scalar computations like trace, det, etc.
//     can panic!() for implementations that wouldnt make sense?
//     need to test all this!
// going to have to relax a lot of the trait bounds if do this
//     might be OK anyway since going reimplement in mechanics/ anyway
//     and can go back to more performant implementations (consider for TensorRank1 as well)
//     just need the inverse and determinant to be in trait defaults

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

impl<const D: usize> std::iter::Sum for TensorRank2<D>
{
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = Self>
    {
        let mut output = TensorRank2::zero();
        iter.for_each(|item|
            output += item
        );
        output
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
where
    Self: TensorRank2Traits<D>
{
    type Output = Self;
    fn mul(self, tensor_rank_2: Self) -> Self::Output
    {
        self.iter().map(|self_i|
            self_i.iter().zip(tensor_rank_2.iter()).map(|(self_ij, tensor_rank_2_j)|
                tensor_rank_2_j * self_ij
            ).sum()
        ).collect()
    }
}

impl<const D: usize> Mul<&Self> for TensorRank2<D>
where
    Self: TensorRank2Traits<D>
{
    type Output = Self;
    fn mul(self, tensor_rank_2: &Self) -> Self::Output
    {
        self.iter().map(|self_i|
            self_i.iter().zip(tensor_rank_2.iter()).map(|(self_ij, tensor_rank_2_j)|
                tensor_rank_2_j * self_ij
            ).sum()
        ).collect()
    }
}

impl<const D: usize> Mul<TensorRank2<D>> for &TensorRank2<D>
where
    TensorRank2<D>: TensorRank2Traits<D>
{
    type Output = TensorRank2<D>;
    fn mul(self, tensor_rank_2: TensorRank2<D>) -> Self::Output
    {
        self.iter().map(|self_i|
            self_i.iter().zip(tensor_rank_2.iter()).map(|(self_ij, tensor_rank_2_j)|
                tensor_rank_2_j * self_ij
            ).sum()
        ).collect()
    }
}

impl<const D: usize> Mul for &TensorRank2<D>
where
    TensorRank2<D>: TensorRank2Traits<D>
{
    type Output = TensorRank2<D>;
    fn mul(self, tensor_rank_2: &TensorRank2<D>) -> Self::Output
    {
        self.iter().map(|self_i|
            self_i.iter().zip(tensor_rank_2.iter()).map(|(self_ij, tensor_rank_2_j)|
                tensor_rank_2_j * self_ij
            ).sum()
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
