#[cfg(test)]
mod test;

pub mod list;
pub mod list_2d;

use std::{
    array::from_fn,
    cmp::Ordering,
    ops::{Add, AddAssign, Div, DivAssign, Index, IndexMut, Mul, MulAssign, Sub, SubAssign},
};

use super::{
    super::write_tensor_rank_0,
    rank_0::TensorRank0,
    rank_1::{list::TensorRank1List, TensorRank1},
    rank_4::TensorRank4,
    Convert, Tensor, Tensors,
};
use list_2d::TensorRank2List2D;
use std::fmt;

/// A *d*-dimensional tensor of rank 2.
///
/// `D` is the dimension, `I`, `J` are the configurations.
#[derive(Debug)]
pub struct TensorRank2<const D: usize, const I: usize, const J: usize>(pub [TensorRank1<D, J>; D]);

/// Display implementation for rank-2 tensors.
impl<const D: usize, const I: usize, const J: usize> fmt::Display for TensorRank2<D, I, J> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "\x1B[s")?;
        write!(f, "[[")?;
        self.iter().enumerate().for_each(|(i, row)| {
            row.iter()
                .for_each(|entry| write_tensor_rank_0(f, entry).unwrap());
            if i + 1 < D {
                writeln!(f, "\x1B[2D],").unwrap();
                write!(f, "\x1B[u").unwrap();
                write!(f, "\x1B[{}B[ ", i + 1).unwrap();
            }
        });
        write!(f, "\x1B[2D]]")?;
        Ok(())
    }
}

/// Inherent implementation of [`TensorRank2`].
impl<const D: usize, const I: usize, const J: usize> TensorRank2<D, I, J> {
    /// Returns the determinant of the rank-2 tensor.
    pub fn determinant(&self) -> TensorRank0 {
        if D == 2 {
            self[0][0] * self[1][1] - self[0][1] * self[1][0]
        } else if D == 3 {
            let c_00 = self[1][1] * self[2][2] - self[1][2] * self[2][1];
            let c_10 = self[1][2] * self[2][0] - self[1][0] * self[2][2];
            let c_20 = self[1][0] * self[2][1] - self[1][1] * self[2][0];
            self[0][0] * c_00 + self[0][1] * c_10 + self[0][2] * c_20
        } else if D == 4 {
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
        } else {
            let (tensor_l, tensor_u) = self.lu_decomposition();
            tensor_l
                .iter()
                .enumerate()
                .zip(tensor_u.iter())
                .map(|((i, tensor_l_i), tensor_u_i)| tensor_l_i[i] * tensor_u_i[i])
                .product()
        }
    }
    /// Returns the deviatoric component of the rank-2 tensor.
    pub fn deviatoric(&self) -> Self {
        Self::identity() * (self.trace() / -(D as TensorRank0)) + self
    }
    /// Returns the deviatoric component and trace of the rank-2 tensor.
    pub fn deviatoric_and_trace(&self) -> (Self, TensorRank0) {
        let trace = self.trace();
        (
            Self::identity() * (trace / -(D as TensorRank0)) + self,
            trace,
        )
    }
    /// Returns a rank-2 tensor constructed from a dyad of the given vectors.
    pub fn dyad(vector_a: &TensorRank1<D, I>, vector_b: &TensorRank1<D, J>) -> Self {
        vector_a
            .iter()
            .map(|vector_a_i| {
                vector_b
                    .iter()
                    .map(|vector_b_j| vector_a_i * vector_b_j)
                    .collect()
            })
            .collect()
    }
    /// Returns the full contraction with another rank-2 tensor.
    pub fn full_contraction(&self, tensor_rank_2: &Self) -> TensorRank0 {
        self.iter()
            .zip(tensor_rank_2.iter())
            .map(|(self_i, tensor_rank_2_i)| {
                self_i
                    .iter()
                    .zip(tensor_rank_2_i.iter())
                    .map(|(self_ij, tensor_rank_2_ij)| self_ij * tensor_rank_2_ij)
                    .sum::<TensorRank0>()
            })
            .sum()
    }
    /// Returns the inverse of the rank-2 tensor.
    pub fn inverse(&self) -> TensorRank2<D, J, I> {
        if D == 2 {
            let mut adjugate = TensorRank2::<D, J, I>::zero();
            adjugate[0][0] = self[1][1];
            adjugate[0][1] = -self[0][1];
            adjugate[1][0] = -self[1][0];
            adjugate[1][1] = self[0][0];
            adjugate / self.determinant()
        } else if D == 3 {
            let mut adjugate = TensorRank2::<D, J, I>::zero();
            let c_00 = self[1][1] * self[2][2] - self[1][2] * self[2][1];
            let c_10 = self[1][2] * self[2][0] - self[1][0] * self[2][2];
            let c_20 = self[1][0] * self[2][1] - self[1][1] * self[2][0];
            adjugate[0][0] = c_00;
            adjugate[0][1] = self[0][2] * self[2][1] - self[0][1] * self[2][2];
            adjugate[0][2] = self[0][1] * self[1][2] - self[0][2] * self[1][1];
            adjugate[1][0] = c_10;
            adjugate[1][1] = self[0][0] * self[2][2] - self[0][2] * self[2][0];
            adjugate[1][2] = self[0][2] * self[1][0] - self[0][0] * self[1][2];
            adjugate[2][0] = c_20;
            adjugate[2][1] = self[0][1] * self[2][0] - self[0][0] * self[2][1];
            adjugate[2][2] = self[0][0] * self[1][1] - self[0][1] * self[1][0];
            adjugate / (self[0][0] * c_00 + self[0][1] * c_10 + self[0][2] * c_20)
        } else if D == 4 {
            let mut adjugate = TensorRank2::<D, J, I>::zero();
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
            adjugate[0][0] = self[1][1] * c5 - self[1][2] * c4 + self[1][3] * c3;
            adjugate[0][1] = self[0][2] * c4 - self[0][1] * c5 - self[0][3] * c3;
            adjugate[0][2] = self[3][1] * s5 - self[3][2] * s4 + self[3][3] * s3;
            adjugate[0][3] = self[2][2] * s4 - self[2][1] * s5 - self[2][3] * s3;
            adjugate[1][0] = self[1][2] * c2 - self[1][0] * c5 - self[1][3] * c1;
            adjugate[1][1] = self[0][0] * c5 - self[0][2] * c2 + self[0][3] * c1;
            adjugate[1][2] = self[3][2] * s2 - self[3][0] * s5 - self[3][3] * s1;
            adjugate[1][3] = self[2][0] * s5 - self[2][2] * s2 + self[2][3] * s1;
            adjugate[2][0] = self[1][0] * c4 - self[1][1] * c2 + self[1][3] * c0;
            adjugate[2][1] = self[0][1] * c2 - self[0][0] * c4 - self[0][3] * c0;
            adjugate[2][2] = self[3][0] * s4 - self[3][1] * s2 + self[3][3] * s0;
            adjugate[2][3] = self[2][1] * s2 - self[2][0] * s4 - self[2][3] * s0;
            adjugate[3][0] = self[1][1] * c1 - self[1][0] * c3 - self[1][2] * c0;
            adjugate[3][1] = self[0][0] * c3 - self[0][1] * c1 + self[0][2] * c0;
            adjugate[3][2] = self[3][1] * s1 - self[3][0] * s3 - self[3][2] * s0;
            adjugate[3][3] = self[2][0] * s3 - self[2][1] * s1 + self[2][2] * s0;
            adjugate / (s0 * c5 - s1 * c4 + s2 * c3 + s3 * c2 - s4 * c1 + s5 * c0)
        } else {
            let (tensor_l_inverse, tensor_u_inverse) = self.lu_decomposition_inverse();
            tensor_u_inverse * tensor_l_inverse
        }
    }
    /// Returns the inverse and determinant of the rank-2 tensor.
    pub fn inverse_and_determinant(&self) -> (TensorRank2<D, J, I>, TensorRank0) {
        if D == 2 {
            let mut adjugate = TensorRank2::<D, J, I>::zero();
            adjugate[0][0] = self[1][1];
            adjugate[0][1] = -self[0][1];
            adjugate[1][0] = -self[1][0];
            adjugate[1][1] = self[0][0];
            let determinant = self.determinant();
            (adjugate / determinant, determinant)
        } else if D == 3 {
            let mut adjugate = TensorRank2::<D, J, I>::zero();
            let c_00 = self[1][1] * self[2][2] - self[1][2] * self[2][1];
            let c_10 = self[1][2] * self[2][0] - self[1][0] * self[2][2];
            let c_20 = self[1][0] * self[2][1] - self[1][1] * self[2][0];
            let determinant = self[0][0] * c_00 + self[0][1] * c_10 + self[0][2] * c_20;
            adjugate[0][0] = c_00;
            adjugate[0][1] = self[0][2] * self[2][1] - self[0][1] * self[2][2];
            adjugate[0][2] = self[0][1] * self[1][2] - self[0][2] * self[1][1];
            adjugate[1][0] = c_10;
            adjugate[1][1] = self[0][0] * self[2][2] - self[0][2] * self[2][0];
            adjugate[1][2] = self[0][2] * self[1][0] - self[0][0] * self[1][2];
            adjugate[2][0] = c_20;
            adjugate[2][1] = self[0][1] * self[2][0] - self[0][0] * self[2][1];
            adjugate[2][2] = self[0][0] * self[1][1] - self[0][1] * self[1][0];
            (adjugate / determinant, determinant)
        } else if D == 4 {
            let mut adjugate = TensorRank2::<D, J, I>::zero();
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
            adjugate[0][0] = self[1][1] * c5 - self[1][2] * c4 + self[1][3] * c3;
            adjugate[0][1] = self[0][2] * c4 - self[0][1] * c5 - self[0][3] * c3;
            adjugate[0][2] = self[3][1] * s5 - self[3][2] * s4 + self[3][3] * s3;
            adjugate[0][3] = self[2][2] * s4 - self[2][1] * s5 - self[2][3] * s3;
            adjugate[1][0] = self[1][2] * c2 - self[1][0] * c5 - self[1][3] * c1;
            adjugate[1][1] = self[0][0] * c5 - self[0][2] * c2 + self[0][3] * c1;
            adjugate[1][2] = self[3][2] * s2 - self[3][0] * s5 - self[3][3] * s1;
            adjugate[1][3] = self[2][0] * s5 - self[2][2] * s2 + self[2][3] * s1;
            adjugate[2][0] = self[1][0] * c4 - self[1][1] * c2 + self[1][3] * c0;
            adjugate[2][1] = self[0][1] * c2 - self[0][0] * c4 - self[0][3] * c0;
            adjugate[2][2] = self[3][0] * s4 - self[3][1] * s2 + self[3][3] * s0;
            adjugate[2][3] = self[2][1] * s2 - self[2][0] * s4 - self[2][3] * s0;
            adjugate[3][0] = self[1][1] * c1 - self[1][0] * c3 - self[1][2] * c0;
            adjugate[3][1] = self[0][0] * c3 - self[0][1] * c1 + self[0][2] * c0;
            adjugate[3][2] = self[3][1] * s1 - self[3][0] * s3 - self[3][2] * s0;
            adjugate[3][3] = self[2][0] * s3 - self[2][1] * s1 + self[2][2] * s0;
            (adjugate / determinant, determinant)
        } else {
            panic!()
        }
    }
    /// Returns the inverse transpose of the rank-2 tensor.
    pub fn inverse_transpose(&self) -> Self {
        if D == 2 {
            let mut adjugate_transpose = TensorRank2::<D, I, J>::zero();
            adjugate_transpose[0][0] = self[1][1];
            adjugate_transpose[0][1] = -self[1][0];
            adjugate_transpose[1][0] = -self[0][1];
            adjugate_transpose[1][1] = self[0][0];
            adjugate_transpose / self.determinant()
        } else if D == 3 {
            let mut adjugate_transpose = TensorRank2::<D, I, J>::zero();
            let c_00 = self[1][1] * self[2][2] - self[1][2] * self[2][1];
            let c_10 = self[1][2] * self[2][0] - self[1][0] * self[2][2];
            let c_20 = self[1][0] * self[2][1] - self[1][1] * self[2][0];
            adjugate_transpose[0][0] = c_00;
            adjugate_transpose[1][0] = self[0][2] * self[2][1] - self[0][1] * self[2][2];
            adjugate_transpose[2][0] = self[0][1] * self[1][2] - self[0][2] * self[1][1];
            adjugate_transpose[0][1] = c_10;
            adjugate_transpose[1][1] = self[0][0] * self[2][2] - self[0][2] * self[2][0];
            adjugate_transpose[2][1] = self[0][2] * self[1][0] - self[0][0] * self[1][2];
            adjugate_transpose[0][2] = c_20;
            adjugate_transpose[1][2] = self[0][1] * self[2][0] - self[0][0] * self[2][1];
            adjugate_transpose[2][2] = self[0][0] * self[1][1] - self[0][1] * self[1][0];
            adjugate_transpose / (self[0][0] * c_00 + self[0][1] * c_10 + self[0][2] * c_20)
        } else if D == 4 {
            let mut adjugate_transpose = TensorRank2::<D, I, J>::zero();
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
            adjugate_transpose[0][0] = self[1][1] * c5 - self[1][2] * c4 + self[1][3] * c3;
            adjugate_transpose[1][0] = self[0][2] * c4 - self[0][1] * c5 - self[0][3] * c3;
            adjugate_transpose[2][0] = self[3][1] * s5 - self[3][2] * s4 + self[3][3] * s3;
            adjugate_transpose[3][0] = self[2][2] * s4 - self[2][1] * s5 - self[2][3] * s3;
            adjugate_transpose[0][1] = self[1][2] * c2 - self[1][0] * c5 - self[1][3] * c1;
            adjugate_transpose[1][1] = self[0][0] * c5 - self[0][2] * c2 + self[0][3] * c1;
            adjugate_transpose[2][1] = self[3][2] * s2 - self[3][0] * s5 - self[3][3] * s1;
            adjugate_transpose[3][1] = self[2][0] * s5 - self[2][2] * s2 + self[2][3] * s1;
            adjugate_transpose[0][2] = self[1][0] * c4 - self[1][1] * c2 + self[1][3] * c0;
            adjugate_transpose[1][2] = self[0][1] * c2 - self[0][0] * c4 - self[0][3] * c0;
            adjugate_transpose[2][2] = self[3][0] * s4 - self[3][1] * s2 + self[3][3] * s0;
            adjugate_transpose[3][2] = self[2][1] * s2 - self[2][0] * s4 - self[2][3] * s0;
            adjugate_transpose[0][3] = self[1][1] * c1 - self[1][0] * c3 - self[1][2] * c0;
            adjugate_transpose[1][3] = self[0][0] * c3 - self[0][1] * c1 + self[0][2] * c0;
            adjugate_transpose[2][3] = self[3][1] * s1 - self[3][0] * s3 - self[3][2] * s0;
            adjugate_transpose[3][3] = self[2][0] * s3 - self[2][1] * s1 + self[2][2] * s0;
            adjugate_transpose / (s0 * c5 - s1 * c4 + s2 * c3 + s3 * c2 - s4 * c1 + s5 * c0)
        } else {
            self.inverse().transpose()
        }
    }
    /// Returns the inverse transpose and determinant of the rank-2 tensor.
    pub fn inverse_transpose_and_determinant(&self) -> (Self, TensorRank0) {
        if D == 2 {
            let mut adjugate_transpose = TensorRank2::<D, I, J>::zero();
            adjugate_transpose[0][0] = self[1][1];
            adjugate_transpose[0][1] = -self[1][0];
            adjugate_transpose[1][0] = -self[0][1];
            adjugate_transpose[1][1] = self[0][0];
            let determinant = self.determinant();
            (adjugate_transpose / determinant, determinant)
        } else if D == 3 {
            let mut adjugate_transpose = TensorRank2::<D, I, J>::zero();
            let c_00 = self[1][1] * self[2][2] - self[1][2] * self[2][1];
            let c_10 = self[1][2] * self[2][0] - self[1][0] * self[2][2];
            let c_20 = self[1][0] * self[2][1] - self[1][1] * self[2][0];
            let determinant = self[0][0] * c_00 + self[0][1] * c_10 + self[0][2] * c_20;
            adjugate_transpose[0][0] = c_00;
            adjugate_transpose[1][0] = self[0][2] * self[2][1] - self[0][1] * self[2][2];
            adjugate_transpose[2][0] = self[0][1] * self[1][2] - self[0][2] * self[1][1];
            adjugate_transpose[0][1] = c_10;
            adjugate_transpose[1][1] = self[0][0] * self[2][2] - self[0][2] * self[2][0];
            adjugate_transpose[2][1] = self[0][2] * self[1][0] - self[0][0] * self[1][2];
            adjugate_transpose[0][2] = c_20;
            adjugate_transpose[1][2] = self[0][1] * self[2][0] - self[0][0] * self[2][1];
            adjugate_transpose[2][2] = self[0][0] * self[1][1] - self[0][1] * self[1][0];
            (adjugate_transpose / determinant, determinant)
        } else if D == 4 {
            let mut adjugate_transpose = TensorRank2::<D, I, J>::zero();
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
            adjugate_transpose[0][0] = self[1][1] * c5 - self[1][2] * c4 + self[1][3] * c3;
            adjugate_transpose[1][0] = self[0][2] * c4 - self[0][1] * c5 - self[0][3] * c3;
            adjugate_transpose[2][0] = self[3][1] * s5 - self[3][2] * s4 + self[3][3] * s3;
            adjugate_transpose[3][0] = self[2][2] * s4 - self[2][1] * s5 - self[2][3] * s3;
            adjugate_transpose[0][1] = self[1][2] * c2 - self[1][0] * c5 - self[1][3] * c1;
            adjugate_transpose[1][1] = self[0][0] * c5 - self[0][2] * c2 + self[0][3] * c1;
            adjugate_transpose[2][1] = self[3][2] * s2 - self[3][0] * s5 - self[3][3] * s1;
            adjugate_transpose[3][1] = self[2][0] * s5 - self[2][2] * s2 + self[2][3] * s1;
            adjugate_transpose[0][2] = self[1][0] * c4 - self[1][1] * c2 + self[1][3] * c0;
            adjugate_transpose[1][2] = self[0][1] * c2 - self[0][0] * c4 - self[0][3] * c0;
            adjugate_transpose[2][2] = self[3][0] * s4 - self[3][1] * s2 + self[3][3] * s0;
            adjugate_transpose[3][2] = self[2][1] * s2 - self[2][0] * s4 - self[2][3] * s0;
            adjugate_transpose[0][3] = self[1][1] * c1 - self[1][0] * c3 - self[1][2] * c0;
            adjugate_transpose[1][3] = self[0][0] * c3 - self[0][1] * c1 + self[0][2] * c0;
            adjugate_transpose[2][3] = self[3][1] * s1 - self[3][0] * s3 - self[3][2] * s0;
            adjugate_transpose[3][3] = self[2][0] * s3 - self[2][1] * s1 + self[2][2] * s0;
            (adjugate_transpose / determinant, determinant)
        } else {
            panic!()
        }
    }
    /// Checks whether the tensor is a diagonal tensor.
    pub fn is_diagonal(&self) -> bool {
        self.iter()
            .enumerate()
            .map(|(i, self_i)| {
                self_i
                    .iter()
                    .enumerate()
                    .map(|(j, self_ij)| (self_ij == &0.0) as u8 * (i != j) as u8)
                    .sum::<u8>()
            })
            .sum::<u8>()
            == (D * D - D) as u8
    }
    /// Checks whether the tensor is the identity tensor.
    pub fn is_identity(&self) -> bool {
        self.iter()
            .enumerate()
            .map(|(i, self_i)| {
                self_i
                    .iter()
                    .enumerate()
                    .map(|(j, self_ij)| (self_ij == &((i == j) as u8 as f64)) as u8)
                    .sum::<u8>()
            })
            .sum::<u8>()
            == (D * D) as u8
    }
    /// Returns the LU decomposition of the rank-2 tensor.
    pub fn lu_decomposition(&self) -> (TensorRank2<D, I, 88>, TensorRank2<D, 88, J>) {
        let mut tensor_l = TensorRank2::zero();
        let mut tensor_u = TensorRank2::zero();
        for i in 0..D {
            for j in 0..D {
                if j >= i {
                    tensor_l[j][i] = self[j][i];
                    for k in 0..i {
                        tensor_l[j][i] -= tensor_l[j][k] * tensor_u[k][i];
                    }
                }
            }
            for j in 0..D {
                match j.cmp(&i) {
                    Ordering::Equal => {
                        tensor_u[i][j] = 1.0;
                    }
                    Ordering::Greater => {
                        tensor_u[i][j] = self[i][j] / tensor_l[i][i];
                        for k in 0..i {
                            tensor_u[i][j] -= (tensor_l[i][k] * tensor_u[k][j]) / tensor_l[i][i];
                        }
                    }
                    Ordering::Less => (),
                }
            }
        }
        (tensor_l, tensor_u)
    }
    /// Returns the inverse of the LU decomposition of the rank-2 tensor.
    pub fn lu_decomposition_inverse(&self) -> (TensorRank2<D, 88, I>, TensorRank2<D, J, 88>) {
        let mut tensor_l = TensorRank2::zero();
        let mut tensor_u = TensorRank2::zero();
        for i in 0..D {
            for j in 0..D {
                if j >= i {
                    tensor_l[j][i] = self[j][i];
                    for k in 0..i {
                        tensor_l[j][i] -= tensor_l[j][k] * tensor_u[k][i];
                    }
                }
            }
            for j in 0..D {
                match j.cmp(&i) {
                    Ordering::Equal => {
                        tensor_u[i][j] = 1.0;
                    }
                    Ordering::Greater => {
                        tensor_u[i][j] = self[i][j] / tensor_l[i][i];
                        for k in 0..i {
                            tensor_u[i][j] -= (tensor_l[i][k] * tensor_u[k][j]) / tensor_l[i][i];
                        }
                    }
                    Ordering::Less => (),
                }
            }
        }
        let mut sum;
        for i in 0..D {
            tensor_l[i][i] = 1.0 / tensor_l[i][i];
            for j in 0..i {
                sum = 0.0;
                for k in j..i {
                    sum += tensor_l[i][k] * tensor_l[k][j];
                }
                tensor_l[i][j] = -sum * tensor_l[i][i];
            }
        }
        for i in 0..D {
            tensor_u[i][i] = 1.0 / tensor_u[i][i];
            for j in 0..i {
                sum = 0.0;
                for k in j..i {
                    sum += tensor_u[j][k] * tensor_u[k][i];
                }
                tensor_u[j][i] = -sum * tensor_u[i][i];
            }
        }
        (tensor_l, tensor_u)
    }
    /// Returns the second invariant of the rank-2 tensor.
    pub fn second_invariant(&self) -> TensorRank0 {
        0.5 * (self.trace().powi(2) - self.squared_trace())
    }
    /// Returns the trace of the rank-2 tensor squared.
    pub fn squared_trace(&self) -> TensorRank0 {
        self.iter()
            .enumerate()
            .map(|(i, self_i)| {
                self_i
                    .iter()
                    .zip(self.iter())
                    .map(|(self_ij, self_j)| self_ij * self_j[i])
                    .sum::<TensorRank0>()
            })
            .sum()
    }
    /// Returns the trace of the rank-2 tensor.
    pub fn trace(&self) -> TensorRank0 {
        (0..D).map(|i| self[i][i]).sum()
    }
    /// Returns the transpose of the rank-2 tensor.
    pub fn transpose(&self) -> TensorRank2<D, J, I> {
        (0..D)
            .map(|i| (0..D).map(|j| self[j][i]).collect())
            .collect()
    }
}

/// Implementation of [`Tensor`] for [`TensorRank2`].
impl<const D: usize, const I: usize, const J: usize> Tensor for TensorRank2<D, I, J> {
    type Array = [[TensorRank0; D]; D];
    type Item = TensorRank1<D, J>;
    fn as_array(&self) -> Self::Array {
        let mut array = [[0.0; D]; D];
        array
            .iter_mut()
            .zip(self.iter())
            .for_each(|(entry, tensor_rank_1)| *entry = tensor_rank_1.as_array());
        array
    }
    fn copy(&self) -> Self {
        self.iter().map(|entry| entry.copy()).collect()
    }
    fn identity() -> Self {
        (0..D)
            .map(|i| (0..D).map(|j| ((i == j) as u8) as TensorRank0).collect())
            .collect()
    }
    fn is_zero(&self) -> bool {
        self.iter()
            .map(|self_i| {
                self_i
                    .iter()
                    .map(|self_ij| (self_ij == &0.0) as u8)
                    .sum::<u8>()
            })
            .sum::<u8>()
            == (D * D) as u8
    }
    fn iter(&self) -> impl Iterator<Item = &Self::Item> {
        self.0.iter()
    }
    fn iter_mut(&mut self) -> impl Iterator<Item = &mut Self::Item> {
        self.0.iter_mut()
    }
    fn new(array: Self::Array) -> Self {
        array
            .iter()
            .map(|array_i| TensorRank1::new(*array_i))
            .collect()
    }
    fn norm_squared(&self) -> TensorRank0 {
        self.iter()
            .map(|self_i| {
                self_i
                    .iter()
                    .map(|self_ij| self_ij.powi(2))
                    .sum::<TensorRank0>()
            })
            .sum()
    }
    fn normalized(&self) -> Self {
        self / self.norm()
    }
    fn zero() -> Self {
        Self(from_fn(|_| TensorRank1::zero()))
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize, const L: usize>
    Convert<TensorRank2<D, K, L>> for TensorRank2<D, I, J>
{
    fn convert(&self) -> TensorRank2<D, K, L> {
        self.iter()
            .map(|self_i| self_i.iter().copied().collect())
            .collect()
    }
}

impl<const D: usize> From<TensorRank2<D, 0, 0>> for TensorRank2<D, 1, 0> {
    fn from(tensor_rank_2: TensorRank2<D, 0, 0>) -> Self {
        tensor_rank_2
            .iter()
            .map(|tensor_rank_1| tensor_rank_1.into())
            .collect()
    }
}

impl<const D: usize> From<TensorRank2<D, 0, 0>> for TensorRank2<D, 1, 1> {
    fn from(tensor_rank_2: TensorRank2<D, 0, 0>) -> Self {
        tensor_rank_2
            .iter()
            .map(|tensor_rank_1| tensor_rank_1.into())
            .collect()
    }
}

impl<const D: usize> From<TensorRank2<D, 0, 1>> for TensorRank2<D, 0, 0> {
    fn from(tensor_rank_2: TensorRank2<D, 0, 1>) -> Self {
        tensor_rank_2
            .iter()
            .map(|tensor_rank_1| tensor_rank_1.into())
            .collect()
    }
}

impl<const D: usize> From<TensorRank2<D, 1, 0>> for TensorRank2<D, 0, 0> {
    fn from(tensor_rank_2: TensorRank2<D, 1, 0>) -> Self {
        tensor_rank_2
            .iter()
            .map(|tensor_rank_1| tensor_rank_1.into())
            .collect()
    }
}

impl<const D: usize> From<TensorRank2<D, 1, 1>> for TensorRank2<D, 1, 0> {
    fn from(tensor_rank_2: TensorRank2<D, 1, 1>) -> Self {
        tensor_rank_2
            .iter()
            .map(|tensor_rank_1| tensor_rank_1.into())
            .collect()
    }
}

impl<const D: usize> From<TensorRank2<D, 1, 2>> for TensorRank2<D, 1, 0> {
    fn from(tensor_rank_2: TensorRank2<D, 1, 2>) -> Self {
        tensor_rank_2
            .iter()
            .map(|tensor_rank_1| tensor_rank_1.into())
            .collect()
    }
}

impl<const D: usize, const I: usize, const J: usize> FromIterator<TensorRank1<D, J>>
    for TensorRank2<D, I, J>
{
    fn from_iter<Ii: IntoIterator<Item = TensorRank1<D, J>>>(into_iterator: Ii) -> Self {
        let mut tensor_rank_2 = Self::zero();
        tensor_rank_2
            .iter_mut()
            .zip(into_iterator)
            .for_each(|(tensor_rank_2_i, value_i)| *tensor_rank_2_i = value_i);
        tensor_rank_2
    }
}

impl<const D: usize, const I: usize, const J: usize> Index<usize> for TensorRank2<D, I, J> {
    type Output = TensorRank1<D, J>;
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl<const D: usize, const I: usize, const J: usize> IndexMut<usize> for TensorRank2<D, I, J> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}

impl<const D: usize, const I: usize, const J: usize> std::iter::Sum for TensorRank2<D, I, J> {
    fn sum<Ii>(iter: Ii) -> Self
    where
        Ii: Iterator<Item = Self>,
    {
        let mut output = Self::zero();
        iter.for_each(|item| output += item);
        output
    }
}

impl<const D: usize, const I: usize, const J: usize> Div<TensorRank0> for TensorRank2<D, I, J> {
    type Output = Self;
    fn div(mut self, tensor_rank_0: TensorRank0) -> Self::Output {
        self /= &tensor_rank_0;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize> Div<TensorRank0> for &TensorRank2<D, I, J> {
    type Output = TensorRank2<D, I, J>;
    fn div(self, tensor_rank_0: TensorRank0) -> Self::Output {
        self.iter().map(|self_i| self_i / tensor_rank_0).collect()
    }
}

impl<const D: usize, const I: usize, const J: usize> Div<&TensorRank0> for TensorRank2<D, I, J> {
    type Output = Self;
    fn div(mut self, tensor_rank_0: &TensorRank0) -> Self::Output {
        self /= tensor_rank_0;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize> Div<&TensorRank0> for &TensorRank2<D, I, J> {
    type Output = TensorRank2<D, I, J>;
    fn div(self, tensor_rank_0: &TensorRank0) -> Self::Output {
        self.iter().map(|self_i| self_i / tensor_rank_0).collect()
    }
}

impl<const D: usize, const I: usize, const J: usize> DivAssign<TensorRank0>
    for TensorRank2<D, I, J>
{
    fn div_assign(&mut self, tensor_rank_0: TensorRank0) {
        self.iter_mut().for_each(|self_i| *self_i /= &tensor_rank_0);
    }
}

impl<const D: usize, const I: usize, const J: usize> DivAssign<&TensorRank0>
    for TensorRank2<D, I, J>
{
    fn div_assign(&mut self, tensor_rank_0: &TensorRank0) {
        self.iter_mut().for_each(|self_i| *self_i /= tensor_rank_0);
    }
}

impl<const D: usize, const I: usize, const J: usize> Mul<TensorRank0> for TensorRank2<D, I, J> {
    type Output = Self;
    fn mul(mut self, tensor_rank_0: TensorRank0) -> Self::Output {
        self *= &tensor_rank_0;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize> Mul<TensorRank0> for &TensorRank2<D, I, J> {
    type Output = TensorRank2<D, I, J>;
    fn mul(self, tensor_rank_0: TensorRank0) -> Self::Output {
        self.iter().map(|self_i| self_i * tensor_rank_0).collect()
    }
}

impl<const D: usize, const I: usize, const J: usize> Mul<&TensorRank0> for TensorRank2<D, I, J> {
    type Output = Self;
    fn mul(mut self, tensor_rank_0: &TensorRank0) -> Self::Output {
        self *= tensor_rank_0;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize> Mul<&TensorRank0> for &TensorRank2<D, I, J> {
    type Output = TensorRank2<D, I, J>;
    fn mul(self, tensor_rank_0: &TensorRank0) -> Self::Output {
        self.iter().map(|self_i| self_i * tensor_rank_0).collect()
    }
}

impl<const D: usize, const I: usize, const J: usize> MulAssign<TensorRank0>
    for TensorRank2<D, I, J>
{
    fn mul_assign(&mut self, tensor_rank_0: TensorRank0) {
        self.iter_mut().for_each(|self_i| *self_i *= &tensor_rank_0);
    }
}

impl<const D: usize, const I: usize, const J: usize> MulAssign<&TensorRank0>
    for TensorRank2<D, I, J>
{
    fn mul_assign(&mut self, tensor_rank_0: &TensorRank0) {
        self.iter_mut().for_each(|self_i| *self_i *= tensor_rank_0);
    }
}

impl<const D: usize, const I: usize, const J: usize> Mul<TensorRank1<D, J>>
    for TensorRank2<D, I, J>
{
    type Output = TensorRank1<D, I>;
    fn mul(self, tensor_rank_1: TensorRank1<D, J>) -> Self::Output {
        self.iter().map(|self_i| self_i * &tensor_rank_1).collect()
    }
}

impl<const D: usize, const I: usize, const J: usize> Mul<&TensorRank1<D, J>>
    for TensorRank2<D, I, J>
{
    type Output = TensorRank1<D, I>;
    fn mul(self, tensor_rank_1: &TensorRank1<D, J>) -> Self::Output {
        self.iter().map(|self_i| self_i * tensor_rank_1).collect()
    }
}

impl<const D: usize, const I: usize, const J: usize> Mul<TensorRank1<D, J>>
    for &TensorRank2<D, I, J>
{
    type Output = TensorRank1<D, I>;
    fn mul(self, tensor_rank_1: TensorRank1<D, J>) -> Self::Output {
        self.iter().map(|self_i| self_i * &tensor_rank_1).collect()
    }
}

impl<const D: usize, const I: usize, const J: usize> Mul<&TensorRank1<D, J>>
    for &TensorRank2<D, I, J>
{
    type Output = TensorRank1<D, I>;
    fn mul(self, tensor_rank_1: &TensorRank1<D, J>) -> Self::Output {
        self.iter().map(|self_i| self_i * tensor_rank_1).collect()
    }
}

impl<const D: usize, const I: usize, const J: usize> Add for TensorRank2<D, I, J> {
    type Output = Self;
    fn add(mut self, tensor_rank_2: Self) -> Self::Output {
        self += tensor_rank_2;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize> Add<&Self> for TensorRank2<D, I, J> {
    type Output = Self;
    fn add(mut self, tensor_rank_2: &Self) -> Self::Output {
        self += tensor_rank_2;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize> Add<TensorRank2<D, I, J>>
    for &TensorRank2<D, I, J>
{
    type Output = TensorRank2<D, I, J>;
    fn add(self, mut tensor_rank_2: TensorRank2<D, I, J>) -> Self::Output {
        tensor_rank_2 += self;
        tensor_rank_2
    }
}

impl<const D: usize, const I: usize, const J: usize> AddAssign for TensorRank2<D, I, J> {
    fn add_assign(&mut self, tensor_rank_2: Self) {
        self.iter_mut()
            .zip(tensor_rank_2.iter())
            .for_each(|(self_i, tensor_rank_2_i)| *self_i += tensor_rank_2_i);
    }
}

impl<const D: usize, const I: usize, const J: usize> AddAssign<&Self> for TensorRank2<D, I, J> {
    fn add_assign(&mut self, tensor_rank_2: &Self) {
        self.iter_mut()
            .zip(tensor_rank_2.iter())
            .for_each(|(self_i, tensor_rank_2_i)| *self_i += tensor_rank_2_i);
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> Mul<TensorRank2<D, J, K>>
    for TensorRank2<D, I, J>
{
    type Output = TensorRank2<D, I, K>;
    fn mul(self, tensor_rank_2: TensorRank2<D, J, K>) -> Self::Output {
        self.iter()
            .map(|self_i| {
                self_i
                    .iter()
                    .zip(tensor_rank_2.iter())
                    .map(|(self_ij, tensor_rank_2_j)| tensor_rank_2_j * self_ij)
                    .sum()
            })
            .collect()
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> Mul<&TensorRank2<D, J, K>>
    for TensorRank2<D, I, J>
{
    type Output = TensorRank2<D, I, K>;
    fn mul(self, tensor_rank_2: &TensorRank2<D, J, K>) -> Self::Output {
        self.iter()
            .map(|self_i| {
                self_i
                    .iter()
                    .zip(tensor_rank_2.iter())
                    .map(|(self_ij, tensor_rank_2_j)| tensor_rank_2_j * self_ij)
                    .sum()
            })
            .collect()
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> Mul<TensorRank2<D, J, K>>
    for &TensorRank2<D, I, J>
{
    type Output = TensorRank2<D, I, K>;
    fn mul(self, tensor_rank_2: TensorRank2<D, J, K>) -> Self::Output {
        self.iter()
            .map(|self_i| {
                self_i
                    .iter()
                    .zip(tensor_rank_2.iter())
                    .map(|(self_ij, tensor_rank_2_j)| tensor_rank_2_j * self_ij)
                    .sum()
            })
            .collect()
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> Mul<&TensorRank2<D, J, K>>
    for &TensorRank2<D, I, J>
{
    type Output = TensorRank2<D, I, K>;
    fn mul(self, tensor_rank_2: &TensorRank2<D, J, K>) -> Self::Output {
        self.iter()
            .map(|self_i| {
                self_i
                    .iter()
                    .zip(tensor_rank_2.iter())
                    .map(|(self_ij, tensor_rank_2_j)| tensor_rank_2_j * self_ij)
                    .sum()
            })
            .collect()
    }
}

impl<const D: usize, const I: usize, const J: usize> Sub for TensorRank2<D, I, J> {
    type Output = Self;
    fn sub(mut self, tensor_rank_2: Self) -> Self::Output {
        self -= tensor_rank_2;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize> Sub<&Self> for TensorRank2<D, I, J> {
    type Output = Self;
    fn sub(mut self, tensor_rank_2: &Self) -> Self::Output {
        self -= tensor_rank_2;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize> SubAssign for TensorRank2<D, I, J> {
    fn sub_assign(&mut self, tensor_rank_2: Self) {
        self.iter_mut()
            .zip(tensor_rank_2.iter())
            .for_each(|(self_i, tensor_rank_2_i)| *self_i -= tensor_rank_2_i);
    }
}

impl<const D: usize, const I: usize, const J: usize> SubAssign<&Self> for TensorRank2<D, I, J> {
    fn sub_assign(&mut self, tensor_rank_2: &Self) {
        self.iter_mut()
            .zip(tensor_rank_2.iter())
            .for_each(|(self_i, tensor_rank_2_i)| *self_i -= tensor_rank_2_i);
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> Mul<TensorRank1List<D, J, W>>
    for TensorRank2<D, I, J>
{
    type Output = TensorRank1List<D, I, W>;
    fn mul(self, tensor_rank_1_list: TensorRank1List<D, J, W>) -> Self::Output {
        tensor_rank_1_list
            .iter()
            .map(|tensor_rank_1| &self * tensor_rank_1)
            .collect()
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> Mul<&TensorRank1List<D, J, W>>
    for TensorRank2<D, I, J>
{
    type Output = TensorRank1List<D, I, W>;
    fn mul(self, tensor_rank_1_list: &TensorRank1List<D, J, W>) -> Self::Output {
        tensor_rank_1_list
            .iter()
            .map(|tensor_rank_1| &self * tensor_rank_1)
            .collect()
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> Mul<TensorRank1List<D, J, W>>
    for &TensorRank2<D, I, J>
{
    type Output = TensorRank1List<D, I, W>;
    fn mul(self, tensor_rank_1_list: TensorRank1List<D, J, W>) -> Self::Output {
        tensor_rank_1_list
            .iter()
            .map(|tensor_rank_1| self * tensor_rank_1)
            .collect()
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> Mul<&TensorRank1List<D, J, W>>
    for &TensorRank2<D, I, J>
{
    type Output = TensorRank1List<D, I, W>;
    fn mul(self, tensor_rank_1_list: &TensorRank1List<D, J, W>) -> Self::Output {
        tensor_rank_1_list
            .iter()
            .map(|tensor_rank_1| self * tensor_rank_1)
            .collect()
    }
}

impl<
        const D: usize,
        const I: usize,
        const J: usize,
        const K: usize,
        const W: usize,
        const X: usize,
    > Mul<TensorRank2List2D<D, J, K, W, X>> for TensorRank2<D, I, J>
{
    type Output = TensorRank2List2D<D, I, K, W, X>;
    fn mul(self, tensor_rank_2_list_2d: TensorRank2List2D<D, J, K, W, X>) -> Self::Output {
        tensor_rank_2_list_2d
            .iter()
            .map(|tensor_rank_2_list_2d_entry| {
                tensor_rank_2_list_2d_entry
                    .iter()
                    .map(|tensor_rank_2| &self * tensor_rank_2)
                    .collect()
            })
            .collect()
    }
}

impl<
        const D: usize,
        const I: usize,
        const J: usize,
        const K: usize,
        const W: usize,
        const X: usize,
    > Mul<TensorRank2List2D<D, J, K, W, X>> for &TensorRank2<D, I, J>
{
    type Output = TensorRank2List2D<D, I, K, W, X>;
    fn mul(self, tensor_rank_2_list_2d: TensorRank2List2D<D, J, K, W, X>) -> Self::Output {
        tensor_rank_2_list_2d
            .iter()
            .map(|tensor_rank_2_list_2d_entry| {
                tensor_rank_2_list_2d_entry
                    .iter()
                    .map(|tensor_rank_2| self * tensor_rank_2)
                    .collect()
            })
            .collect()
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl<const I: usize, const J: usize, const K: usize, const L: usize> Div<TensorRank4<3, I, J, K, L>>
    for TensorRank2<3, I, J>
{
    type Output = TensorRank2<3, K, L>;
    fn div(self, tensor_rank_4: TensorRank4<3, I, J, K, L>) -> Self::Output {
        let mut tensor_rank_1 = TensorRank1::<9, 88>::zero();
        self.iter().enumerate().for_each(|(i, self_i)| {
            self_i
                .iter()
                .enumerate()
                .for_each(|(j, self_ij)| tensor_rank_1[3 * i + j] = *self_ij)
        });
        let mut tensor_rank_2 = TensorRank2::<9, 88, 99>::zero();
        tensor_rank_4
            .iter()
            .enumerate()
            .for_each(|(i, tensor_rank_4_i)| {
                tensor_rank_4_i
                    .iter()
                    .enumerate()
                    .for_each(|(j, tensor_rank_4_ij)| {
                        tensor_rank_4_ij
                            .iter()
                            .enumerate()
                            .for_each(|(k, tensor_rank_4_ijk)| {
                                tensor_rank_4_ijk.iter().enumerate().for_each(
                                    |(l, tensor_rank_4_ijkl)| {
                                        tensor_rank_2[3 * i + j][3 * k + l] = *tensor_rank_4_ijkl
                                    },
                                )
                            })
                    })
            });
        let output_tensor_rank_1 = tensor_rank_2.inverse() * tensor_rank_1;
        let mut output = TensorRank2::zero();
        output.iter_mut().enumerate().for_each(|(i, output_i)| {
            output_i
                .iter_mut()
                .enumerate()
                .for_each(|(j, output_ij)| *output_ij = output_tensor_rank_1[3 * i + j])
        });
        output
    }
}
