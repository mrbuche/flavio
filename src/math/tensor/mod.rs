#[cfg(test)]
pub mod test;

pub mod rank_0;
pub mod rank_1;
pub mod rank_2;
pub mod rank_3;
pub mod rank_4;

use rank_0::TensorRank0;
use std::{
    fmt::{Debug, Display},
    ops::{Add, AddAssign, Div, Mul, Sub, SubAssign},
};

/// A value-to-value conversion that does not consume the input value.
///
/// This is as opposed to [`Into`](https://doc.rust-lang.org/std/convert/trait.Into.html), which consumes the input value.
pub trait Convert<T> {
    /// Converts this type into the (usually inferred) input type.
    fn convert(&self) -> T;
}

/// Possible errors for tensors.
#[derive(Debug)]
pub enum TensorError {
    NotPositiveDefinite,
}

impl PartialEq for TensorError {
    fn eq(&self, other: &Self) -> bool {
        match self {
            Self::NotPositiveDefinite => match other {
                Self::NotPositiveDefinite => true,
            },
        }
    }
}

/// Common methods for tensors.
pub trait Tensor
where
    for<'a> Self: Sized
        + Debug
        + Display
        + Add<Self, Output = Self>
        + Add<&'a Self, Output = Self>
        + AddAssign
        + AddAssign<&'a Self>
        + Div<TensorRank0, Output = Self>
        + Mul<TensorRank0, Output = Self>
        + Sub<Self, Output = Self>
        + Sub<&'a Self, Output = Self>
        + SubAssign
        + SubAssign<&'a Self>,
    Self::Item: Tensor,
{
    type Array;
    type Item;
    /// Returns the tensor as an array.
    fn as_array(&self) -> Self::Array;
    /// Returns a copy.
    ///
    /// This method was implemented instead of the Copy trait to avoid unintended copy creations.
    fn copy(&self) -> Self;
    /// Returns the full contraction with another tensor.
    fn full_contraction(&self, tensor: &Self) -> TensorRank0 {
        self.iter()
            .zip(tensor.iter())
            .map(|(self_entry, tensor_entry)| self_entry.full_contraction(tensor_entry))
            .sum()
    }
    /// Returns a reference to the entry at the specified indices.
    fn get_at(&self, _indices: &[usize]) -> &TensorRank0 {
        panic!("Need to implement get_at() for {:?}.", self)
    }
    /// Returns a mutable reference to the entry at the specified indices.
    fn get_at_mut(&mut self, _indices: &[usize]) -> &mut TensorRank0 {
        panic!("Need to implement get_at_mut() for {:?}.", self)
    }
    /// Returns the identity tensor.
    fn identity() -> Self;
    /// Checks whether the tensor is positive-definite.
    fn is_positive_definite(&self) -> bool {
        panic!("Need to implement is_positive_definite() for {:?}.", self)
    }
    /// Checks whether the tensor is the zero tensor.
    fn is_zero(&self) -> bool {
        self.iter().filter(|entry| !entry.is_zero()).count() == 0
    }
    /// Returns an iterator.
    ///
    /// The iterator yields all items from start to end. [Read more](https://doc.rust-lang.org/std/iter/)
    fn iter(&self) -> impl Iterator<Item = &Self::Item>;
    /// Returns an iterator that allows modifying each value.
    ///
    /// The iterator yields all items from start to end. [Read more](https://doc.rust-lang.org/std/iter/)
    fn iter_mut(&mut self) -> impl Iterator<Item = &mut Self::Item>;
    /// Returns a tensor given an array.
    fn new(array: Self::Array) -> Self;
    /// Returns the tensor norm.
    fn norm(&self) -> TensorRank0 {
        self.norm_squared().sqrt()
    }
    /// Returns the tensor norm squared.
    fn norm_squared(&self) -> TensorRank0 {
        self.full_contraction(self)
    }
    /// Returns the tensor normalized.
    fn normalized(&self) -> Self {
        panic!("Need to implement normalized() for {:?}.", self)
    }
    /// Returns the zero tensor.
    fn zero() -> Self;
}

// do Vector, Matrix, and MatrixSym as the 1D and 2D Vec types
// MatrixSym so you can do Cholesky and store ~1/2 the data
// try to keep MatrixSym from impl less-efficient things that don't use things like Cholesky

// benchmark this stuff before merging!
