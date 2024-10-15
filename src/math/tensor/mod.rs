#[cfg(test)]
pub mod test;

pub mod rank_0;
pub mod rank_1;
pub mod rank_2;
pub mod rank_3;
pub mod rank_4;

use rank_0::TensorRank0;
use std::{
    fmt,
    ops::{Add, AddAssign, Div, Mul, Sub},
};

/// A value-to-value conversion that does not consume the input value.
///
/// This is as opposed to [`Into`](https://doc.rust-lang.org/std/convert/trait.Into.html), which consumes the input value.
pub trait Convert<T> {
    /// Converts this type into the (usually inferred) input type.
    fn convert(&self) -> T;
}

/// Common methods for tensors.
pub trait Tensor
where
    for<'a> Self: Add<Self, Output = Self>
        + Add<&'a Self, Output = Self>
        + AddAssign
        + Div<TensorRank0, Output = Self>
        + Mul<TensorRank0, Output = Self>
        + Sub<Self, Output = Self>
        + Sub<&'a Self, Output = Self>
        + Sized
        + fmt::Display,
{
    type Array;
    type Item;
    /// Returns the tensor as an array.
    fn as_array(&self) -> Self::Array;
    /// Returns a copy.
    ///
    /// This method was implemented instead of the Copy trait to avoid unintended copy creations.
    fn copy(&self) -> Self;
    /// Returns the identity tensor.
    fn identity() -> Self;
    /// Checks whether the tensor is the zero tensor.
    fn is_zero(&self) -> bool;
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
    fn norm_squared(&self) -> TensorRank0;
    /// Returns the tensor normalized.
    fn normalized(&self) -> Self;
    /// Returns the zero tensor.
    fn zero() -> Self;
}

/// Common methods for lists of tensors.
pub trait Tensors {
    type Array;
    type Item;
    /// Returns the list of tensors as an array.
    fn as_array(&self) -> Self::Array;
    /// Returns a copy.
    ///
    /// This method was implemented instead of the Copy trait to avoid unintended copy creations.
    fn copy(&self) -> Self;
    /// Returns an iterator.
    ///
    /// The iterator yields all items from start to end. [Read more](https://doc.rust-lang.org/std/iter/)
    fn iter(&self) -> impl Iterator<Item = &Self::Item>;
    /// Returns an iterator that allows modifying each value.
    ///
    /// The iterator yields all items from start to end. [Read more](https://doc.rust-lang.org/std/iter/)
    fn iter_mut(&mut self) -> impl Iterator<Item = &mut Self::Item>;
    /// Returns a list of tensors given an array.
    fn new(array: Self::Array) -> Self;
    /// Returns a list of zero tensors.
    fn zero() -> Self;
}
