pub mod rank_0;
pub mod rank_1;
pub mod rank_2;
pub mod rank_3;
pub mod rank_4;

use rank_0::TensorRank0;
use std::ops::{Add, Div, Mul, Sub};

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
        + Div<TensorRank0, Output = Self>
        + Mul<TensorRank0, Output = Self>
        + Sub<&'a Self, Output = Self>
        + Sized,
{
    type Array;
    /// Returns the tensor as an array.
    fn as_array(&self) -> Self::Array;
    /// Returns a copy.
    ///
    /// This method was implemented instead of the Copy trait to avoid unintended copy creations.
    fn copy(&self) -> Self;
    /// Returns a tensor given an array.
    fn new(array: Self::Array) -> Self;
    /// Returns the tensor norm.
    fn norm(&self) -> TensorRank0;
    /// Returns the zero tensor.
    fn zero() -> Self;
}

/// Common methods for lists of tensors.
pub trait Tensors {
    type Array;
    type Item;
    /// Returns the list of tensors as an array.
    fn as_array(&self) -> Self::Array;
    /// Returns the sum of the full dot product of each tensor in each list.
    fn dot(&self, tensors: &Self) -> TensorRank0;
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
