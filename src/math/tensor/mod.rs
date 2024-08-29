pub mod rank_0;
pub mod rank_1;
pub mod rank_2;
pub mod rank_3;
pub mod rank_4;

use rank_0::TensorRank0;

/// A value-to-value conversion that does not consume the input value.
///
/// This is as opposed to [`Into`](https://doc.rust-lang.org/std/convert/trait.Into.html), which consumes the input value.
pub trait Convert<T> {
    /// Converts this type into the (usually inferred) input type.
    fn convert(&self) -> T;
}

/// Implements common methods for tensors.
pub trait Tensor {
    /// Returns the tensor norm.
    fn normy(&self) -> TensorRank0;
    /// Returns the zero tensor.
    fn zeroy() -> Self;
}

/// Implements common methods for lists of tensors.
pub trait Tensors {
    type Item;
    /// Returns an iterator that allows modifying each value.
    ///
    /// The iterator yields all items from start to end. [Read more](https://doc.rust-lang.org/std/iter/)
    fn iter_muty(&mut self) -> impl Iterator<Item = &mut Self::Item>;
    /// Returns a list of zero tensors.
    fn zeroy() -> Self;
}
