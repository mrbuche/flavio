#[cfg(test)]
mod test;

use super::
{
    TensorRank0,
    TensorRank1,
    TensorRank1Trait,
    list::
    {
        TensorRank1List,
        TensorRank1ListTrait
    }
};

/// A 2D list of *d*-dimensional tensors of rank 1.
///
/// `D` is the dimension, `I` is the configuration, `W` and `X` are the list lengths.
pub struct TensorRank1List2D<const D: usize, const I: usize, const W: usize, const X: usize>
(
    [TensorRank1List<D, I, W>; X]
);

/// Inherent implementation of [`TensorRank1List2D`].
impl<const D: usize, const I: usize, const W: usize, const X: usize> TensorRank1List2D<D, I, W, X>
{
    /// Returns an iterator.
    ///
    /// The iterator yields all items from start to end. [Read more](https://doc.rust-lang.org/std/iter/)
    pub fn iter(&self) -> impl Iterator<Item = &TensorRank1List<D, I, W>>
    {
        self.0.iter()
    }
    /// Returns an iterator that allows modifying each value.
    ///
    /// The iterator yields all items from start to end. [Read more](https://doc.rust-lang.org/std/iter/)
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut TensorRank1List<D, I, W>>
    {
        self.0.iter_mut()
    }
}

/// Required methods for 2D rank-2 tensor lists.
pub trait TensorRank1List2DTrait<const D: usize, const W: usize, const X: usize>
{
    /// Returns a list of rank-1 zero tensors.
    fn zero() -> Self;
}

/// Implementation of [`TensorRank1List2DTrait`] for [`TensorRank1List2D`].
impl<const D: usize, const I: usize, const W: usize, const X: usize> TensorRank1List2DTrait<D, W, X> for TensorRank1List2D<D, I, W, X>
{
    fn zero() -> Self
    {
        Self(std::array::from_fn(|_| TensorRank1List::zero()))
    }
}