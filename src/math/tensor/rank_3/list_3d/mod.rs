#[cfg(test)]
pub mod test;

use super::
{
    TensorRank0,
    list_2d::
    {
        TensorRank3List2D,
        TensorRank3List2DTrait
    }
};

type MakeClippyHappy<const D: usize> = [[[TensorRank0; D]; D]; D];

/// A 3D list of *d*-dimensional tensors of rank 3.
///
/// `D` is the dimension, `I`, `J`, `K` are the configurations `W` and `X` are the list lengths.
pub struct TensorRank3List3D<const D: usize, const I: usize, const J: usize, const K: usize, const W: usize, const X: usize>
(
    [TensorRank3List2D<D, I, J, K, W>; X]
);

/// Inherent implementation of [`TensorRank3List3D`].
impl<const D: usize, const I: usize, const J: usize, const K: usize, const W: usize, const X: usize> TensorRank3List3D<D, I, J, K, W, X>
{
    /// Returns an iterator.
    ///
    /// The iterator yields all items from start to end. [Read more](https://doc.rust-lang.org/std/iter/)
    pub fn iter(&self) -> impl Iterator<Item = &TensorRank3List2D<D, I, J, K, W>>
    {
        self.0.iter()
    }
    /// Returns an iterator that allows modifying each value.
    ///
    /// The iterator yields all items from start to end. [Read more](https://doc.rust-lang.org/std/iter/)
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut TensorRank3List2D<D, I, J, K, W>>
    {
        self.0.iter_mut()
    }
}

/// Required methods for 3D rank-3 tensor lists.
pub trait TensorRank3List3DTrait<const D: usize, const W: usize, const X: usize>
{
    /// Returns a list of rank-3 tensors given an array.
    fn new(array: [[[MakeClippyHappy<D>; W]; W]; X]) -> Self;
    /// Returns a list of rank-3 zero tensors.
    fn zero() -> Self;
}

/// Implementation of [`TensorRank3List3DTrait`] for [`TensorRank3List3D`].
impl<const D: usize, const I: usize, const J: usize, const K: usize, const W: usize, const X: usize> TensorRank3List3DTrait<D, W, X> for TensorRank3List3D<D, I, J, K, W, X>
{
    fn new(array: [[[MakeClippyHappy<D>; W]; W]; X]) -> Self
    {
        array.iter().map(|array_i|
            TensorRank3List2D::new(*array_i)
        ).collect()
    }
    fn zero() -> Self
    {
        Self(std::array::from_fn(|_| TensorRank3List2D::zero()))
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize, const W: usize, const X: usize> FromIterator<TensorRank3List2D<D, I, J, K, W>> for TensorRank3List3D<D, I, J, K, W, X>
{
    fn from_iter<Ii: IntoIterator<Item=TensorRank3List2D<D, I, J, K, W>>>(into_iterator: Ii) -> Self
    {
        let mut tensor_rank_3_list_3d = Self::zero();
        tensor_rank_3_list_3d.iter_mut().zip(into_iterator).for_each(|(tensor_rank_3_list_2d, entry)|
            *tensor_rank_3_list_2d = entry
        );
        tensor_rank_3_list_3d
    }
}