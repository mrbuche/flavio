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
        TensorRank1Trait
    },
    rank_2::TensorRank2
};

/// A *d*-dimensional tensor of rank 3.
///
/// `D` is the dimension, `I`, `J`, `K` are the configurations.
pub struct TensorRank3<const D: usize, const I: usize, const J: usize, const K: usize>
(
    [TensorRank2<D, J, K>; D]
);

/// Inherent implementation of [`TensorRank3`].
impl<const D: usize, const I: usize, const J: usize, const K: usize> TensorRank3<D, I, J, K>
{
    /// Returns an iterator.
    ///
    /// The iterator yields all items from start to end. [Read more](https://doc.rust-lang.org/std/iter/)
    pub fn iter(&self) -> impl Iterator<Item=&TensorRank2<D, J, K>>
    {
        self.0.iter()
    }
    /// Returns an iterator that allows modifying each value.
    ///
    /// The iterator yields all items from start to end. [Read more](https://doc.rust-lang.org/std/iter/)
    pub fn iter_mut(&mut self) -> impl Iterator<Item=&mut TensorRank2<D, J, K>>
    {
        self.0.iter_mut()
    }
    /// Returns the rank-3 zero tensor.
    pub fn zero() -> Self
    {
        Self(std::array::from_fn(|_| TensorRank2::zero()))
    }
}

/// Required methods for rank-3 tensors.
pub trait TensorRank3Trait<const D: usize>
{
    /// Returns a rank-3 tensor given an array.
    fn new(array: [[[TensorRank0; D]; D]; D]) -> Self;
}

/// Implementation of [`TensorRank3Trait`] for [`TensorRank3`].
impl<const D: usize, const I: usize, const J: usize, const K: usize> TensorRank3Trait<D> for TensorRank3<D, I, J, K>
{
    fn new(array: [[[TensorRank0; D]; D]; D]) -> Self
    {
        array.iter().map(|array_i|
            array_i.iter().map(|array_ij|
                TensorRank1::new(*array_ij)
            ).collect()
        ).collect()
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> FromIterator<TensorRank2<D, J, K>> for TensorRank3<D, I, J, K>
{
    fn from_iter<Ii: IntoIterator<Item=TensorRank2<D, J, K>>>(into_iterator: Ii) -> Self
    {
        let mut tensor_rank_3 = Self::zero();
        tensor_rank_3.iter_mut().zip(into_iterator).for_each(|(tensor_rank_3_i, value_i)|
            *tensor_rank_3_i = value_i
        );
        tensor_rank_3
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> Index<usize> for TensorRank3<D, I, J, K>
{
    type Output = TensorRank2<D, J, K>;
    fn index(&self, index: usize) -> &Self::Output
    {
        &self.0[index]
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> IndexMut<usize> for TensorRank3<D, I, J, K>
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output
    {
        &mut self.0[index]
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> Div<TensorRank0> for TensorRank3<D, I, J, K>
{
    type Output = Self;
    fn div(mut self, tensor_rank_0: TensorRank0) -> Self::Output
    {
        self /= &tensor_rank_0;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> Div<&TensorRank0> for TensorRank3<D, I, J, K>
{
    type Output = Self;
    fn div(mut self, tensor_rank_0: &TensorRank0) -> Self::Output
    {
        self /= tensor_rank_0;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> DivAssign<TensorRank0> for TensorRank3<D, I, J, K>
{
    fn div_assign(&mut self, tensor_rank_0: TensorRank0)
    {
        self.iter_mut().for_each(|self_i|
            *self_i /= &tensor_rank_0
        );
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> DivAssign<&TensorRank0> for TensorRank3<D, I, J, K>
{
    fn div_assign(&mut self, tensor_rank_0: &TensorRank0)
    {
        self.iter_mut().for_each(|self_i|
            *self_i /= tensor_rank_0
        );
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> Mul<TensorRank0> for TensorRank3<D, I, J, K>
{
    type Output = Self;
    fn mul(mut self, tensor_rank_0: TensorRank0) -> Self::Output
    {
        self *= &tensor_rank_0;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> Mul<&TensorRank0> for TensorRank3<D, I, J, K>
{
    type Output = Self;
    fn mul(mut self, tensor_rank_0: &TensorRank0) -> Self::Output
    {
        self *= tensor_rank_0;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> MulAssign<TensorRank0> for TensorRank3<D, I, J, K>
{
    fn mul_assign(&mut self, tensor_rank_0: TensorRank0)
    {
        self.iter_mut().for_each(|self_i|
            *self_i *= &tensor_rank_0
        );
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> MulAssign<&TensorRank0> for TensorRank3<D, I, J, K>
{
    fn mul_assign(&mut self, tensor_rank_0: &TensorRank0)
    {
        self.iter_mut().for_each(|self_i|
            *self_i *= tensor_rank_0
        );
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> Add for TensorRank3<D, I, J, K>
{
    type Output = Self;
    fn add(mut self, tensor_rank_3: Self) -> Self::Output
    {
        self += tensor_rank_3;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> Add<&Self> for TensorRank3<D, I, J, K>
{
    type Output = Self;
    fn add(mut self, tensor_rank_3: &Self) -> Self::Output
    {
        self += tensor_rank_3;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> Add<TensorRank3<D, I, J, K>> for &TensorRank3<D, I, J, K>
{
    type Output = TensorRank3<D, I, J, K>;
    fn add(self, mut tensor_rank_3: TensorRank3<D, I, J, K>) -> Self::Output
    {
        tensor_rank_3 += self;
        tensor_rank_3
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> AddAssign for TensorRank3<D, I, J, K>
{
    fn add_assign(&mut self, tensor_rank_3: Self)
    {
        self.iter_mut().zip(tensor_rank_3.iter()).for_each(|(self_i, tensor_rank_3_i)|
            *self_i += tensor_rank_3_i
        );
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> AddAssign<&Self> for TensorRank3<D, I, J, K>
{
    fn add_assign(&mut self, tensor_rank_3: &Self)
    {
        self.iter_mut().zip(tensor_rank_3.iter()).for_each(|(self_i, tensor_rank_3_i)|
            *self_i += tensor_rank_3_i
        );
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> Sub for TensorRank3<D, I, J, K>
{
    type Output = Self;
    fn sub(mut self, tensor_rank_3: Self) -> Self::Output
    {
        self -= tensor_rank_3;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> Sub<&Self> for TensorRank3<D, I, J, K>
{
    type Output = Self;
    fn sub(mut self, tensor_rank_3: &Self) -> Self::Output
    {
        self -= tensor_rank_3;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> SubAssign for TensorRank3<D, I, J, K>
{
    fn sub_assign(&mut self, tensor_rank_3: Self)
    {
        self.iter_mut().zip(tensor_rank_3.iter()).for_each(|(self_i, tensor_rank_3_i)|
            *self_i -= tensor_rank_3_i
        );
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize> SubAssign<&Self> for TensorRank3<D, I, J, K>
{
    fn sub_assign(&mut self, tensor_rank_3: &Self)
    {
        self.iter_mut().zip(tensor_rank_3.iter()).for_each(|(self_i, tensor_rank_3_i)|
            *self_i -= tensor_rank_3_i
        );
    }
}