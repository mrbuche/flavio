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
    Mul
};

use super::
{
    TensorRank0,
    TensorRank1Trait,
    TensorRank2,
    list::
    {
        TensorRank2List,
        TensorRank2ListTrait
    }
};

/// A 2D list of *d*-dimensional tensors of rank 2.
///
/// `D` is the dimension, `I`, `J` are the configurations, `W` is the list length.
pub struct TensorRank2List2D<const D: usize, const I: usize, const J: usize, const W: usize>
(
    [TensorRank2List<D, I, J, W>; W]
);

/// Inherent implementation of [`TensorRank2List2D`].
impl<const D: usize, const I: usize, const J: usize, const W: usize> TensorRank2List2D<D, I, J, W>
{
    /// Returns an iterator.
    ///
    /// The iterator yields all items from start to end. [Read more](https://doc.rust-lang.org/std/iter/)
    pub fn iter(&self) -> impl Iterator<Item = &TensorRank2List<D, I, J, W>>
    {
        self.0.iter()
    }
    /// Returns an iterator that allows modifying each value.
    ///
    /// The iterator yields all items from start to end. [Read more](https://doc.rust-lang.org/std/iter/)
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut TensorRank2List<D, I, J, W>>
    {
        self.0.iter_mut()
    }
}

/// Required methods for 2D rank-2 tensor lists.
pub trait TensorRank2List2DTrait<const D: usize, const W: usize>
{
    /// Returns the 2D rank-2 tensor list as an array.
    fn as_array(&self) -> [[[[TensorRank0; D]; D]; W]; W];
    /// Returns a list of rank-2 tensors given an array.
    fn new(array: [[[[TensorRank0; D]; D]; W]; W]) -> Self;
    /// Returns a list of rank-2 zero tensors.
    fn zero() -> Self;
}

/// Implementation of [`TensorRank2List2DTrait`] for [`TensorRank2List2D`].
impl<const D: usize, const I: usize, const J: usize, const W: usize> TensorRank2List2DTrait<D, W> for TensorRank2List2D<D, I, J, W>
{
    fn as_array(&self) -> [[[[TensorRank0; D]; D]; W]; W]
    {
        let mut array = [[[[0.0; D]; D]; W]; W];
        array.iter_mut()
        .zip(self.iter())
        .for_each(|(entry, array_entry)|
            entry.iter_mut()
            .zip(array_entry.iter())
            .for_each(|(entry_rank_2, tensor_rank_2)|
                entry_rank_2.iter_mut()
                .zip(tensor_rank_2.iter())
                .for_each(|(entry_rank_1, tensor_rank_1)|
                    *entry_rank_1 = tensor_rank_1.as_array()
                )
            )
        );
        array
    }
    fn new(array: [[[[TensorRank0; D]; D]; W]; W]) -> Self
    {
        array.iter().map(|array_i|
            TensorRank2List::new(*array_i)
        ).collect()
    }
    fn zero() -> Self
    {
        Self(std::array::from_fn(|_| TensorRank2List::zero()))
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> FromIterator<TensorRank2List<D, I, J, W>> for TensorRank2List2D<D, I, J, W>
{
    fn from_iter<Ii: IntoIterator<Item=TensorRank2List<D, I, J, W>>>(into_iterator: Ii) -> Self
    {
        let mut tensor_rank_2_list_2d = Self::zero();
        tensor_rank_2_list_2d.iter_mut().zip(into_iterator).for_each(|(tensor_rank_2_list, entry)|
            *tensor_rank_2_list = entry
        );
        tensor_rank_2_list_2d
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> Index<usize> for TensorRank2List2D<D, I, J, W>
{
    type Output = TensorRank2List<D, I, J, W>;
    fn index(&self, index: usize) -> &Self::Output
    {
        &self.0[index]
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> IndexMut<usize> for TensorRank2List2D<D, I, J, W>
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output
    {
        &mut self.0[index]
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> std::iter::Sum for TensorRank2List2D<D, I, J, W>
{
    fn sum<Ii>(iter: Ii) -> Self
    where
        Ii: Iterator<Item = Self>
    {
        let mut output = TensorRank2List2D::zero();
        iter.for_each(|item|
            output += item
        );
        output
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> Div<TensorRank0> for TensorRank2List2D<D, I, J, W>
{
    type Output = Self;
    fn div(mut self, tensor_rank_0: TensorRank0) -> Self::Output
    {
        self /= &tensor_rank_0;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> Div<TensorRank0> for &TensorRank2List2D<D, I, J, W>
{
    type Output = TensorRank2List2D<D, I, J, W>;
    fn div(self, tensor_rank_0: TensorRank0) -> Self::Output
    {
        self.iter().map(|self_i|
            self_i / tensor_rank_0
        ).collect()
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> Div<&TensorRank0> for TensorRank2List2D<D, I, J, W>
{
    type Output = Self;
    fn div(mut self, tensor_rank_0: &TensorRank0) -> Self::Output
    {
        self /= tensor_rank_0;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> Div<&TensorRank0> for &TensorRank2List2D<D, I, J, W>
{
    type Output = TensorRank2List2D<D, I, J, W>;
    fn div(self, tensor_rank_0: &TensorRank0) -> Self::Output
    {
        self.iter().map(|self_i|
            self_i / tensor_rank_0
        ).collect()
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> DivAssign<TensorRank0> for TensorRank2List2D<D, I, J, W>
{
    fn div_assign(&mut self, tensor_rank_0: TensorRank0)
    {
        self.iter_mut().for_each(|entry|
            *entry /= &tensor_rank_0
        );
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> DivAssign<&TensorRank0> for TensorRank2List2D<D, I, J, W>
{
    fn div_assign(&mut self, tensor_rank_0: &TensorRank0)
    {
        self.iter_mut().for_each(|entry|
            *entry /= tensor_rank_0
        );
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize, const W: usize> Mul<TensorRank2<D, J, K>> for TensorRank2List2D<D, I, J, W>
{
    type Output = TensorRank2List2D<D, I, K, W>;
    fn mul(self, tensor_rank_2: TensorRank2<D, J, K>) -> Self::Output
    {
        self.iter().map(|self_entry|
            self_entry.iter().map(|self_tensor_rank_2|
                self_tensor_rank_2 * &tensor_rank_2
            ).collect()
        ).collect()
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize, const W: usize> Mul<&TensorRank2<D, J, K>> for TensorRank2List2D<D, I, J, W>
{
    type Output = TensorRank2List2D<D, I, K, W>;
    fn mul(self, tensor_rank_2: &TensorRank2<D, J, K>) -> Self::Output
    {
        self.iter().map(|self_entry|
            self_entry.iter().map(|self_tensor_rank_2|
                self_tensor_rank_2 * tensor_rank_2
            ).collect()
        ).collect()
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> Add for TensorRank2List2D<D, I, J, W>
{
    type Output = Self;
    fn add(mut self, tensor_rank_2_list_2d: Self) -> Self::Output
    {
        self += tensor_rank_2_list_2d;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> Add<&Self> for TensorRank2List2D<D, I, J, W>
{
    type Output = Self;
    fn add(mut self, tensor_rank_2_list_2d: &Self) -> Self::Output
    {
        self += tensor_rank_2_list_2d;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> Add<TensorRank2List2D<D, I, J, W>> for &TensorRank2List2D<D, I, J, W>
{
    type Output = TensorRank2List2D<D, I, J, W>;
    fn add(self, mut tensor_rank_2_list_2d: TensorRank2List2D<D, I, J, W>) -> Self::Output
    {
        tensor_rank_2_list_2d += self;
        tensor_rank_2_list_2d
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> AddAssign for TensorRank2List2D<D, I, J, W>
{
    fn add_assign(&mut self, tensor_rank_2_list_2d: Self)
    {
        self.iter_mut().zip(tensor_rank_2_list_2d.iter())
        .for_each(|(self_entry, tensor_rank_2_list)|
            *self_entry += tensor_rank_2_list
        );
    }
}

impl<const D: usize, const I: usize, const J: usize, const W: usize> AddAssign<&Self> for TensorRank2List2D<D, I, J, W>
{
    fn add_assign(&mut self, tensor_rank_2_list_2d: &Self)
    {
        self.iter_mut().zip(tensor_rank_2_list_2d.iter())
        .for_each(|(self_entry, tensor_rank_2_list)|
            *self_entry += tensor_rank_2_list
        );
    }
}