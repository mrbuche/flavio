#[cfg(test)]
mod test;

use std::ops::
{
    Index,
    IndexMut,
    Mul,
    MulAssign
};

use super::
{
    TensorRank0
};

/// A list of tensors of rank 0 (a list of scalars).
///
/// `W` is the list length.
pub struct TensorRank0List<const W: usize>
(
    [TensorRank0; W]
);

impl<const W: usize> TensorRank0List<W>
{
    /// Returns an iterator.
    ///
    /// The iterator yields all items from start to end. [Read more](https://doc.rust-lang.org/std/iter/)
    pub fn iter(&self) -> impl Iterator<Item = &TensorRank0>
    {
        self.0.iter()
    }
    /// Returns an iterator that allows modifying each value.
    ///
    /// The iterator yields all items from start to end. [Read more](https://doc.rust-lang.org/std/iter/)
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut TensorRank0>
    {
        self.0.iter_mut()
    }
}

/// Required methods for rank-0 tensor lists.
pub trait TensorRank0ListTrait<const W: usize>
{
    /// Returns a list of rank-0 tensors given an array.
    fn new(array: [TensorRank0; W]) -> Self;
    /// Returns a list of rank-0 zero tensors.
    fn zero() -> Self;
}

/// Implementation of [`TensorRank0ListTrait`] for [`TensorRank0List`].
impl<const W: usize> TensorRank0ListTrait<W> for TensorRank0List<W>
{
    fn new(array: [TensorRank0; W]) -> Self
    {
        Self(array)
    }
    fn zero() -> Self
    {
        Self([0.0; W])
    }
}

impl<const W: usize> FromIterator<TensorRank0> for TensorRank0List<W>
{
    fn from_iter<Ii: IntoIterator<Item=TensorRank0>>(into_iterator: Ii) -> Self
    {
        let mut tensor_rank_0_list = Self::zero();
        tensor_rank_0_list.iter_mut().zip(into_iterator).for_each(|(tensor_rank_0, entry)|
            *tensor_rank_0 = entry
        );
        tensor_rank_0_list
    }
}

impl<const W: usize> Index<usize> for TensorRank0List<W>
{
    type Output = TensorRank0;
    fn index(&self, index: usize) -> &Self::Output
    {
        &self.0[index]
    }
}

impl<const W: usize> IndexMut<usize> for TensorRank0List<W>
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output
    {
        &mut self.0[index]
    }
}

impl<const W: usize> Mul<TensorRank0> for TensorRank0List<W>
{
    type Output = Self;
    fn mul(mut self, tensor_rank_0: TensorRank0) -> Self::Output
    {
        self *= tensor_rank_0;
        self
    }
}

impl<const W: usize> Mul<TensorRank0> for &TensorRank0List<W>
{
    type Output = TensorRank0List<W>;
    fn mul(self, tensor_rank_0: TensorRank0) -> Self::Output
    {
        self.iter().map(|self_i|
            self_i * tensor_rank_0
        ).collect()
    }
}

impl<const W: usize> Mul<&TensorRank0> for TensorRank0List<W>
{
    type Output = Self;
    fn mul(mut self, tensor_rank_0: &TensorRank0) -> Self::Output
    {
        self *= tensor_rank_0;
        self
    }
}

impl<const W: usize> Mul<&TensorRank0> for &TensorRank0List<W>
{
    type Output = TensorRank0List<W>;
    fn mul(self, tensor_rank_0: &TensorRank0) -> Self::Output
    {
        self.iter().map(|self_i|
            self_i * tensor_rank_0
        ).collect()
    }
}

impl<const W: usize> MulAssign<TensorRank0> for TensorRank0List<W>
{
    fn mul_assign(&mut self, tensor_rank_0: TensorRank0)
    {
        self.iter_mut().for_each(|self_i|
            *self_i *= &tensor_rank_0
        );
    }
}

impl<const W: usize> MulAssign<&TensorRank0> for TensorRank0List<W>
{
    fn mul_assign(&mut self, tensor_rank_0: &TensorRank0)
    {
        self.iter_mut().for_each(|self_i|
            *self_i *= tensor_rank_0
        );
    }
}

impl<const W: usize> Mul for TensorRank0List<W>
{
    type Output = TensorRank0;
    fn mul(self, tensor_rank_0_list: Self) -> Self::Output
    {
        self.iter()
        .zip(tensor_rank_0_list.iter())
        .map(|(self_entry, tensor_rank_0_list_entry)|
            self_entry * tensor_rank_0_list_entry
        ).sum()
    }
}

impl<const W: usize> Mul<&Self> for TensorRank0List<W>
{
    type Output = TensorRank0;
    fn mul(self, tensor_rank_0_list: &Self) -> Self::Output
    {
        self.iter()
        .zip(tensor_rank_0_list.iter())
        .map(|(self_entry, tensor_rank_0_list_entry)|
            self_entry * tensor_rank_0_list_entry
        ).sum()
    }
}

impl<const W: usize> Mul<TensorRank0List<W>> for &TensorRank0List<W>
{
    type Output = TensorRank0;
    fn mul(self, tensor_rank_0_list: TensorRank0List<W>) -> Self::Output
    {
        self.iter()
        .zip(tensor_rank_0_list.iter())
        .map(|(self_entry, tensor_rank_0_list_entry)|
            self_entry * tensor_rank_0_list_entry
        ).sum()
    }
}

impl<const W: usize> Mul for &TensorRank0List<W>
{
    type Output = TensorRank0;
    fn mul(self, tensor_rank_0_list: Self) -> Self::Output
    {
        self.iter()
        .zip(tensor_rank_0_list.iter())
        .map(|(self_entry, tensor_rank_0_list_entry)|
            self_entry * tensor_rank_0_list_entry
        ).sum()
    }
}
