#[cfg(test)]
pub mod test;

#[cfg(test)]
use super::super::{test::TensorError, Tensor};

use super::{super::Tensors, list::TensorRank3List, TensorRank0};
use std::{
    array::from_fn,
    fmt::{Display, Formatter, Result},
};

/// A 2D list of *d*-dimensional tensors of rank 3.
///
/// `D` is the dimension, `I`, `J`, `K` are the configurations `W` and `X` are the list lengths.
pub struct TensorRank3List2D<
    const D: usize,
    const I: usize,
    const J: usize,
    const K: usize,
    const W: usize,
    const X: usize,
>([TensorRank3List<D, I, J, K, W>; X]);

impl<
        const D: usize,
        const I: usize,
        const J: usize,
        const K: usize,
        const W: usize,
        const X: usize,
    > Display for TensorRank3List2D<D, I, J, K, W, X>
{
    fn fmt(&self, _f: &mut Formatter) -> Result {
        Ok(())
    }
}

#[cfg(test)]
impl<
        const D: usize,
        const I: usize,
        const J: usize,
        const K: usize,
        const W: usize,
        const X: usize,
    > TensorError for TensorRank3List2D<D, I, J, K, W, X>
{
    fn error(
        &self,
        comparator: &Self,
        tol_abs: &TensorRank0,
        tol_rel: &TensorRank0,
    ) -> Option<usize> {
        let error_count = self
            .iter()
            .zip(comparator.iter())
            .map(|(self_a, comparator_a)| {
                self_a
                    .iter()
                    .zip(comparator_a.iter())
                    .map(|(self_ab, comparator_ab)| {
                        self_ab
                            .iter()
                            .zip(comparator_ab.iter())
                            .map(|(self_ab_i, comparator_ab_i)| {
                                self_ab_i
                                    .iter()
                                    .zip(comparator_ab_i.iter())
                                    .map(|(self_ab_ij, comparator_ab_ij)| {
                                        self_ab_ij
                                            .iter()
                                            .zip(comparator_ab_ij.iter())
                                            .filter(|(&self_ab_ijk, &comparator_ab_ijk)| {
                                                &(self_ab_ijk - comparator_ab_ijk).abs() >= tol_abs
                                                    && &(self_ab_ijk / comparator_ab_ijk - 1.0)
                                                        .abs()
                                                        >= tol_rel
                                            })
                                            .count()
                                    })
                                    .sum::<usize>()
                            })
                            .sum::<usize>()
                    })
                    .sum::<usize>()
            })
            .sum();
        if error_count > 0 {
            Some(error_count)
        } else {
            None
        }
    }
    fn error_fd(&self, comparator: &Self, epsilon: &TensorRank0) -> Option<usize> {
        let error_count = self
            .iter()
            .zip(comparator.iter())
            .map(|(self_a, comparator_a)| {
                self_a
                    .iter()
                    .zip(comparator_a.iter())
                    .map(|(self_ab, comparator_ab)| {
                        self_ab
                            .iter()
                            .zip(comparator_ab.iter())
                            .map(|(self_ab_i, comparator_ab_i)| {
                                self_ab_i
                                    .iter()
                                    .zip(comparator_ab_i.iter())
                                    .map(|(self_ab_ij, comparator_ab_ij)| {
                                        self_ab_ij
                                            .iter()
                                            .zip(comparator_ab_ij.iter())
                                            .filter(|(&self_ab_ijk, &comparator_ab_ijk)| {
                                                &(self_ab_ijk / comparator_ab_ijk - 1.0).abs()
                                                    >= epsilon
                                                    && (&self_ab_ijk.abs() >= epsilon
                                                        || &comparator_ab_ijk.abs() >= epsilon)
                                            })
                                            .count()
                                    })
                                    .sum::<usize>()
                            })
                            .sum::<usize>()
                    })
                    .sum::<usize>()
            })
            .sum();
        if error_count > 0 {
            Some(error_count)
        } else {
            None
        }
    }
}

impl<
        const D: usize,
        const I: usize,
        const J: usize,
        const K: usize,
        const W: usize,
        const X: usize,
    > Tensors for TensorRank3List2D<D, I, J, K, W, X>
{
    type Array = [[[[[TensorRank0; D]; D]; D]; W]; X];
    type Item = TensorRank3List<D, I, J, K, W>;
    fn as_array(&self) -> Self::Array {
        let mut array = [[[[[0.0; D]; D]; D]; W]; X];
        array
            .iter_mut()
            .zip(self.iter())
            .for_each(|(entry_rank_4_list, tensor_rank_4_list)| {
                *entry_rank_4_list = tensor_rank_4_list.as_array()
            });
        array
    }
    fn copy(&self) -> Self {
        self.iter().map(|entry| entry.copy()).collect()
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
            .map(|array_i| TensorRank3List::new(*array_i))
            .collect()
    }
    fn zero() -> Self {
        Self(from_fn(|_| TensorRank3List::zero()))
    }
}

impl<
        const D: usize,
        const I: usize,
        const J: usize,
        const K: usize,
        const W: usize,
        const X: usize,
    > FromIterator<TensorRank3List<D, I, J, K, W>> for TensorRank3List2D<D, I, J, K, W, X>
{
    fn from_iter<Ii: IntoIterator<Item = TensorRank3List<D, I, J, K, W>>>(
        into_iterator: Ii,
    ) -> Self {
        let mut tensor_rank_3_list_2d = Self::zero();
        tensor_rank_3_list_2d
            .iter_mut()
            .zip(into_iterator)
            .for_each(|(tensor_rank_3_list, entry)| *tensor_rank_3_list = entry);
        tensor_rank_3_list_2d
    }
}
