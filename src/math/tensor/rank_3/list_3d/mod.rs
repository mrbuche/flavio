#[cfg(test)]
pub mod test;

#[cfg(test)]
use super::super::{test::TensorError, Tensor};

use super::{super::Tensors, list_2d::TensorRank3List2D, TensorRank0};
use std::{
    array::from_fn,
    fmt::{Display, Formatter, Result},
};

type MakeClippyHappy<const D: usize> = [[[TensorRank0; D]; D]; D];

/// A 3D list of *d*-dimensional tensors of rank 3.
///
/// `D` is the dimension, `I`, `J`, `K` are the configurations `W`, `X`, and `Y` are the list lengths.
pub struct TensorRank3List3D<
    const D: usize,
    const I: usize,
    const J: usize,
    const K: usize,
    const W: usize,
    const X: usize,
    const Y: usize,
>([TensorRank3List2D<D, I, J, K, W, X>; Y]);

impl<
        const D: usize,
        const I: usize,
        const J: usize,
        const K: usize,
        const W: usize,
        const X: usize,
        const Y: usize,
    > Display for TensorRank3List3D<D, I, J, K, W, X, Y>
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
        const Y: usize,
    > TensorError for TensorRank3List3D<D, I, J, K, W, X, Y>
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
                            .map(|(self_abc, comparator_abc)| {
                                self_abc
                                    .iter()
                                    .zip(comparator_abc.iter())
                                    .map(|(self_abc_i, comparator_abc_i)| {
                                        self_abc_i
                                            .iter()
                                            .zip(comparator_abc_i.iter())
                                            .map(|(self_abc_ij, comparator_abc_ij)| {
                                                self_abc_ij
                                                    .iter()
                                                    .zip(comparator_abc_ij.iter())
                                                    .filter(
                                                        |(&self_abc_ijk, &comparator_abc_ijk)| {
                                                            &(self_abc_ijk - comparator_abc_ijk)
                                                                .abs()
                                                                >= tol_abs
                                                                && &(self_abc_ijk
                                                                    / comparator_abc_ijk
                                                                    - 1.0)
                                                                    .abs()
                                                                    >= tol_rel
                                                        },
                                                    )
                                                    .count()
                                            })
                                            .sum::<usize>()
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
    fn error_fd(&self, comparator: &Self, epsilon: &TensorRank0) -> Option<(bool, usize)> {
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
                            .map(|(self_abc, comparator_abc)| {
                                self_abc
                                    .iter()
                                    .zip(comparator_abc.iter())
                                    .map(|(self_abc_i, comparator_abc_i)| {
                                        self_abc_i
                                            .iter()
                                            .zip(comparator_abc_i.iter())
                                            .map(|(self_abc_ij, comparator_abc_ij)| {
                                                self_abc_ij
                                                    .iter()
                                                    .zip(comparator_abc_ij.iter())
                                                    .filter(
                                                        |(&self_abc_ijk, &comparator_abc_ijk)| {
                                                            &(self_abc_ijk / comparator_abc_ijk
                                                                - 1.0)
                                                                .abs()
                                                                >= epsilon
                                                                && (&self_abc_ijk.abs() >= epsilon
                                                                    || &comparator_abc_ijk.abs()
                                                                        >= epsilon)
                                                        },
                                                    )
                                                    .count()
                                            })
                                            .sum::<usize>()
                                    })
                                    .sum::<usize>()
                            })
                            .sum::<usize>()
                    })
                    .sum::<usize>()
            })
            .sum();
        if error_count > 0 {
            Some((true, error_count))
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
        const Y: usize,
    > Tensors for TensorRank3List3D<D, I, J, K, W, X, Y>
{
    type Array = [[[MakeClippyHappy<D>; W]; X]; Y];
    type Item = TensorRank3List2D<D, I, J, K, W, X>;
    fn as_array(&self) -> Self::Array {
        let mut array = [[[[[[0.0; D]; D]; D]; W]; X]; Y];
        array.iter_mut().zip(self.iter()).for_each(
            |(entry_rank_4_list_2d, tensor_rank_4_list_2d)| {
                *entry_rank_4_list_2d = tensor_rank_4_list_2d.as_array()
            },
        );
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
            .map(|array_i| TensorRank3List2D::new(*array_i))
            .collect()
    }
    fn zero() -> Self {
        Self(from_fn(|_| TensorRank3List2D::zero()))
    }
}

impl<
        const D: usize,
        const I: usize,
        const J: usize,
        const K: usize,
        const W: usize,
        const X: usize,
        const Y: usize,
    > FromIterator<TensorRank3List2D<D, I, J, K, W, X>> for TensorRank3List3D<D, I, J, K, W, X, Y>
{
    fn from_iter<Ii: IntoIterator<Item = TensorRank3List2D<D, I, J, K, W, X>>>(
        into_iterator: Ii,
    ) -> Self {
        let mut tensor_rank_3_list_3d = Self::zero();
        tensor_rank_3_list_3d
            .iter_mut()
            .zip(into_iterator)
            .for_each(|(tensor_rank_3_list_2d, entry)| *tensor_rank_3_list_2d = entry);
        tensor_rank_3_list_3d
    }
}
