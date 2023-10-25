#[cfg(test)]
mod test;

use crate::math::
{
    TensorRank2,
    TensorRank2Trait
};

pub struct Tensor<const D: usize, const I: usize, const J: usize>
(
    TensorRank2<D>
);

pub trait TensorTrait<const D: usize>
{}