#[cfg(test)]
mod test;

use crate::math::
{
    TensorRank2,
    TensorRank2Traits
};

pub struct Tensor<const D: usize, const I: usize, const J: usize>
(
    pub TensorRank2<D>
);

pub trait TensorTraits<const D: usize>
{}