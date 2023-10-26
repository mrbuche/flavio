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

// impl Mul<Tensor<D, J, K> for Tensor<D, I, J>
// type Output = Tensor<D, I, K>

// fn inverse(&self) -> Tensor<D, J, I> // for Tensor<D, I, J>

// fn from_dyad(vector_a: Vector<D, I>, vector_b: Vector<D, J>) -> Tensor<D, I, J>
