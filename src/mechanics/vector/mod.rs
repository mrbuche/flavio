#[cfg(test)]
mod test;

use crate::math::
{
    TensorRank1,
    TensorRank1Trait
};

pub struct Vector<const D: usize, const I: usize>
(
    TensorRank1<D>
);

pub trait VectorTrait<const D: usize>
{}

// traits should be independent of I and J
// conversion methods can explicitly state output I/J types since defined by name

// using 0 for reference
// 1 for current
// 2 for intermediate (like F_2)
// then do things like:
// impl Mul<Tensor<D, J, K> for Tensor<D, I, J>
// type Output = Tensor<D, I, K>
// fn inverse(&self) -> Tensor<D, J, I> // for Tensor<D, I, J>
// fn from_dyad(vector_a: Vector<D, I>, vector_b: Vector<D, J>) -> Tensor<D, I, J>
// and also
// impl TensorRank2Trait<D> for Tensor<D, I, J>
