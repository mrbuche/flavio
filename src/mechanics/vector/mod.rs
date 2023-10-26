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

use crate::math::
{
    TensorRank1,
    TensorRank1Trait
};

use super::
{
    scalar::Scalar
};

pub struct Vector<const D: usize, const I: usize>
(
    TensorRank1<D>
);

impl<const D: usize, const I: usize> Vector<D, I>
{
    pub fn iter(&self) -> impl Iterator<Item = &Scalar>
    {
        self.0.iter()
    }
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut Scalar>
    {
        self.0.iter_mut()
    }
}

pub trait VectorTrait<const D: usize>
{
    fn new(array: [Scalar; D]) -> Self;
    fn norm(&self) -> Scalar;
    fn to_current_config(self) -> Vector<D, 1>;
    fn to_intermediate_config(self) -> Vector<D, 2>;
    fn to_reference_config(self) -> Vector<D, 0>;
    fn zero() -> Self;
}

impl<const D: usize, const I: usize> VectorTrait<D> for Vector<D, I>
{
    fn new(array: [Scalar; D]) -> Self
    {
        Self(TensorRank1::new(array))
    }
    fn norm(&self) -> Scalar
    {
        (self * self).sqrt()
    }
    fn to_current_config(self) -> Vector<D, 1>
    {
        Vector(self.0)
    }
    fn to_intermediate_config(self) -> Vector<D, 2>
    {
        Vector(self.0)
    }
    fn to_reference_config(self) -> Vector<D, 0>
    {
        Vector(self.0)
    }
    fn zero() -> Self
    {
        Self(TensorRank1::zero())
    }
}

impl<const D: usize, const I: usize> FromIterator<Scalar> for Vector<D, I>
{
    fn from_iter<II: IntoIterator<Item=Scalar>>(into_iterator: II) -> Self
    {
        let mut vector = Self::zero();
        vector.iter_mut().zip(into_iterator).for_each(|(vector_i, value_i)|
            *vector_i = value_i
        );
        vector
    }
}

impl<const D: usize, const I: usize> Mul for Vector<D, I>
{
    type Output = Scalar;
    fn mul(self, vector: Self) -> Self::Output
    {
        self.0 * vector.0
    }
}

impl<const D: usize, const I: usize> Mul<&Self> for Vector<D, I>
{
    type Output = Scalar;
    fn mul(self, vector: &Self) -> Self::Output
    {
        self.0 * &vector.0
    }
}

impl<const D: usize, const I: usize> Mul<Vector<D, I>> for &Vector<D, I>
{
    type Output = Scalar;
    fn mul(self, vector: Vector<D, I>) -> Self::Output
    {
        &self.0 * vector.0
    }
}

impl<const D: usize, const I: usize> Mul for &Vector<D, I>
{
    type Output = Scalar;
    fn mul(self, vector: Self) -> Self::Output
    {
        &self.0 * &vector.0
    }
}