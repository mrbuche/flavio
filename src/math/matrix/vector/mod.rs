use crate::math::{Tensor, TensorRank0, TensorVec, write_tensor_rank_0};
use std::{
    fmt::{Display, Formatter, Result},
    ops::{Add, AddAssign, Div, DivAssign, Index, IndexMut, Mul, MulAssign, Sub, SubAssign},
};

/// A vector.
#[derive(Debug)]
pub struct Vector(Vec<TensorRank0>);

impl Display for Vector {
    fn fmt(&self, f: &mut Formatter) -> Result {
        write!(f, "\x1B[s")?;
        write!(f, "[")?;
        self.0.chunks(5).enumerate().try_for_each(|(i, chunk)| {
            chunk
                .iter()
                .try_for_each(|entry| write_tensor_rank_0(f, entry))?;
            if (i + 1) * 5 < self.len() {
                writeln!(f, "\x1B[2D,")?;
                write!(f, "\x1B[u")?;
                write!(f, "\x1B[{}B ", i + 1)?;
            }
            Ok(())
        })?;
        write!(f, "\x1B[2D]")?;
        Ok(())
    }
}

impl PartialEq for Vector {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl FromIterator<TensorRank0> for Vector {
    fn from_iter<Ii: IntoIterator<Item = TensorRank0>>(into_iterator: Ii) -> Self {
        Self(Vec::from_iter(into_iterator))
    }
}

impl Index<usize> for Vector {
    type Output = TensorRank0;
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl IndexMut<usize> for Vector {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}

impl Tensor for Vector {
    type Item = TensorRank0;
    fn copy(&self) -> Self {
        self.iter().map(|entry| entry.copy()).collect()
    }
    fn iter(&self) -> impl Iterator<Item = &Self::Item> {
        self.0.iter()
    }
    fn iter_mut(&mut self) -> impl Iterator<Item = &mut Self::Item> {
        self.0.iter_mut()
    }
}

impl<'a> TensorVec<'a> for Vector {
    type Item = TensorRank0;
    type Slice = &'a [TensorRank0];
    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
    fn len(&self) -> usize {
        self.0.len()
    }
    fn new(slice: Self::Slice) -> Self {
        slice.iter().copied().collect()
    }
    fn zero(len: usize) -> Self {
        Self(vec![0.0; len])
    }
}
