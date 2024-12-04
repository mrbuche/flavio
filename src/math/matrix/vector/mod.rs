use crate::math::{Tensor, TensorRank0, TensorVec};

/// ???
pub struct Vector(Vec<TensorRank0>);

impl FromIterator<TensorRank0> for Vector {
    fn from_iter<Ii: IntoIterator<Item = TensorRank0>>(into_iterator: Ii) -> Self {
        Self(Vec::from_iter(into_iterator))
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
