use crate::math::{Tensor, TensorRank0, TensorVec};

/// ???
pub struct Vector(Vec<TensorRank0>);

impl FromIterator<TensorRank0> for Vector {
    fn from_iter<Ii: IntoIterator<Item = TensorRank0>>(into_iterator: Ii) -> Self {
        Self(Vec::from_iter(into_iterator))
    }
}

// impl Tensor for Vector {
//     type Item = TensorRank0;
//     fn copy(&self) -> Self {
//         self.iter().map(|entry| entry.copy()).collect()
//     }
//     fn get_at(&self, indices: &[usize]) -> &TensorRank0 {
//         &self[indices[0]]
//     }
//     fn get_at_mut(&mut self, indices: &[usize]) -> &mut TensorRank0 {
//         &mut self[indices[0]]
//     }
//     fn iter(&self) -> impl Iterator<Item = &Self::Item> {
//         self.0.iter()
//     }
//     fn iter_mut(&mut self) -> impl Iterator<Item = &mut Self::Item> {
//         self.0.iter_mut()
//     }
// }

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
