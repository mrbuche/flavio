#[cfg(feature = "32")]
type Float = f32;

#[cfg(feature = "64")]
type Float = f64;

#[cfg(feature = "math")]
pub mod math;

// implement tensors as Tensor<D, I, J>
// similarly for vectors, etc.
// using 0 for reference
// 1 for current
// 2 for intermediate (like F_2)
// then do things like:
// impl Mul<Tensor<D, J, K> for Tensor<D, I, J>
// type Output = Tensor<D, I, K>
// fn inverse(&self) -> Tensor<D, J, I> // for Tensor<D, I, J>
// fn from_dyad(vector_a: Vector<D, I>, vector_b: Vector<D, J>) -> Tensor<D, I, J>
