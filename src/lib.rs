#[cfg(test)]
mod test;

#[cfg(feature = "math")]
pub mod math;

#[cfg(feature = "32")]
type Float = f32;

#[cfg(feature = "64")]
type Float = f64;

#[cfg(feature = "32")]
pub const ABS_TOL: Float = 1e-4;

#[cfg(feature = "64")]
pub const ABS_TOL: Float = 1e-13;

#[cfg(feature = "32")]
pub const REL_TOL: Float = 1e-4;

#[cfg(feature = "64")]
pub const REL_TOL: Float = 1e-13;

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
// and also
// impl TensorRank2Traits<D> for Tensor<D, I, J>