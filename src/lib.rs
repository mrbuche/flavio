#![doc = include_str!("../README.md")]

#[cfg(test)]
mod test;

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

#[cfg(test)]
#[cfg(feature = "32")]
pub const EPSILON: Float = 1e-2;

#[cfg(test)]
#[cfg(feature = "64")]
pub const EPSILON: Float = 1e-6;

#[cfg(feature = "constitutive")]
pub mod constitutive;

#[cfg(feature = "math")]
pub mod math;

#[cfg(feature = "mechanics")]
pub mod mechanics;
