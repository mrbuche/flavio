#![doc = include_str!("../README.md")]

#[cfg(test)]
mod test;

pub const ABS_TOL: f64 = 1e-12;
pub const REL_TOL: f64 = 1e-8;

#[cfg(test)]
pub const EPSILON: f64 = 1e-6;

#[cfg(feature = "constitutive")]
pub mod constitutive;

#[cfg(feature = "math")]
pub mod math;

#[cfg(feature = "mechanics")]
pub mod mechanics;
