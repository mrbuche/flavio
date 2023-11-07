#![doc = include_str!("../README.md")]

#[cfg(test)]
mod test;

/// Absolute tolerance.
pub const ABS_TOL: f64 = 1e-12;

/// Relative tolerance.
pub const REL_TOL: f64 = 1e-12;

#[cfg(test)]
/// A perturbation.
pub const EPSILON: f64 = 1e-6;

#[cfg(feature = "constitutive")]
/// Constitutive model library.
pub mod constitutive;

#[cfg(feature = "fem")]
/// Finite element library.
pub mod fem;

#[cfg(feature = "math")]
/// Mathematics library.
pub mod math;

#[cfg(feature = "mechanics")]
/// Mechanics library.
pub mod mechanics;
