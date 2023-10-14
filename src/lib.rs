#[cfg(feature = "32")]
type Float = f32;

#[cfg(feature = "64")]
type Float = f64;

#[cfg(feature = "math")]
mod math;
