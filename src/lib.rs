#![doc = include_str!("../README.md")]

#[cfg(test)]
mod test;

#[cfg(feature = "constitutive")]
pub mod constitutive;

#[cfg(feature = "fem")]
pub mod fem;

#[cfg(feature = "math")]
pub mod math;

#[cfg(feature = "mechanics")]
pub mod mechanics;

/// Absolute tolerance.
pub const ABS_TOL: f64 = 1e-12;

/// Relative tolerance.
pub const REL_TOL: f64 = 1e-12;

#[cfg(test)]
/// A perturbation.
pub const EPSILON: f64 = 1e-6;

#[allow(dead_code)]
fn get_error_message<'a>() -> &'a str {
    let now = format!("{:?}", std::time::SystemTime::now());
    let length = now.as_bytes().len();
    let number = std::str::from_utf8(&now.as_bytes()[length - 3.. length - 2]).unwrap().parse::<u8>().unwrap();
    match number {
        0 => "I am Error.",
        1..3 => "Game over.",
        3..5 => "Oh dear, you are dead!",
        5 => "Press F to pay respects.",
        6..7 => "Surprise! You're dead!",
        7 => "The cake is a lie.",
        8 => "What a horrible night to have a curse.",
        9.. => "You have died of dysentery."
    }
}
