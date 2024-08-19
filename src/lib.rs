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
    let number = &now[length - 3..length - 2].parse::<u8>().unwrap();
    match number {
        0 => "Game over.",
        1 => "I am Error.",
        2 => "Oh dear, you are dead!",
        3 => "Press F to pay respects.",
        4 => "Surprise! You're dead!",
        5 => "The cake is a lie.",
        6 => "This is not your grave, but you are welcome in it.",
        7 => "What a horrible night to have a curse.",
        8 => "You have died of dysentery.",
        9.. => "You've met with a terrible fate, haven't you?",
    }
}
