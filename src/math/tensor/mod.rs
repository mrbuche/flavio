pub mod rank_0;
pub mod rank_1;
pub mod rank_2;
pub mod rank_3;
pub mod rank_4;

/// A value-to-value conversion that does not consume the input value.
///
/// This is as opposed to [`Into`](https://doc.rust-lang.org/std/convert/trait.Into.html), which consumes the input value.
pub trait Convert<T> {
    /// Converts this type into the (usually inferred) input type.
    fn convert(&self) -> T;
}
