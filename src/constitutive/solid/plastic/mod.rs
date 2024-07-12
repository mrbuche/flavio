//! Plastic constitutive models.

#[cfg(test)]
pub mod test;

use super::*;

/// Required methods for plastic constitutive models.
pub trait Plastic<'a>
where
    Self: Solid<'a>
{
    // assume that Fp is an ISV here too since can't use these by themselves?
    // so then you make getter/setter methods here?
}
