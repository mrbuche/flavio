#[cfg(test)]
mod test;

pub mod triangle;

use super::*;

pub trait CompositeSurfaceElement<'a, C, const G: usize, const M: usize, const N: usize, const O: usize, const P: usize, const Q: usize>
where
    C: Constitutive<'a>,
    Self: CompositeElement<'a, C, G, M, N, O, P, Q>
{}