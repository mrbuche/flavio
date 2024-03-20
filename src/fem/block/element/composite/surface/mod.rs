#[cfg(test)]
mod test;

pub mod triangle;

use super::*;

pub trait CompositeSurfaceElement<'a, C, const G: usize, const M: usize, const N: usize, const O: usize, const P: usize, const Q: usize>
where
    C: Constitutive<'a>,
    Self: CompositeElement<'a, C, G, M, N, O, P, Q>
{
    fn calculate_normals<const I: usize>(nodal_coordinates: &Coordinates<I, O>) -> Normals<I, P>
    {
        todo!("Need to calculate and store reference normals for each subtriangle during initialization.")
    }
    fn get_reference_normals(&self) -> &ReferenceNormals<P>;
}

macro_rules! composite_surface_element_boilerplate
{
    ($element: ident) =>
    {
        impl<'a, C> CompositeSurfaceElement<'a, C, G, M, N, O, P, Q> for $element<C>
        where
            C: Constitutive<'a>
        {
            fn get_reference_normals(&self) -> &ReferenceNormals<P>
            {
                &self.reference_normals
            }
        }
    }
}
pub(crate) use composite_surface_element_boilerplate;