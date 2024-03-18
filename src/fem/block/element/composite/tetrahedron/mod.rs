#[cfg(test)]
mod test;

use super::*;

const G: usize = 4;
const M: usize = 3;
const N: usize = 9;
const O: usize = 9;

pub struct Tetrahedron<C>
{
    constitutive_models: [C; G],
    integration_weights: IntegrationWeights<G>,
    projected_gradient_vectors: ProjectedGradientVectors<G, N>
}

impl<'a, C> FiniteElement<'a, C, G, N> for Tetrahedron<C>
where
    C: Constitutive<'a>
{
    fn new(constitutive_model_parameters: Parameters<'a>, reference_nodal_coordinates: ReferenceNodalCoordinates<N>) -> Self
    {
        Self
        {
            constitutive_models: std::array::from_fn(|_| <C>::new(constitutive_model_parameters)),
            integration_weights: IntegrationWeights::new([1.0/24.0; G]),
            projected_gradient_vectors: Self::calculate_projected_gradient_vectors(&reference_nodal_coordinates)
        }
    }
}

impl<'a, C> CompositeElement<'a, C, G, M, N, O> for Tetrahedron<C>
where
    C: Constitutive<'a>
{
    fn calculate_projected_gradient_vectors(reference_nodal_coordinates: &ReferenceNodalCoordinates<N>) -> ProjectedGradientVectors<G, N>
    {
        todo!()
    }
    fn get_constitutive_models(&self) -> &[C; G]
    {
        &self.constitutive_models
    }
    fn get_integration_weights(&self) -> &IntegrationWeights<G>
    {
        &self.integration_weights
    }
    fn get_projected_gradient_vectors(&self) -> &ProjectedGradientVectors<G, N>
    {
        &self.projected_gradient_vectors
    }
}

super::composite_element_boilerplate!(Tetrahedron);