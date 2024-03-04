#[cfg(test)]
mod test;

use super::*;

const G: usize = 1;
const M: usize = 2;
const N: usize = 3;
const O: usize = 3;

pub struct Triangle<C>
{
    constitutive_models: [C; G],
    gradient_vectors: GradientVectors<N>,
    reference_normal: ReferenceNormal
}

impl<'a, C> FiniteElement<'a, C, G, N> for Triangle<C>
where
    C: Constitutive<'a>
{
    fn get_constitutive_models(&self) -> &[C; G]
    {
        &self.constitutive_models
    }
    fn get_integration_weights(&self) -> IntegrationWeights<G>
    {
        IntegrationWeights::new([1.0; G])
    }
    fn new(constitutive_model_parameters: Parameters<'a>, reference_nodal_coordinates: ReferenceNodalCoordinates<N>) -> Self
    {
        Self
        {
            constitutive_models: std::array::from_fn(|_| <C>::new(constitutive_model_parameters)),
            gradient_vectors: Self::calculate_gradient_vectors(&reference_nodal_coordinates),
            reference_normal: Self::calculate_normal(&reference_nodal_coordinates)
        }
    }
}

impl<'a, C> LinearElement<'a, C, G, M, N, O> for Triangle<C>
where
    C: Constitutive<'a>
{
    fn calculate_deformation_gradient(&self, nodal_coordinates: &NodalCoordinates<N>) -> DeformationGradient
    {
        self.calculate_deformation_gradient_linear_surface_element(nodal_coordinates)
    }
    fn calculate_deformation_gradient_rate(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> DeformationGradientRate
    {
        self.calculate_deformation_gradient_rate_linear_surface_element(nodal_coordinates, nodal_velocities)
    }
    fn calculate_gradient_vectors(reference_nodal_coordinates: &ReferenceNodalCoordinates<N>) -> GradientVectors<N>
    {
        Self::calculate_gradient_vectors_linear_surface_element(reference_nodal_coordinates)
    }
    fn calculate_standard_gradient_operator() -> StandardGradientOperator<M, N>
    {
        StandardGradientOperator::new([
            [-1.0, -1.0],
            [ 1.0,  0.0],
            [ 0.0,  1.0]
        ])
    }
    fn get_gradient_vectors(&self) -> &GradientVectors<N>
    {
        &self.gradient_vectors
    }
}

impl<'a, C> ElasticFiniteElement<'a, C, G, N> for Triangle<C>
where
    C: Elastic<'a>
{
    fn calculate_nodal_forces(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalForces<N>
    {
        self.calculate_nodal_forces_linear_surface_element(nodal_coordinates)
    }
    fn calculate_nodal_stiffnesses(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalStiffnesses<N>
    {
        self.calculate_nodal_stiffnesses_linear_surface_element(nodal_coordinates)
    }
}

super::linear_surface_element_boilerplate!(Triangle);