#[cfg(test)]
mod test;

use super::*;

const G: usize = 1;
const M: usize = 2;
const N: usize = 3;
const O: usize = 3;

const INTEGRATION_WEIGHT: Scalar = 1.0/2.0;

pub struct Triangle<C>
{
    constitutive_model: C,
    gradient_vectors: GradientVectors<N>,
    reference_normal: ReferenceNormal
}

impl<'a, C> FiniteElement<'a, C, G, N> for Triangle<C>
where
    C: Constitutive<'a>
{
    fn new(constitutive_model_parameters: Parameters<'a>, reference_nodal_coordinates: ReferenceNodalCoordinates<N>) -> Self
    {
        Self
        {
            constitutive_model: <C>::new(constitutive_model_parameters),
            gradient_vectors: Self::calculate_gradient_vectors(&reference_nodal_coordinates),
            reference_normal: Self::calculate_reference_normal(&reference_nodal_coordinates)
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
    linear_surface_element_boilerplate_inner!{}
}

super::linear_surface_element_boilerplate!(Triangle);