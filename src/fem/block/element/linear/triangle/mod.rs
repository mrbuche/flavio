#[cfg(test)]
mod test;

use super::*;

const G: usize = 1;
const M: usize = 2;
const N: usize = 3;

pub struct Triangle<'a, C>
{
    constitutive_models: [C; G],
    gradient_vectors: GradientVectors<N>,
    reference_normal: ReferenceNormal,
    thickness: &'a Scalar
}

impl<'a, C> FiniteElement<'a, C, G, N> for Triangle<'a, C>
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
            reference_normal: Self::calculate_reference_normal(&Self::calculate_reference_dual_basis_vectors(&reference_nodal_coordinates)),
            thickness: &1.23
            // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            // use a set_thickness() thing?
        }
    }
}

impl<'a, C> LinearSurfaceElement<'a, C, G, M, N> for Triangle<'a, C>
where
    C: Constitutive<'a>
{
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
    fn get_reference_normal(&self) -> &ReferenceNormal
    {
        &self.reference_normal
    }
    fn get_thickness(&self) -> &Scalar
    {
        self.thickness
    }
}

impl<'a, C> ElasticLinearSurfaceElement<'a, C, G, M, N> for Triangle<'a, C>
where
    C: Elastic<'a>
{}