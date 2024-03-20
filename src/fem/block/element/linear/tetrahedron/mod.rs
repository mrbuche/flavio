#[cfg(test)]
mod test;

use super::*;

const G: usize = 1;
const M: usize = 3;
const N: usize = 4;
const O: usize = 4;

const INTEGRATION_WEIGHT: Scalar = 1.0/6.0;

pub struct Tetrahedron<C>
{
    constitutive_model: C,
    gradient_vectors: GradientVectors<N>
}

impl<'a, C> FiniteElement<'a, C, G, N> for Tetrahedron<C>
where
    C: Constitutive<'a>
{
    fn new(constitutive_model_parameters: Parameters<'a>, reference_nodal_coordinates: ReferenceNodalCoordinates<N>) -> Self
    {
        Self
        {
            constitutive_model: <C>::new(constitutive_model_parameters),
            gradient_vectors: Self::calculate_gradient_vectors(&reference_nodal_coordinates)
        }
    }
}

impl<'a, C> LinearElement<'a, C, G, M, N, O> for Tetrahedron<C>
where
    C: Constitutive<'a>
{
    fn calculate_gradient_vectors(reference_nodal_coordinates: &ReferenceNodalCoordinates<N>) -> GradientVectors<N>
    {
        let standard_gradient_operator = Self::calculate_standard_gradient_operator();
        (reference_nodal_coordinates * &standard_gradient_operator).inverse_transpose() * standard_gradient_operator
    }
    fn calculate_standard_gradient_operator() -> StandardGradientOperator<M, O>
    {
        StandardGradientOperator::new([
            [-1.0, -1.0, -1.0],
            [ 1.0,  0.0,  0.0],
            [ 0.0,  1.0,  0.0],
            [ 0.0,  0.0,  1.0]
        ])
    }
    fn get_constitutive_model(&self) -> &C
    {
        &self.constitutive_model
    }
    fn get_gradient_vectors(&self) -> &GradientVectors<N>
    {
        &self.gradient_vectors
    }
    fn get_integration_weight(&self) -> &Scalar
    {
        &INTEGRATION_WEIGHT
    }
}

super::linear_element_boilerplate!(Tetrahedron);