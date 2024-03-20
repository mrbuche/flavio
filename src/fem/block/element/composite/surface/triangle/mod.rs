#[cfg(test)]
mod test;

use super::*;

const G: usize = 3;
const M: usize = 2;
const N: usize = 6;
const O: usize = 6;
const P: usize = 4;
const Q: usize = 3;

const INTEGRATION_WEIGHT: Scalar = 1.0/6.0;

pub struct Triangle<C>
{
    constitutive_models: [C; G],
    projected_gradient_vectors: ProjectedGradientVectors<G, N>,
    reference_normals: ReferenceNormals<P>,
    scaled_composite_jacobians: Scalars<G>
}

impl<'a, C> FiniteElement<'a, C, G, N> for Triangle<C>
where
    C: Constitutive<'a>
{
    fn new(constitutive_model_parameters: Parameters<'a>, reference_nodal_coordinates: ReferenceNodalCoordinates<N>) -> Self
    {
        Self
        {
            constitutive_models: std::array::from_fn(|_| <C>::new(constitutive_model_parameters)),
            projected_gradient_vectors: Self::calculate_projected_gradient_vectors(&reference_nodal_coordinates),
            reference_normals: Self::calculate_normals(&reference_nodal_coordinates),
            scaled_composite_jacobians: Self::calculate_scaled_composite_jacobian_at_integration_points(&reference_nodal_coordinates)
        }
    }
}

impl<'a, C> CompositeElement<'a, C, G, M, N, O, P, Q> for Triangle<C>
where
    C: Constitutive<'a>
{
    fn calculate_inverse_normalized_projection_matrix() -> NormalizedProjectionMatrix<Q>
    {
        todo!("Be careful and see what you can make a default implementation in composite/mod.rs to avoid copying code.")
    }
    fn calculate_jacobians_and_parametric_gradient_operators(reference_nodal_coordinates: &ReferenceNodalCoordinates<N>) -> (Scalars<P>, ParametricGradientOperators<P>)
    {
        todo!("Be careful and see what you can make a default implementation in composite/mod.rs to avoid copying code.")
    }
    fn calculate_projected_gradient_vectors(reference_nodal_coordinates: &ReferenceNodalCoordinates<N>) -> ProjectedGradientVectors<G, N>
    {
        todo!("Be careful and see what you can make a default implementation in composite/mod.rs to avoid copying code.")
    }
    fn calculate_scaled_composite_jacobian_at_integration_points(reference_nodal_coordinates: &ReferenceNodalCoordinates<N>) -> Scalars<G>
    {
         todo!("Be careful and see what you can make a default implementation in composite/mod.rs to avoid copying code.")
    }
    fn calculate_shape_function_integrals() -> ShapeFunctionIntegrals<P, Q>
    {
         todo!("Be careful and see what you can make a default implementation in composite/mod.rs to avoid copying code.")
    }
    fn calculate_shape_function_integrals_products() -> ShapeFunctionIntegralsProducts<P, Q>
    {
         todo!("Be careful and see what you can make a default implementation in composite/mod.rs to avoid copying code.")
    }
    fn calculate_shape_functions_at_integration_points() -> ShapeFunctionsAtIntegrationPoints<G, Q>
    {
         todo!("Be careful and see what you can make a default implementation in composite/mod.rs to avoid copying code.")
    }
    fn calculate_standard_gradient_operators() -> StandardGradientOperators<M, O, P>
    {
         todo!("Be careful and see what you can make a default implementation in composite/mod.rs to avoid copying code.")
    }
    fn calculate_standard_gradient_operators_transposed() -> StandardGradientOperatorsTransposed<M, O, P>
    {
         todo!("Be careful and see what you can make a default implementation in composite/mod.rs to avoid copying code.")
    }
    fn get_constitutive_models(&self) -> &[C; G]
    {
        &self.constitutive_models
    }
    fn get_projected_gradient_vectors(&self) -> &ProjectedGradientVectors<G, N>
    {
        &self.projected_gradient_vectors
    }
    fn get_scaled_composite_jacobians(&self) -> &Scalars<G>
    {
        &self.scaled_composite_jacobians
    }
}

composite_surface_element_boilerplate!(Triangle);