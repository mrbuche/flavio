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
    scaled_composite_jacobians: Scalars<G>,
    scaled_reference_normals: ScaledReferenceNormals<G, P>
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
            scaled_composite_jacobians: Self::calculate_scaled_composite_jacobian_at_integration_points(&reference_nodal_coordinates),
            scaled_reference_normals: Self::calculate_scaled_reference_normals(&reference_nodal_coordinates)
        }
    }
}

impl<'a, C> CompositeElement<'a, C, G, M, N, O, P, Q> for Triangle<C>
where
    C: Constitutive<'a>
{
    fn calculate_deformation_gradients(&self, nodal_coordinates: &NodalCoordinates<N>) -> DeformationGradients<G>
    {
        self.calculate_deformation_gradients_composite_surface_element(nodal_coordinates)
    }
    fn calculate_deformation_gradient_rates(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> DeformationGradientRates<G>
    {
        self.calculate_deformation_gradient_rates_composite_surface_element(nodal_coordinates, nodal_velocities)
    }
    fn calculate_projected_gradient_vectors(reference_nodal_coordinates: &ReferenceNodalCoordinates<N>) -> ProjectedGradientVectors<G, N>
    {
        let jacobians = Self::calculate_jacobians(reference_nodal_coordinates);
        let dual_bases = Self::calculate_dual_bases(reference_nodal_coordinates);
        let inverse_projection_matrix = Self::calculate_inverse_projection_matrix(&jacobians);
        Self::calculate_shape_functions_at_integration_points().iter()
        .map(|shape_functions_at_integration_point|
            Self::calculate_standard_gradient_operators_transposed().iter()
            .map(|standard_gradient_operators_a|
                Self::calculate_shape_function_integrals().iter()
                .zip(standard_gradient_operators_a.iter()
                .zip(dual_bases.iter()
                .zip(jacobians.iter())))
                .map(|(shape_function_integral, (standard_gradient_operator, (dual_basis_vectors, jacobian)))|
                    dual_basis_vectors.iter()
                    .zip(standard_gradient_operator.iter())
                    .map(|(dual_basis_vector, standard_gradient_operator_mu)|
                        dual_basis_vector * standard_gradient_operator_mu
                    ).sum::<Vector<0>>() * jacobian * (
                        shape_functions_at_integration_point * (&inverse_projection_matrix * shape_function_integral)
                    )
                ).sum()
            ).collect()
        ).collect()
    }
    composite_surface_element_boilerplate_inner!{}
}

composite_surface_element_boilerplate!(Triangle);