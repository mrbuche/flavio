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
    fn calculate_inverse_normalized_projection_matrix() -> NormalizedProjectionMatrix<Q>
    {
        let diag: Scalar = 3.0/64.0;
        let off: Scalar = -1.0/64.0;
        NormalizedProjectionMatrix::new([
            [diag,  off,  off],
            [ off, diag,  off],
            [ off,  off, diag]
        ])
    }
    fn calculate_jacobians(reference_nodal_coordinates: &ReferenceNodalCoordinates<O>) -> Scalars<P>
    {
        Self::calculate_bases(reference_nodal_coordinates).iter()
        .map(|basis_vectors|
            basis_vectors[0].cross(&basis_vectors[1]).norm()
        ).collect()
    }
    fn calculate_projected_gradient_vectors(reference_nodal_coordinates: &ReferenceNodalCoordinates<N>) -> ProjectedGradientVectors<G, N>
    {
        let jacobians = Self::calculate_jacobians(reference_nodal_coordinates);
        let dual_bases = Self::calculate_dual_bases(reference_nodal_coordinates);
        let inverse_projection_matrix =
        Self::calculate_shape_function_integrals_products().iter()
        .zip(jacobians.iter())
        .map(|(shape_function_integrals_products, jacobian)|
            shape_function_integrals_products * jacobian
        ).sum::<ProjectionMatrix<Q>>().inverse();
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
    fn calculate_shape_function_integrals() -> ShapeFunctionIntegrals<P, Q>
    {
        ShapeFunctionIntegrals::new([
            [32.0,  8.0,  8.0],
            [ 8.0, 32.0,  8.0],
            [ 8.0,  8.0, 32.0],
            [16.0, 16.0, 16.0]
        ])
    }
    fn calculate_shape_function_integrals_products() -> ShapeFunctionIntegralsProducts<P, Q>
    {
        ShapeFunctionIntegralsProducts::new([[
            [22.0,  5.0,  5.0],
            [ 5.0,  2.0,  1.0],
            [ 5.0,  1.0,  2.0]
        ], [
            [ 2.0,  5.0,  1.0],
            [ 5.0, 22.0,  5.0],
            [ 1.0,  5.0,  2.0]
        ], [
            [ 2.0,  1.0,  5.0],
            [ 1.0,  2.0,  5.0],
            [ 5.0,  5.0, 22.0]
        ], [
            [ 6.0,  5.0,  5.0],
            [ 5.0,  6.0,  5.0],
            [ 5.0,  5.0,  6.0]
        ]])
    }
    fn calculate_shape_functions_at_integration_points() -> ShapeFunctionsAtIntegrationPoints<G, Q>
    {
        let diag: Scalar = 0.666_666_666_666_666_6;
        let off: Scalar = 0.166_666_666_666_666_7;
        ShapeFunctionsAtIntegrationPoints::new([
            [diag,  off,  off],
            [ off, diag,  off],
            [ off,  off, diag]
        ])
    }
    fn calculate_standard_gradient_operators() -> StandardGradientOperators<M, O, P>
    {
        StandardGradientOperators::new([[
            [ 2.0,  0.0],
            [ 0.0,  0.0],
            [ 0.0,  0.0],
            [ 0.0,  2.0],
            [ 0.0,  0.0],
            [-2.0, -2.0]
        ], [
            [ 0.0,  0.0],
            [ 0.0,  2.0],
            [ 0.0,  0.0],
            [ 2.0,  0.0],
            [-2.0, -2.0],
            [ 0.0,  0.0]
        ], [
            [ 0.0,  0.0],
            [ 0.0,  0.0],
            [-2.0, -2.0],
            [ 0.0,  0.0],
            [ 0.0,  2.0],
            [ 2.0,  0.0]
        ], [
            [ 0.0,  0.0],
            [ 0.0,  0.0],
            [ 0.0,  0.0],
            [ 2.0,  2.0],
            [-2.0,  0.0],
            [ 0.0, -2.0]
        ]])
    }
    fn calculate_standard_gradient_operators_transposed() -> StandardGradientOperatorsTransposed<M, O, P>
    {
        let standard_gradient_operators = Self::calculate_standard_gradient_operators();
        let mut standard_gradient_operators_transposed = StandardGradientOperatorsTransposed::zero();
        standard_gradient_operators_transposed.iter_mut().enumerate()
        .for_each(|(n, standard_gradient_operators_transposed_n)|
            standard_gradient_operators_transposed_n.iter_mut()
            .zip(standard_gradient_operators.iter())
            .for_each(|(standard_gradient_operators_transposed_n_e, standard_gradient_operators_e)|
                standard_gradient_operators_transposed_n_e.iter_mut()
                .zip(standard_gradient_operators_e[n].iter())
                .for_each(|(standard_gradient_operators_transposed_n_e_i, standard_gradient_operators_e_n_i)|
                    *standard_gradient_operators_transposed_n_e_i = *standard_gradient_operators_e_n_i
                )
            )
        );
        standard_gradient_operators_transposed
    }
    fn get_constitutive_models(&self) -> &[C; G]
    {
        &self.constitutive_models
    }
    fn get_integration_weight() -> Scalar
    {
        INTEGRATION_WEIGHT
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