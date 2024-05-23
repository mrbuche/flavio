#[cfg(test)]
mod test;

use super::*;

const G: usize = 4;
const M: usize = 3;
const N: usize = 10;
const O: usize = 10;
const P: usize = 12;
const Q: usize = 4;

const INTEGRATION_WEIGHT: Scalar = 1.0/24.0;

pub struct Tetrahedron<C>
{
    constitutive_models: [C; G],
    integration_weights: Scalars<G>,
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
            integration_weights: Self::calculate_reference_jacobians(&reference_nodal_coordinates) * INTEGRATION_WEIGHT,
            projected_gradient_vectors: Self::calculate_projected_gradient_vectors(&reference_nodal_coordinates)
        }
    }
}

impl<'a, C> CompositeElement<'a, C, G, M, N, O, P, Q> for Tetrahedron<C>
where
    C: Constitutive<'a>
{
    fn calculate_inverse_normalized_projection_matrix() -> NormalizedProjectionMatrix<Q>
    {
        let diag: Scalar = 4.0/640.0;
        let off: Scalar = -1.0/640.0;
        NormalizedProjectionMatrix::new([
            [diag,  off,  off,  off],
            [ off, diag,  off,  off],
            [ off,  off, diag,  off],
            [ off,  off,  off, diag]
        ])
    }
    fn calculate_projected_gradient_vectors(reference_nodal_coordinates: &ReferenceNodalCoordinates<O>) -> ProjectedGradientVectors<G, N>
    {
        let parametric_gradient_operators =
        Self::calculate_standard_gradient_operators().iter()
        .map(|standard_gradient_operator|
            reference_nodal_coordinates * standard_gradient_operator
        ).collect::<ParametricGradientOperators<P>>();
        let reference_jacobians_subelements = Self::calculate_reference_jacobians_subelements(reference_nodal_coordinates);
        let inverse_projection_matrix = Self::calculate_inverse_projection_matrix(&reference_jacobians_subelements);
        Self::calculate_shape_functions_at_integration_points().iter()
        .map(|shape_functions_at_integration_point|
            Self::calculate_standard_gradient_operators_transposed().iter()
            .map(|standard_gradient_operators_a|
                Self::calculate_shape_function_integrals().iter()
                .zip(standard_gradient_operators_a.iter()
                .zip(parametric_gradient_operators.iter()
                .zip(reference_jacobians_subelements.iter())))
                .map(|(shape_function_integral, (standard_gradient_operator, (parametric_gradient_operator, reference_jacobian_subelement)))|
                    (parametric_gradient_operator.inverse_transpose() * standard_gradient_operator) * reference_jacobian_subelement
                    * (shape_functions_at_integration_point * (&inverse_projection_matrix * shape_function_integral))
                ).sum()
            ).collect()
        ).collect()
    }
    fn calculate_reference_jacobians_subelements(reference_nodal_coordinates: &ReferenceNodalCoordinates<N>) -> Scalars<P>
    {
        Self::calculate_standard_gradient_operators().iter()
        .map(|standard_gradient_operator|
            reference_nodal_coordinates * standard_gradient_operator
        ).collect::<ParametricGradientOperators<P>>().iter()
        .map(|parametric_gradient_operator|
            parametric_gradient_operator.determinant()
        ).collect()
    }
    fn calculate_shape_function_integrals() -> ShapeFunctionIntegrals<P, Q>
    {
        ShapeFunctionIntegrals::new([
            [200.0,  40.0,  40.0,  40.0],
            [ 40.0, 200.0,  40.0,  40.0],
            [ 40.0,  40.0, 200.0,  40.0],
            [ 40.0,  40.0,  40.0, 200.0],
            [ 30.0,  70.0,  30.0,  30.0],
            [ 10.0,  50.0,  50.0,  50.0],
            [ 30.0,  30.0,  30.0,  70.0],
            [ 50.0,  50.0,  10.0,  50.0],
            [ 50.0,  50.0,  50.0,  10.0],
            [ 30.0,  30.0,  70.0,  30.0],
            [ 50.0,  10.0,  50.0,  50.0],
            [ 70.0,  30.0,  30.0,  30.0]
        ])
    }
    fn calculate_shape_function_integrals_products() -> ShapeFunctionIntegralsProducts<P, Q>
    {
        ShapeFunctionIntegralsProducts::new([[
            [128.0,  24.0,  24.0,  24.0],
            [ 24.0,   8.0,   4.0,   4.0],
            [ 24.0,   4.0,   8.0,   4.0],
            [ 24.0,   4.0,   4.0,   8.0]
        ], [
            [  8.0,  24.0,   4.0,   4.0],
            [ 24.0, 128.0,  24.0,  24.0],
            [  4.0,  24.0,   8.0,   4.0],
            [  4.0,  24.0,   4.0,   8.0]
        ], [
            [  8.0,   4.0,  24.0,   4.0],
            [  4.0,   8.0,  24.0,   4.0],
            [ 24.0,  24.0, 128.0,  24.0],
            [  4.0,   4.0,  24.0,   8.0]
        ], [
            [  8.0,   4.0,   4.0,  24.0],
            [  4.0,   8.0,   4.0,  24.0],
            [  4.0,   4.0,   8.0,  24.0],
            [ 24.0,  24.0,  24.0, 128.0]
        ], [
            [  7.0,  13.0,   5.0,   5.0],
            [ 13.0,  31.0,  13.0,  13.0],
            [  5.0,  13.0,   7.0,   5.0],
            [  5.0,  13.0,   5.0,   7.0]
        ], [
            [  1.0,   3.0,   3.0,   3.0],
            [  3.0,  17.0,  15.0,  15.0],
            [  3.0,  15.0,  17.0,  15.0],
            [  3.0,  15.0,  15.0,  17.0]
        ], [
            [  7.0,   5.0,   5.0,  13.0],
            [  5.0,   7.0,   5.0,  13.0],
            [  5.0,   5.0,   7.0,  13.0],
            [ 13.0,  13.0,  13.0,  31.0]
        ], [
            [ 17.0,  15.0,   3.0,  15.0],
            [ 15.0,  17.0,   3.0,  15.0],
            [  3.0,   3.0,   1.0,   3.0],
            [ 15.0,  15.0,   3.0,  17.0]
        ], [
            [ 17.0,  15.0,  15.0,   3.0],
            [ 15.0,  17.0,  15.0,   3.0],
            [ 15.0,  15.0,  17.0,   3.0],
            [  3.0,   3.0,   3.0,   1.0]
        ], [
            [  7.0,   5.0,  13.0,   5.0],
            [  5.0,   7.0,  13.0,   5.0],
            [ 13.0,  13.0,  31.0,  13.0],
            [  5.0,   5.0,  13.0,   7.0]
        ], [
            [ 17.0,   3.0,  15.0,  15.0],
            [  3.0,   1.0,   3.0,   3.0],
            [ 15.0,   3.0,  17.0,  15.0],
            [ 15.0,   3.0,  15.0,  17.0]
        ], [
            [ 31.0,  13.0,  13.0,  13.0],
            [ 13.0,   7.0,   5.0,   5.0],
            [ 13.0,   5.0,   7.0,   5.0],
            [ 13.0,   5.0,   5.0,   7.0]
        ]])
    }
    fn calculate_shape_functions_at_integration_points() -> ShapeFunctionsAtIntegrationPoints<G, Q>
    {
        let diag: Scalar = 0.585_410_196_624_968_5;
        let off: Scalar = 0.138_196_601_125_010_5;
        ShapeFunctionsAtIntegrationPoints::new([
            [diag,  off,  off,  off],
            [ off, diag,  off,  off],
            [ off,  off, diag,  off],
            [ off,  off,  off, diag]
        ])
    }
    fn calculate_standard_gradient_operators() -> StandardGradientOperators<M, O, P>
    {
        let mut standard_gradient_operators = StandardGradientOperators::zero();
        let standard_gradient_operators_transposed = Self::calculate_standard_gradient_operators_transposed();
        standard_gradient_operators.iter_mut().enumerate()
        .for_each(|(e, standard_gradient_operator_e)|
            standard_gradient_operator_e.iter_mut()
            .zip(standard_gradient_operators_transposed.iter())
            .for_each(|(standard_gradient_operator_e_n, standard_gradient_operator_transposed_n)|
                standard_gradient_operator_e_n.iter_mut()
                .zip(standard_gradient_operator_transposed_n[e].iter())
                .for_each(|(standard_gradient_operator_e_n_i, standard_gradient_operator_transposed_n_e_i)|
                    *standard_gradient_operator_e_n_i = *standard_gradient_operator_transposed_n_e_i
                )
            )
        );
        standard_gradient_operators
    }
    fn calculate_standard_gradient_operators_transposed() -> StandardGradientOperatorsTransposed<M, O, P>
    {
        let n23: Scalar = 2.0/3.0;
        let n43: Scalar = 4.0/3.0;
        StandardGradientOperatorsTransposed::new([[
            [-2.0, -2.0, -2.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0]
        ], [
            [ 0.0,  0.0,  0.0],
            [ 2.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0]
        ], [
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  2.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0]
        ], [
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  2.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0]
        ], [
            [ 2.0,  0.0,  0.0],
            [-2.0, -2.0, -2.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [-n23, -2.0, -2.0],
            [-n23, -n23, -n23],
            [ 0.0,  0.0, -n23],
            [ 0.0, -n43, -2.0],
            [ 0.0, -2.0, -n43],
            [ 0.0, -n23,  0.0],
            [ n23,  0.0,  0.0],
            [ n23, -n43, -n43]
        ], [
            [ 0.0,  0.0,  0.0],
            [ 0.0,  2.0,  0.0],
            [ 2.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ n43,  2.0,  0.0],
            [ n43,  n43, -n23],
            [ 0.0,  0.0, -n23],
            [ 0.0,  n23,  0.0],
            [ 2.0,  2.0,  n23],
            [ 2.0,  n43,  0.0],
            [ n23,  0.0,  0.0],
            [ n23,  n23,  n23]
        ], [
            [ 0.0,  2.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [-2.0, -2.0, -2.0],
            [ 0.0,  0.0,  0.0],
            [-n23,  0.0,  0.0],
            [-n23, -n23, -n23],
            [ 0.0,  0.0, -n23],
            [ 0.0,  n23,  0.0],
            [-2.0,  0.0, -n43],
            [-2.0, -n23, -2.0],
            [-n43,  0.0, -2.0],
            [-n43,  n23, -n43]
        ], [
            [ 0.0,  0.0,  2.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [-2.0, -2.0, -2.0],
            [-n23,  0.0,  0.0],
            [-n23, -n23, -n23],
            [-2.0, -2.0, -n23],
            [-2.0, -n43,  0.0],
            [ 0.0,  0.0,  n23],
            [ 0.0, -n23,  0.0],
            [-n43, -2.0,  0.0],
            [-n43, -n43,  n23]
        ], [
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  2.0],
            [ 0.0,  0.0,  0.0],
            [ 2.0,  0.0,  0.0],
            [ n43,  0.0,  2.0],
            [ n43, -n23,  n43],
            [ 2.0,  0.0,  n43],
            [ 2.0,  n23,  2.0],
            [ 0.0,  0.0,  n23],
            [ 0.0, -n23,  0.0],
            [ n23,  0.0,  0.0],
            [ n23,  n23,  n23]
        ], [
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  2.0],
            [ 0.0,  2.0,  0.0],
            [-n23,  0.0,  0.0],
            [-n23,  n43,  n43],
            [ 0.0,  2.0,  n43],
            [ 0.0,  n23,  0.0],
            [ 0.0,  0.0,  n23],
            [ 0.0,  n43,  2.0],
            [ n23,  2.0,  2.0],
            [ n23,  n23,  n23]
        ]])
    }
    fn get_constitutive_models(&self) -> &[C; G]
    {
        &self.constitutive_models
    }
    fn get_integration_weights(&self) -> &Scalars<G>
    {
        &self.integration_weights
    }
    fn get_projected_gradient_vectors(&self) -> &ProjectedGradientVectors<G, N>
    {
        &self.projected_gradient_vectors
    }
}

impl<'a, C> ElasticFiniteElement<'a, C, G, N> for Tetrahedron<C>
where
    C: Elastic<'a>
{
    fn calculate_nodal_forces(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalForces<N>
    {
        self.calculate_nodal_forces_composite_element(nodal_coordinates)
    }
    fn calculate_nodal_stiffnesses(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalStiffnesses<N>
    {
        self.calculate_nodal_stiffnesses_composite_element(nodal_coordinates)
    }
}

impl<'a, C> ElasticCompositeElement<'a, C, G, M, N, O, P, Q> for Tetrahedron<C>
where
    C: Elastic<'a>
{}

impl<'a, C> HyperelasticFiniteElement<'a, C, G, N> for Tetrahedron<C>
where
    C: Hyperelastic<'a>
{
    fn calculate_helmholtz_free_energy(&self, nodal_coordinates: &NodalCoordinates<N>) -> Scalar
    {
        self.calculate_helmholtz_free_energy_composite_element(nodal_coordinates)
    }
}

impl<'a, C> HyperelasticCompositeElement<'a, C, G, M, N, O, P, Q> for Tetrahedron<C>
where
    C: Hyperelastic<'a>
{}

impl<'a, C> ViscoelasticFiniteElement<'a, C, G, N> for Tetrahedron<C>
where
    C: Viscoelastic<'a>
{
    fn calculate_nodal_forces(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> NodalForces<N>
    {
        self.calculate_nodal_forces_composite_element(nodal_coordinates, nodal_velocities)
    }
    fn calculate_nodal_stiffnesses(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> NodalStiffnesses<N>
    {
        self.calculate_nodal_stiffnesses_composite_element(nodal_coordinates, nodal_velocities)
    }
}

impl<'a, C> ViscoelasticCompositeElement<'a, C, G, M, N, O, P, Q> for Tetrahedron<C>
where
    C: Viscoelastic<'a>
{}

impl<'a, C> ElasticHyperviscousFiniteElement<'a, C, G, N> for Tetrahedron<C>
where
    C: ElasticHyperviscous<'a>
{
    fn calculate_viscous_dissipation(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> Scalar
    {
        self.calculate_viscous_dissipation_composite_element(nodal_coordinates, nodal_velocities)
    }
    fn calculate_dissipation_potential(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> Scalar
    {
        self.calculate_dissipation_potential_composite_element(nodal_coordinates, nodal_velocities)
    }
}

impl<'a, C> ElasticHyperviscousCompositeElement<'a, C, G, M, N, O, P, Q> for Tetrahedron<C>
where
    C: ElasticHyperviscous<'a>
{}

impl<'a, C> HyperviscoelasticFiniteElement<'a, C, G, N> for Tetrahedron<C>
where
    C: Hyperviscoelastic<'a>
{
    fn calculate_helmholtz_free_energy(&self, nodal_coordinates: &NodalCoordinates<N>) -> Scalar
    {
        self.calculate_helmholtz_free_energy_composite_element(nodal_coordinates)
    }
}

impl<'a, C> HyperviscoelasticCompositeElement<'a, C, G, M, N, O, P, Q> for Tetrahedron<C>
where
    C: Hyperviscoelastic<'a>
{}