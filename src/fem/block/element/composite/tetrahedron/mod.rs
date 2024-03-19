#[cfg(test)]
mod test;

use super::*;

const G: usize = 4;
const M: usize = 3;
const N: usize = 10;
const O: usize = 10;
const P: usize = 12;
const Q: usize = 4;

// is it 1/4, or 1/24?
// has implications in other elements
// i.e. linear tetrahedron would have weight 1/6, linear tri and localiztion wedge 1/2
// "The weights of a quadrature rule always sum to the volume of the reference element."
const INTEGRATION_WEIGHT: Scalar = 1.0/24.0;

pub struct Tetrahedron<C>
{
    constitutive_models: [C; G],
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
            projected_gradient_vectors: Self::calculate_projected_gradient_vectors(&reference_nodal_coordinates)
        }
    }
}

impl<'a, C> CompositeElement<'a, C, G, M, N, O, P, Q> for Tetrahedron<C>
where
    C: Constitutive<'a>
{
    fn calculate_projected_gradient_vectors(reference_nodal_coordinates: &ReferenceNodalCoordinates<N>) -> ProjectedGradientVectors<G, N>
    {
        let parametric_gradient_operators =
        Self::calculate_standard_gradient_operators().iter()
        .map(|standard_gradient_operator|
            reference_nodal_coordinates * standard_gradient_operator
        ).collect::<ParametricGradientOperators<P>>();
        let jacobians =
        parametric_gradient_operators.iter()
        .map(|parametric_gradient_operator|
            parametric_gradient_operator.determinant()
        ).collect::<Scalars<P>>();
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
                .zip(parametric_gradient_operators.iter()
                .zip(jacobians.iter())))
                .map(|(shape_function_integral, (standard_gradient_operator, (parametric_gradient_operator, jacobian)))|
                    (parametric_gradient_operator.inverse_transpose() * standard_gradient_operator) * jacobian
                    * (shape_functions_at_integration_point * (&inverse_projection_matrix * shape_function_integral))
                ).sum()
            ).collect()
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
            [ off,  off,  off, diag],
        ])
    }
    fn calculate_standard_gradient_operators() -> StandardGradientOperators<M, O, P>
    {
        let n23: Scalar = 2.0/3.0;
        let n43: Scalar = 4.0/3.0;
        StandardGradientOperators::new([[
            [-2.0, -2.0, -2.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 2.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  2.0,  0.0],
            [ 0.0,  0.0,  2.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0]
        ], [
            [ 0.0,  0.0,  0.0],
            [ 2.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [-2.0, -2.0, -2.0],
            [ 0.0,  2.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  2.0],
            [ 0.0,  0.0,  0.0]
        ], [
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  2.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 2.0,  0.0,  0.0],
            [-2.0, -2.0, -2.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  2.0]
        ], [
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  2.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [-2.0, -2.0, -2.0],
            [ 2.0,  0.0,  0.0],
            [ 0.0,  2.0,  0.0]
        ], [
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [-n23, -2.0, -2.0],
            [ n43,  2.0,  0.0],
            [-n23,  0.0,  0.0],
            [-n23,  0.0,  0.0],
            [ n43,  0.0,  2.0],
            [-n23,  0.0,  0.0]
        ], [
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [-n23, -n23, -n23],
            [ n43,  n43, -n23],
            [-n23, -n23, -n23],
            [-n23, -n23, -n23],
            [ n43, -n23,  n43],
            [-n23,  n43,  n43]
        ], [
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0, -n23],
            [ 0.0,  0.0, -n23],
            [ 0.0,  0.0, -n23],
            [-2.0, -2.0, -n23],
            [ 2.0,  0.0,  n43],
            [ 0.0,  2.0,  n43]
        ], [
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0, -n43, -2.0],
            [ 0.0,  n23,  0.0],
            [ 0.0,  n23,  0.0],
            [-2.0, -n43,  0.0],
            [ 2.0,  n23,  2.0],
            [ 0.0,  n23,  0.0]
        ], [
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0, -2.0, -n43],
            [ 2.0,  2.0,  n23],
            [-2.0,  0.0, -n43],
            [ 0.0,  0.0,  n23],
            [ 0.0,  0.0,  n23],
            [ 0.0,  0.0,  n23]
        ], [
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0, -n23,  0.0],
            [ 2.0,  n43,  0.0],
            [-2.0, -n23, -2.0],
            [ 0.0, -n23,  0.0],
            [ 0.0, -n23,  0.0],
            [ 0.0,  n43,  2.0]
        ], [
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ n23,  0.0,  0.0],
            [ n23,  0.0,  0.0],
            [-n43,  0.0, -2.0],
            [-n43, -2.0,  0.0],
            [ n23,  0.0,  0.0],
            [ n23,  2.0,  2.0]
        ], [
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ n23, -n43, -n43],
            [ n23,  n23,  n23],
            [-n43,  n23, -n43],
            [-n43, -n43,  n23],
            [ n23,  n23,  n23],
            [ n23,  n23,  n23]
        ]])
    }
    fn calculate_standard_gradient_operators_transposed() -> StandardGradientOperatorsTransposed<M, O, P>
    {
        let n23: Scalar = 2.0/3.0;
        let n43: Scalar = 4.0/3.0;
        StandardGradientOperators::new([[
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
            [ 2.0,  0.0,  0.0],
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
            [-2.0, -2.0, -2.0],
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
            [ 0.0,  0.0,  0.0],
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
    fn get_integration_weight(&self) -> &Scalar
    {
        &INTEGRATION_WEIGHT
    }
    fn get_projected_gradient_vectors(&self) -> &ProjectedGradientVectors<G, N>
    {
        &self.projected_gradient_vectors
    }
}

super::composite_element_boilerplate!(Tetrahedron);