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
        let standard_gradient_operators = Self::calculate_standard_gradient_operators();
        let parametric_gradient_operators =
        standard_gradient_operators.iter()
        .map(|standard_gradient_operator|
            reference_nodal_coordinates * standard_gradient_operator
        ).collect::<ParametricGradientOperators<P>>();
        let jacobians =
        parametric_gradient_operators.iter()
        .map(|parametric_gradient_operator|
            parametric_gradient_operator.determinant()
        ).collect::<Scalars<P>>();
        //
        // need to loop/collect over _a of standard_gradient_operators as well
        //
        // Self::calculate_shape_function_integrals().iter()
        // .zip(standard_gradient_operators.iter()
        // .zip(parametric_gradient_operators.iter()
        // .zip(jacobians.iter())))
        // .map(|(shape_function_integral, (standard_gradient_operator, (parametric_gradient_operator, jacobian)))|
        //     (standard_gradient_operator * parametric_gradient_operator.inverse_transpose()) * jacobian
        // ).sum();
        // let gradient_operators = 
        // Self::calculate_standard_gradient_operators().iter()
        // .map(|standard_gradient_operator|
        //     (reference_nodal_coordinates * standard_gradient_operator).inverse_transpose() * standard_gradient_operator
        // ).sum();
        // (0..G).map(|_|
        //     gradient_operator
        // ).collect()
        todo!()
        // evaluating the lambda shape function at each integration point gives you the <G>
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
            [1966080.0, 368640.0, 368640.0, 368640.0],
            [ 368640.0, 122880.0,  61440.0,  61440.0],
            [ 368640.0,  61440.0, 122880.0,  61440.0],
            [ 368640.0,  61440.0,  61440.0, 122880.0]
        ], [
            [122880.0,  368640.0,  61440.0,  61440.0],
            [368640.0, 1966080.0, 368640.0, 368640.0],
            [ 61440.0,  368640.0, 122880.0,  61440.0],
            [ 61440.0,  368640.0,  61440.0, 122880.0]
        ], [
            [122880.0,  61440.0,  368640.0,  61440.0],
            [ 61440.0, 122880.0,  368640.0,  61440.0],
            [368640.0, 368640.0, 1966080.0, 368640.0],
            [ 61440.0,  61440.0,  368640.0, 122880.0]
        ], [
            [122880.0,  61440.0,  61440.0,  368640.0],
            [ 61440.0, 122880.0,  61440.0,  368640.0],
            [ 61440.0,  61440.0, 122880.0,  368640.0],
            [368640.0, 368640.0, 368640.0, 1966080.0]
        ], [
            [107520.0, 199680.0,  76800.0,  76800.0],
            [199680.0, 476160.0, 199680.0, 199680.0],
            [ 76800.0, 199680.0, 107520.0,  76800.0],
            [ 76800.0, 199680.0,  76800.0, 107520.0]
        ], [
            [15360.0,  46080.0,  46080.0,  46080.0],
            [46080.0, 261120.0, 230400.0, 230400.0],
            [46080.0, 230400.0, 261120.0, 230400.0],
            [46080.0, 230400.0, 230400.0, 261120.0]
        ], [
            [107520.0,  76800.0,  76800.0, 199680.0],
            [ 76800.0, 107520.0,  76800.0, 199680.0],
            [ 76800.0,  76800.0, 107520.0, 199680.0],
            [199680.0, 199680.0, 199680.0, 476160.0]
        ], [
            [261120.0, 230400.0, 46080.0, 230400.0],
            [230400.0, 261120.0, 46080.0, 230400.0],
            [ 46080.0,  46080.0, 15360.0,  46080.0],
            [230400.0, 230400.0, 46080.0, 261120.0]
        ], [
            [261120.0, 230400.0, 230400.0, 46080.0],
            [230400.0, 261120.0, 230400.0, 46080.0],
            [230400.0, 230400.0, 261120.0, 46080.0],
             [46080.0,  46080.0,  46080.0, 15360.0]
        ], [
            [107520.0,  76800.0, 199680.0,  76800.0],
            [ 76800.0, 107520.0, 199680.0,  76800.0],
            [199680.0, 199680.0, 476160.0, 199680.0],
            [ 76800.0,  76800.0, 199680.0, 107520.0]
        ], [
            [261120.0, 46080.0, 230400.0, 230400.0],
            [ 46080.0, 15360.0,  46080.0,  46080.0],
            [230400.0, 46080.0, 261120.0, 230400.0],
            [230400.0, 46080.0, 230400.0, 261120.0]
        ], [
            [476160.0, 199680.0, 199680.0, 199680.0],
            [199680.0, 107520.0,  76800.0,  76800.0],
            [199680.0,  76800.0, 107520.0,  76800.0],
            [199680.0,  76800.0,  76800.0, 107520.0]
        ]])
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