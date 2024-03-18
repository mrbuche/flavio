#[cfg(test)]
mod test;

use super::*;

const G: usize = 4;
const M: usize = 3;
const N: usize = 10;
const O: usize = 10;
const P: usize = 12;

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

impl<'a, C> CompositeElement<'a, C, G, M, N, O, P> for Tetrahedron<C>
where
    C: Constitutive<'a>
{
    fn calculate_projected_gradient_vectors(reference_nodal_coordinates: &ReferenceNodalCoordinates<N>) -> ProjectedGradientVectors<G, N>
    {
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