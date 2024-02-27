#[cfg(test)]
mod test;

use crate::math::TensorRank2;

use super::*;

const G: usize = 1;
const M: usize = 2;
const N: usize = 3;

pub struct Triangle<C>
{
    constitutive_models: [C; G],
    gradient_vectors: GradientVectors<N>
}

impl<'a, C> FiniteElement<'a, C, G, N> for Triangle<C>
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
            gradient_vectors: Self::calculate_gradient_vectors(&reference_nodal_coordinates)
        }
    }
}

impl<'a, C> LinearFiniteElement<'a, C, G, M, N> for Triangle<C>
where
    C: Constitutive<'a>
{
    fn calculate_deformation_gradient(&self, nodal_coordinates: &NodalCoordinates<N>) -> DeformationGradient
    {
        nodal_coordinates.iter()
        .zip(self.get_gradient_vectors().iter())
        .map(|(current_nodal_coordinate, gradient_vector)|
            DeformationGradient::dyad(current_nodal_coordinate, gradient_vector)
        ).sum::<DeformationGradient>()

        + DeformationGradient::new([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    }
    fn calculate_gradient_vectors(reference_nodal_coordinates: &ReferenceNodalCoordinates<N>) -> GradientVectors<N>
    {
        let standard_gradient_operator = Self::calculate_standard_gradient_operator();
        let basis_vectors: Vectors<0, 2> =
            standard_gradient_operator.iter()
            .zip(reference_nodal_coordinates.iter())
            .map(|(standard_gradient_operator_a, reference_nodal_coordinates_a)|
                standard_gradient_operator_a.iter()
                .map(|standard_gradient_operator_a_m|
                    reference_nodal_coordinates_a.iter()
                    .map(|reference_nodal_coordinates_a_i|
                        reference_nodal_coordinates_a_i*standard_gradient_operator_a_m
                    ).collect()
                ).collect()
            ).sum();
        let metric_tensor: TensorRank2<2, 0, 0> =
            basis_vectors.iter()
            .map(|basis_vector_m|
                basis_vectors.iter()
                .map(|basis_vector_n|
                    basis_vector_m*basis_vector_n
                ).collect()
            ).collect::<TensorRank2<2, 0, 0>>()
            .inverse();
        let dual_basis_vectors: Vectors<0, 2> =
            metric_tensor.iter()
            .map(|metric_tensor_m|
                metric_tensor_m.iter()
                .zip(basis_vectors.iter())
                .map(|(metric_tensor_mn, basis_vector_n)|
                    basis_vector_n*metric_tensor_mn
                ).sum()
            ).collect();
        standard_gradient_operator.iter()
        .map(|standard_gradient_operator_a|
            standard_gradient_operator_a.iter()
            .zip(dual_basis_vectors.iter())
            .map(|(standard_gradient_operator_a_m, dual_basis_vector_m)|
                dual_basis_vector_m*standard_gradient_operator_a_m
            ).sum()
        ).collect()
    }
    fn calculate_standard_gradient_operator() -> StandardGradientOperator<M, N>
    {
        StandardGradientOperator::new([
            [ 1.0,  0.0],
            [ 0.0,  1.0],
            [-1.0, -1.0],
        ])
    }
    fn get_gradient_vectors(&self) -> &GradientVectors<N>
    {
        &self.gradient_vectors
    }
}