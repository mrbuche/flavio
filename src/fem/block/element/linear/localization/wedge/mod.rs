#[cfg(test)]
mod test;

use super::*;

const G: usize = 1;
const M: usize = 2;
const N: usize = 6;
const O: usize = 3;

pub struct Wedge<C>
{
    constitutive_models: [C; G],
    gradient_vectors: GradientVectors<N>,
    reference_normal: ReferenceNormal
}

impl<'a, C> FiniteElement<'a, C, G, N> for Wedge<C>
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
            reference_normal: Self::calculate_normal(&Self::calculate_midplane(&reference_nodal_coordinates))
        }
    }
}

impl<'a, C> LinearElement<'a, C, G, M, N, O> for Wedge<C>
where
    C: Constitutive<'a>
{
    fn calculate_deformation_gradient(&self, nodal_coordinates: &NodalCoordinates<N>) -> DeformationGradient
    {
        self.calculate_deformation_gradient_linear_localization_element(nodal_coordinates)
    }
    fn calculate_deformation_gradient_rate(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> DeformationGradientRate
    {
        self.calculate_deformation_gradient_rate_linear_localization_element(nodal_coordinates, nodal_velocities)
    }
    fn calculate_gradient_vectors(reference_nodal_coordinates: &ReferenceNodalCoordinates<N>) -> GradientVectors<N>
    {
        let reference_nodal_coordinates_midplane = Self::calculate_midplane(&reference_nodal_coordinates);
        let reference_dual_basis_vectors = Self::calculate_dual_basis(&reference_nodal_coordinates_midplane);
        let reference_normal = Self::calculate_normal(&reference_nodal_coordinates_midplane);
        let gradient_vectors_midplane = Self::calculate_standard_gradient_operator().iter()
        .map(|standard_gradient_operator_a|
            standard_gradient_operator_a.iter()
            .zip(reference_dual_basis_vectors.iter())
            .map(|(standard_gradient_operator_a_m, dual_reference_basis_vector_m)|
                dual_reference_basis_vector_m*standard_gradient_operator_a_m
            ).sum()
        ).collect::<GradientVectors<O>>();
        let mut gradient_vectors = GradientVectors::zero();
        gradient_vectors.iter_mut().take(O)
        .zip(gradient_vectors_midplane.iter())
        .for_each(|(gradient_vector_a, gradient_vector_midplane_a)|
            *gradient_vector_a = gradient_vector_midplane_a * 0.5 - &reference_normal / 3.0
        );
        gradient_vectors.iter_mut().skip(O)
        .zip(gradient_vectors_midplane.iter())
        .for_each(|(gradient_vector_a, gradient_vector_midplane_a)|
            *gradient_vector_a = gradient_vector_midplane_a * 0.5 + &reference_normal / 3.0
        );
        gradient_vectors
    }
    fn calculate_standard_gradient_operator() -> StandardGradientOperator<M, O>
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
}

impl<'a, C> LinearLocalizationElement<'a, C, G, M, N, O> for Wedge<C>
where
    C: Constitutive<'a>
{
    fn calculate_midplane<const I: usize>(nodal_coordinates: &Coordinates<I, N>) -> Coordinates<I, O>
    {
        nodal_coordinates.iter().skip(O)
        .zip(nodal_coordinates.iter().take(O))
        .map(|(nodal_coordinates_top, nodal_coordinates_bottom)|
            nodal_coordinates_top.iter()
            .zip(nodal_coordinates_bottom.iter())
            .map(|(nodal_coordinates_top_i, nodal_coordinates_bottom_i)|
                (nodal_coordinates_top_i + nodal_coordinates_bottom_i) * 0.5
            ).collect()
        ).collect()
    }
}

super::linear_surface_element_boilerplate!(Wedge);