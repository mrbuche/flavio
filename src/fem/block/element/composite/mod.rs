#[cfg(test)]
mod test;

pub mod tetrahedron;

use super::*;

type ProjectedGradientVectors<const G: usize, const N: usize> = Vectors2D<0, N, G>;

pub trait CompositeElement<'a, C, const G: usize, const M: usize, const N: usize, const O: usize>
where
    C: Constitutive<'a>,
    Self: FiniteElement<'a, C, G, N>
{
    fn calculate_deformation_gradients(&self, nodal_coordinates: &NodalCoordinates<N>) -> DeformationGradients<G>
    {
        self.get_projected_gradient_vectors().iter()
        .map(|projected_gradient_vectors_g|
            nodal_coordinates.iter()
            .zip(projected_gradient_vectors_g.iter())
            .map(|(nodal_coordinate_a, projected_gradient_vector_g_a)|
                DeformationGradient::dyad(nodal_coordinate_a, projected_gradient_vector_g_a)
            ).sum()
        ).collect()
    }
    fn calculate_deformation_gradient_rates(&self, _: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> DeformationGradientRates<G>
    {
        self.get_projected_gradient_vectors().iter()
        .map(|projected_gradient_vectors_g|
            nodal_velocities.iter()
            .zip(projected_gradient_vectors_g.iter())
            .map(|(nodal_velocity_a, projected_gradient_vector_g_a)|
                DeformationGradientRate::dyad(nodal_velocity_a, projected_gradient_vector_g_a)
            ).sum()
        ).collect()
    }
    fn calculate_projected_gradient_vectors(reference_nodal_coordinates: &ReferenceNodalCoordinates<N>) -> ProjectedGradientVectors<G, N>;
    fn get_constitutive_models(&self) -> &[C; G];
    fn get_integration_weights(&self) -> &IntegrationWeights<G>;
    fn get_projected_gradient_vectors(&self) -> &ProjectedGradientVectors<G, N>;
}

pub trait ElasticCompositeElement<'a, C, const G: usize, const M: usize, const N: usize, const O: usize>
where
    C: Elastic<'a>,
    Self: CompositeElement<'a, C, G, M, N, O>
{
    fn calculate_nodal_forces_composite_element(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalForces<N>
    {
        self.get_constitutive_models().iter()
        .zip(self.calculate_deformation_gradients(nodal_coordinates).iter()
        .zip(self.get_integration_weights().iter()))
        .map(|(constitutive_model, (deformation_gradient, integration_weight))|
            constitutive_model.calculate_first_piola_kirchoff_stress(
                deformation_gradient
            ) * integration_weight
        ).collect::<FirstPiolaKirchoffStresses<G>>()
        .iter()
        .zip(self.get_projected_gradient_vectors().iter())
        .map(|(first_piola_kirchoff_stress, projected_gradient_vectors)|
            projected_gradient_vectors.iter()
            .map(|projected_gradient_vector|
                first_piola_kirchoff_stress * projected_gradient_vector
            ).collect()
        ).sum()
        // let first_piola_kirchoff_stress = self.get_constitutive_models()[0]
        // .calculate_first_piola_kirchoff_stress(
        //     &self.calculate_deformation_gradient(nodal_coordinates)
        // );
        // self.get_gradient_vectors().iter()
        // .map(|gradient_vector|
        //     &first_piola_kirchoff_stress * gradient_vector
        // ).collect()
    }
    fn calculate_nodal_stiffnesses_composite_element(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalStiffnesses<N>
    {
        todo!()
        // let first_piola_kirchoff_tangent_stiffness = self.get_constitutive_models()[0]
        // .calculate_first_piola_kirchoff_tangent_stiffness(
        //     &self.calculate_deformation_gradient(nodal_coordinates)
        // );
        // let gradient_vectors = self.get_gradient_vectors();
        // gradient_vectors.iter()
        // .map(|gradient_vector_a|
        //     gradient_vectors.iter()
        //     .map(|gradient_vector_b|
        //         first_piola_kirchoff_tangent_stiffness
        //         .contract_second_fourth_indices_with_first_indices_of(
        //             gradient_vector_a, gradient_vector_b
        //         )
        //     ).collect()
        // ).collect()
    }
}

macro_rules! composite_element_boilerplate
{
    ($element: ident) =>
    {
        impl<'a, C> ElasticFiniteElement<'a, C, G, N> for $element<C>
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
        impl<'a, C> ElasticCompositeElement<'a, C, G, M, N, O> for $element<C>
        where
            C: Elastic<'a>
        {}
    }
}
pub(crate) use composite_element_boilerplate;