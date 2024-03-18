#[cfg(test)]
mod test;

pub mod tetrahedron;

use super::*;

pub trait CompositeElement<'a, C, const G: usize, const M: usize, const N: usize, const O: usize, const P: usize>
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
    fn calculate_standard_gradient_operators() -> StandardGradientOperators<M, O, P>;
    fn get_constitutive_models(&self) -> &[C; G];
    fn get_integration_weight(&self) -> &Scalar;
    fn get_projected_gradient_vectors(&self) -> &ProjectedGradientVectors<G, N>;
}

pub trait ElasticCompositeElement<'a, C, const G: usize, const M: usize, const N: usize, const O: usize, const P: usize>
where
    C: Elastic<'a>,
    Self: CompositeElement<'a, C, G, M, N, O, P>
{
    fn calculate_nodal_forces_composite_element(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalForces<N>
    {
        self.get_constitutive_models().iter()
        .zip(self.calculate_deformation_gradients(nodal_coordinates).iter())
        .map(|(constitutive_model, deformation_gradient)|
            constitutive_model.calculate_first_piola_kirchoff_stress(
                deformation_gradient
            ) / self.get_integration_weight()
        ).collect::<FirstPiolaKirchoffStresses<G>>()
        .iter()
        .zip(self.get_projected_gradient_vectors().iter())
        .map(|(first_piola_kirchoff_stress, projected_gradient_vectors)|
            projected_gradient_vectors.iter()
            .map(|projected_gradient_vector|
                first_piola_kirchoff_stress * projected_gradient_vector
            ).collect()
        ).sum()
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

pub trait HyperelasticCompositeElement<'a, C, const G: usize, const M: usize, const N: usize, const O: usize, const P: usize>
where
    C: Hyperelastic<'a>,
    Self: ElasticCompositeElement<'a, C, G, M, N, O, P>
{
    fn calculate_helmholtz_free_energy_composite_element(&self, nodal_coordinates: &NodalCoordinates<N>) -> Scalar
    {
        self.get_constitutive_models().iter()
        .zip(self.calculate_deformation_gradients(nodal_coordinates).iter())
        .map(|(constitutive_model, deformation_gradient)|
            constitutive_model.calculate_helmholtz_free_energy_density(
                deformation_gradient
            )
        ).sum::<Scalar>() / self.get_integration_weight()
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
        impl<'a, C> ElasticCompositeElement<'a, C, G, M, N, O, P> for $element<C>
        where
            C: Elastic<'a>
        {}
        impl<'a, C> HyperelasticFiniteElement<'a, C, G, N> for $element<C>
        where
            C: Hyperelastic<'a>
        {
            fn calculate_helmholtz_free_energy(&self, nodal_coordinates: &NodalCoordinates<N>) -> Scalar
            {
                self.calculate_helmholtz_free_energy_composite_element(nodal_coordinates)
            }
        }
        impl<'a, C> HyperelasticCompositeElement<'a, C, G, M, N, O, P> for $element<C>
        where
            C: Hyperelastic<'a>
        {}
    }
}
pub(crate) use composite_element_boilerplate;