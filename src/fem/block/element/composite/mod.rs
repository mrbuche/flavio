#[cfg(test)]
mod test;

pub mod tetrahedron;

use super::*;

pub trait CompositeElement<'a, C, const G: usize, const M: usize, const N: usize, const O: usize, const P: usize, const Q: usize>
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
    fn calculate_shape_function_integrals() -> ShapeFunctionIntegrals<P, Q>;
    fn calculate_shape_function_integrals_products() -> ShapeFunctionIntegralsProducts<P, Q>;
    fn calculate_shape_functions_at_integration_points() -> ShapeFunctionsAtIntegrationPoints<G, Q>;
    fn calculate_standard_gradient_operators() -> StandardGradientOperators<M, O, P>;
    fn calculate_standard_gradient_operators_transposed() -> StandardGradientOperatorsTransposed<M, O, P>;
    fn get_constitutive_models(&self) -> &[C; G];
    fn get_integration_weight(&self) -> &Scalar;
    fn get_projected_gradient_vectors(&self) -> &ProjectedGradientVectors<G, N>;
}

pub trait ElasticCompositeElement<'a, C, const G: usize, const M: usize, const N: usize, const O: usize, const P: usize, const Q: usize>
where
    C: Elastic<'a>,
    Self: CompositeElement<'a, C, G, M, N, O, P, Q>
{
    fn calculate_nodal_forces_composite_element(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalForces<N>
    {
        self.get_constitutive_models().iter()
        .zip(self.calculate_deformation_gradients(nodal_coordinates).iter())
        .map(|(constitutive_model, deformation_gradient)|
            constitutive_model.calculate_first_piola_kirchoff_stress(deformation_gradient)
        ).collect::<FirstPiolaKirchoffStresses<G>>().iter()
        .zip(self.get_projected_gradient_vectors().iter())
        .map(|(first_piola_kirchoff_stress, projected_gradient_vectors)|
            projected_gradient_vectors.iter()
            .map(|projected_gradient_vector|
                first_piola_kirchoff_stress * projected_gradient_vector
            ).collect()
        ).sum::<NodalForces<N>>() / self.get_integration_weight()
    }
    fn calculate_nodal_stiffnesses_composite_element(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalStiffnesses<N>
    {
        self.get_constitutive_models().iter()
        .zip(self.calculate_deformation_gradients(nodal_coordinates).iter())
        .map(|(constitutive_model, deformation_gradient)|
            constitutive_model.calculate_first_piola_kirchoff_tangent_stiffness(deformation_gradient)
        ).collect::<FirstPiolaKirchoffTangentStiffnesses<G>>().iter()
        .zip(self.get_projected_gradient_vectors().iter()
        .zip(self.get_projected_gradient_vectors().iter()))
        .map(|(first_piola_kirchoff_tangent_stiffness, (projected_gradient_vectors_a, projected_gradient_vectors_b))|
            projected_gradient_vectors_a.iter()
            .map(|projected_gradient_vector_a|
                projected_gradient_vectors_b.iter()
                .map(|projected_gradient_vector_b|
                    first_piola_kirchoff_tangent_stiffness
                    .contract_second_fourth_indices_with_first_indices_of(
                        projected_gradient_vector_a, projected_gradient_vector_b
                    )
                ).collect()
            ).collect()
        ).sum::<NodalStiffnesses<N>>() / self.get_integration_weight()
    }
}

pub trait HyperelasticCompositeElement<'a, C, const G: usize, const M: usize, const N: usize, const O: usize, const P: usize, const Q: usize>
where
    C: Hyperelastic<'a>,
    Self: ElasticCompositeElement<'a, C, G, M, N, O, P, Q>
{
    fn calculate_helmholtz_free_energy_composite_element(&self, nodal_coordinates: &NodalCoordinates<N>) -> Scalar
    {
        self.get_constitutive_models().iter()
        .zip(self.calculate_deformation_gradients(nodal_coordinates).iter())
        .map(|(constitutive_model, deformation_gradient)|
            constitutive_model.calculate_helmholtz_free_energy_density(deformation_gradient)
        ).sum::<Scalar>() / self.get_integration_weight()
    }
}

pub trait ViscoelasticCompositeElement<'a, C, const G: usize, const M: usize, const N: usize, const O: usize, const P: usize, const Q: usize>
where
    C: Viscoelastic<'a>,
    Self: CompositeElement<'a, C, G, M, N, O, P, Q>
{
    fn calculate_nodal_forces_composite_element(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> NodalForces<N>
    {
        self.get_constitutive_models().iter()
        .zip(self.calculate_deformation_gradients(nodal_coordinates).iter()
        .zip(self.calculate_deformation_gradient_rates(nodal_coordinates, nodal_velocities).iter()))
        .map(|(constitutive_model, (deformation_gradient, deformation_gradient_rate))|
            constitutive_model.calculate_first_piola_kirchoff_stress(deformation_gradient, deformation_gradient_rate)
        ).collect::<FirstPiolaKirchoffStresses<G>>().iter()
        .zip(self.get_projected_gradient_vectors().iter())
        .map(|(first_piola_kirchoff_stress, projected_gradient_vectors)|
            projected_gradient_vectors.iter()
            .map(|projected_gradient_vector|
                first_piola_kirchoff_stress * projected_gradient_vector
            ).collect()
        ).sum::<NodalForces<N>>() / self.get_integration_weight()
    }
    fn calculate_nodal_stiffnesses_composite_element(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> NodalStiffnesses<N>
    {
        self.get_constitutive_models().iter()
        .zip(self.calculate_deformation_gradients(nodal_coordinates).iter()
        .zip(self.calculate_deformation_gradient_rates(nodal_coordinates, nodal_velocities).iter()))
        .map(|(constitutive_model, (deformation_gradient, deformation_gradient_rate))|
            constitutive_model.calculate_first_piola_kirchoff_rate_tangent_stiffness(deformation_gradient, deformation_gradient_rate)
        ).collect::<FirstPiolaKirchoffRateTangentStiffnesses<G>>().iter()
        .zip(self.get_projected_gradient_vectors().iter()
        .zip(self.get_projected_gradient_vectors().iter()))
        .map(|(first_piola_kirchoff_rate_tangent_stiffness, (projected_gradient_vectors_a, projected_gradient_vectors_b))|
            projected_gradient_vectors_a.iter()
            .map(|projected_gradient_vector_a|
                projected_gradient_vectors_b.iter()
                .map(|projected_gradient_vector_b|
                    first_piola_kirchoff_rate_tangent_stiffness
                    .contract_second_fourth_indices_with_first_indices_of(
                        projected_gradient_vector_a, projected_gradient_vector_b
                    )
                ).collect()
            ).collect()
        ).sum::<NodalStiffnesses<N>>() / self.get_integration_weight()
    }
}

pub trait ElasticHyperviscousCompositeElement<'a, C, const G: usize, const M: usize, const N: usize, const O: usize, const P: usize, const Q: usize>
where
    C: ElasticHyperviscous<'a>,
    Self: ViscoelasticCompositeElement<'a, C, G, M, N, O, P, Q>
{
    fn calculate_viscous_dissipation_composite_element(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> Scalar
    {
        self.get_constitutive_models().iter()
        .zip(self.calculate_deformation_gradients(nodal_coordinates).iter()
        .zip(self.calculate_deformation_gradient_rates(nodal_coordinates, nodal_velocities).iter()))
        .map(|(constitutive_model, (deformation_gradient, deformation_gradient_rate))|
            constitutive_model.calculate_viscous_dissipation(deformation_gradient, deformation_gradient_rate)
        ).sum::<Scalar>() / self.get_integration_weight()
    }
    fn calculate_dissipation_potential_composite_element(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> Scalar
    {
        self.get_constitutive_models().iter()
        .zip(self.calculate_deformation_gradients(nodal_coordinates).iter()
        .zip(self.calculate_deformation_gradient_rates(nodal_coordinates, nodal_velocities).iter()))
        .map(|(constitutive_model, (deformation_gradient, deformation_gradient_rate))|
            constitutive_model.calculate_dissipation_potential(deformation_gradient, deformation_gradient_rate)
        ).sum::<Scalar>() / self.get_integration_weight()
    }
}

pub trait HyperviscoelasticCompositeElement<'a, C, const G: usize, const M: usize, const N: usize, const O: usize, const P: usize, const Q: usize>
where
    C: Hyperviscoelastic<'a>,
    Self: ElasticHyperviscousCompositeElement<'a, C, G, M, N, O, P, Q>
{
    fn calculate_helmholtz_free_energy_composite_element(&self, nodal_coordinates: &NodalCoordinates<N>) -> Scalar
    {
        self.get_constitutive_models().iter()
        .zip(self.calculate_deformation_gradients(nodal_coordinates).iter())
        .map(|(constitutive_model, deformation_gradient)|
            constitutive_model.calculate_helmholtz_free_energy_density(deformation_gradient)
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
        impl<'a, C> ElasticCompositeElement<'a, C, G, M, N, O, P, Q> for $element<C>
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
        impl<'a, C> HyperelasticCompositeElement<'a, C, G, M, N, O, P, Q> for $element<C>
        where
            C: Hyperelastic<'a>
        {}
        impl<'a, C> ViscoelasticFiniteElement<'a, C, G, N> for $element<C>
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
        impl<'a, C> ViscoelasticCompositeElement<'a, C, G, M, N, O, P, Q> for $element<C>
        where
            C: Viscoelastic<'a>
        {}
        impl<'a, C> ElasticHyperviscousFiniteElement<'a, C, G, N> for $element<C>
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
        impl<'a, C> ElasticHyperviscousCompositeElement<'a, C, G, M, N, O, P, Q> for $element<C>
        where
            C: ElasticHyperviscous<'a>
        {}
        impl<'a, C> HyperviscoelasticFiniteElement<'a, C, G, N> for $element<C>
        where
            C: Hyperviscoelastic<'a>
        {
            fn calculate_helmholtz_free_energy(&self, nodal_coordinates: &NodalCoordinates<N>) -> Scalar
            {
                self.calculate_helmholtz_free_energy_composite_element(nodal_coordinates)
            }
        }
        impl<'a, C> HyperviscoelasticCompositeElement<'a, C, G, M, N, O, P, Q> for $element<C>
        where
            C: Hyperviscoelastic<'a>
        {}
    }
}
pub(crate) use composite_element_boilerplate;