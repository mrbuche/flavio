#[cfg(test)]
mod test;

pub mod surface;
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
        .map(|projected_gradient_vectors|
            nodal_coordinates.iter()
            .zip(projected_gradient_vectors.iter())
            .map(|(nodal_coordinate, projected_gradient_vector)|
                DeformationGradient::dyad(nodal_coordinate, projected_gradient_vector)
            ).sum()
        ).collect()
    }
    fn calculate_deformation_gradient_rates(&self, _: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> DeformationGradientRates<G>
    {
        self.get_projected_gradient_vectors().iter()
        .map(|projected_gradient_vectors|
            nodal_velocities.iter()
            .zip(projected_gradient_vectors.iter())
            .map(|(nodal_velocity, projected_gradient_vector)|
                DeformationGradientRate::dyad(nodal_velocity, projected_gradient_vector)
            ).sum()
        ).collect()
    }
    fn calculate_inverse_normalized_projection_matrix() -> NormalizedProjectionMatrix<Q>;
    fn calculate_jacobians_and_parametric_gradient_operators(reference_nodal_coordinates: &ReferenceNodalCoordinates<N>) -> (Scalars<P>, ParametricGradientOperators<P>);
    fn calculate_projected_gradient_vectors(reference_nodal_coordinates: &ReferenceNodalCoordinates<N>) -> ProjectedGradientVectors<G, N>;
    fn calculate_scaled_composite_jacobian_at_integration_points(reference_nodal_coordinates: &ReferenceNodalCoordinates<N>) -> Scalars<G>
    {
        let (jacobians, _) = Self::calculate_jacobians_and_parametric_gradient_operators(reference_nodal_coordinates);
        let vector = Self::calculate_inverse_normalized_projection_matrix() *
        Self::calculate_shape_function_integrals().iter()
        .zip(jacobians.iter())
        .map(|(shape_function_integral, jacobian)|
            shape_function_integral * jacobian
        ).sum::<TensorRank1<Q, 9>>();
        Self::calculate_shape_functions_at_integration_points().iter()
        .map(|shape_functions_at_integration_point|
            (shape_functions_at_integration_point * &vector) * Self::get_integration_weight()
        ).collect()
    }
    fn calculate_shape_function_integrals() -> ShapeFunctionIntegrals<P, Q>;
    fn calculate_shape_function_integrals_products() -> ShapeFunctionIntegralsProducts<P, Q>;
    fn calculate_shape_functions_at_integration_points() -> ShapeFunctionsAtIntegrationPoints<G, Q>;
    fn calculate_standard_gradient_operators() -> StandardGradientOperators<M, O, P>;
    fn calculate_standard_gradient_operators_transposed() -> StandardGradientOperatorsTransposed<M, O, P>;
    fn get_constitutive_models(&self) -> &[C; G];
    fn get_integration_weight() -> Scalar;
    fn get_projected_gradient_vectors(&self) -> &ProjectedGradientVectors<G, N>;
    fn get_scaled_composite_jacobians(&self) -> &Scalars<G>;
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
        .zip(self.get_projected_gradient_vectors().iter()
        .zip(self.get_scaled_composite_jacobians().iter()))
        .map(|(first_piola_kirchoff_stress, (projected_gradient_vectors, scaled_composite_jacobian))|
            projected_gradient_vectors.iter()
            .map(|projected_gradient_vector|
                (first_piola_kirchoff_stress * projected_gradient_vector) * scaled_composite_jacobian
            ).collect()
        ).sum::<NodalForces<N>>()
    }
    fn calculate_nodal_stiffnesses_composite_element(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalStiffnesses<N>
    {
        self.get_constitutive_models().iter()
        .zip(self.calculate_deformation_gradients(nodal_coordinates).iter())
        .map(|(constitutive_model, deformation_gradient)|
            constitutive_model.calculate_first_piola_kirchoff_tangent_stiffness(deformation_gradient)
        ).collect::<FirstPiolaKirchoffTangentStiffnesses<G>>().iter()
        .zip(self.get_projected_gradient_vectors().iter()
        .zip(self.get_projected_gradient_vectors().iter()
        .zip(self.get_scaled_composite_jacobians().iter())))
        .map(|(first_piola_kirchoff_tangent_stiffness, (projected_gradient_vectors_a, (projected_gradient_vectors_b, scaled_composite_jacobian)))|
            projected_gradient_vectors_a.iter()
            .map(|projected_gradient_vector_a|
                projected_gradient_vectors_b.iter()
                .map(|projected_gradient_vector_b|
                    first_piola_kirchoff_tangent_stiffness
                    .contract_second_fourth_indices_with_first_indices_of(
                        projected_gradient_vector_a, projected_gradient_vector_b
                    ) * scaled_composite_jacobian
                ).collect()
            ).collect()
        ).sum::<NodalStiffnesses<N>>()
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
        .zip(self.calculate_deformation_gradients(nodal_coordinates).iter()
        .zip(self.get_scaled_composite_jacobians().iter()))
        .map(|(constitutive_model, (deformation_gradient, scaled_composite_jacobian))|
            constitutive_model.calculate_helmholtz_free_energy_density(deformation_gradient) * scaled_composite_jacobian
        ).sum::<Scalar>()
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
        .zip(self.get_projected_gradient_vectors().iter()
        .zip(self.get_scaled_composite_jacobians().iter()))
        .map(|(first_piola_kirchoff_stress, (projected_gradient_vectors, scaled_composite_jacobian))|
            projected_gradient_vectors.iter()
            .map(|projected_gradient_vector|
                (first_piola_kirchoff_stress * projected_gradient_vector) * scaled_composite_jacobian
            ).collect()
        ).sum::<NodalForces<N>>()
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
        .zip(self.get_projected_gradient_vectors().iter()
        .zip(self.get_scaled_composite_jacobians().iter())))
        .map(|(first_piola_kirchoff_rate_tangent_stiffness, (projected_gradient_vectors_a, (projected_gradient_vectors_b, scaled_composite_jacobian)))|
            projected_gradient_vectors_a.iter()
            .map(|projected_gradient_vector_a|
                projected_gradient_vectors_b.iter()
                .map(|projected_gradient_vector_b|
                    first_piola_kirchoff_rate_tangent_stiffness
                    .contract_second_fourth_indices_with_first_indices_of(
                        projected_gradient_vector_a, projected_gradient_vector_b
                    ) * scaled_composite_jacobian
                ).collect()
            ).collect()
        ).sum::<NodalStiffnesses<N>>()
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
        .zip(self.calculate_deformation_gradient_rates(nodal_coordinates, nodal_velocities).iter()
        .zip(self.get_scaled_composite_jacobians().iter())))
        .map(|(constitutive_model, (deformation_gradient, (deformation_gradient_rate, scaled_composite_jacobian)))|
            constitutive_model.calculate_viscous_dissipation(deformation_gradient, deformation_gradient_rate) * scaled_composite_jacobian
        ).sum::<Scalar>()
    }
    fn calculate_dissipation_potential_composite_element(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> Scalar
    {
        self.get_constitutive_models().iter()
        .zip(self.calculate_deformation_gradients(nodal_coordinates).iter()
        .zip(self.calculate_deformation_gradient_rates(nodal_coordinates, nodal_velocities).iter()
        .zip(self.get_scaled_composite_jacobians().iter())))
        .map(|(constitutive_model, (deformation_gradient, (deformation_gradient_rate, scaled_composite_jacobian)))|
            constitutive_model.calculate_dissipation_potential(deformation_gradient, deformation_gradient_rate) * scaled_composite_jacobian
        ).sum::<Scalar>()
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
        .zip(self.calculate_deformation_gradients(nodal_coordinates).iter()
        .zip(self.get_scaled_composite_jacobians().iter()))
        .map(|(constitutive_model, (deformation_gradient, scaled_composite_jacobian))|
            constitutive_model.calculate_helmholtz_free_energy_density(deformation_gradient) * scaled_composite_jacobian
        ).sum::<Scalar>()
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