#[cfg(test)]
pub mod test;

pub mod localization;
pub mod surface;
pub mod tetrahedron;

use super::*;

pub trait LinearElement<'a, C, const G: usize, const M: usize, const N: usize, const O: usize>
where
    C: Constitutive<'a>
{
    fn calculate_deformation_gradient(&self, nodal_coordinates: &NodalCoordinates<N>) -> DeformationGradient
    {
        nodal_coordinates.iter()
        .zip(self.get_gradient_vectors().iter())
        .map(|(nodal_coordinate, gradient_vector)|
            DeformationGradient::dyad(nodal_coordinate, gradient_vector)
        ).sum()
    }
    fn calculate_deformation_gradient_rate(&self, _: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> DeformationGradientRate
    {
        nodal_velocities.iter()
        .zip(self.get_gradient_vectors().iter())
        .map(|(nodal_velocity, gradient_vector)|
            DeformationGradientRate::dyad(nodal_velocity, gradient_vector)
        ).sum()
    }
    fn calculate_gradient_vectors(_reference_nodal_coordinates: &ReferenceNodalCoordinates<O>) -> GradientVectors<N>
    {
        panic!()
    }
    fn calculate_reference_jacobian(reference_nodal_coordinates: &ReferenceNodalCoordinates<O>) -> Scalar;
    fn calculate_standard_gradient_operator() -> StandardGradientOperator<M, O>;
    fn get_constitutive_model(&self) -> &C;
    fn get_gradient_vectors(&self) -> &GradientVectors<N>;
    fn get_integration_weight(&self) -> &Scalar;
}

pub trait ElasticLinearElement<'a, C, const G: usize, const M: usize, const N: usize, const O: usize>
where
    C: Elastic<'a>,
    Self: LinearElement<'a, C, G, M, N, O>
{
    fn calculate_nodal_forces_linear_element(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalForces<N>
    {
        let first_piola_kirchoff_stress = self.get_constitutive_model()
        .calculate_first_piola_kirchoff_stress(
            &self.calculate_deformation_gradient(nodal_coordinates)
        );
        self.get_gradient_vectors().iter()
        .map(|gradient_vector|
            (&first_piola_kirchoff_stress * gradient_vector) * self.get_integration_weight()
        ).collect()
    }
    fn calculate_nodal_stiffnesses_linear_element(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalStiffnesses<N>
    {
        let first_piola_kirchoff_tangent_stiffness = self.get_constitutive_model()
        .calculate_first_piola_kirchoff_tangent_stiffness(
            &self.calculate_deformation_gradient(nodal_coordinates)
        );
        let gradient_vectors = self.get_gradient_vectors();
        gradient_vectors.iter()
        .map(|gradient_vector_a|
            gradient_vectors.iter()
            .map(|gradient_vector_b|
                first_piola_kirchoff_tangent_stiffness
                .contract_second_fourth_indices_with_first_indices_of(
                    gradient_vector_a, gradient_vector_b
                ) * self.get_integration_weight()
            ).collect()
        ).collect()
    }
}

pub trait HyperelasticLinearElement<'a, C, const G: usize, const M: usize, const N: usize, const O: usize>
where
    C: Hyperelastic<'a>,
    Self: ElasticLinearElement<'a, C, G, M, N, O>
{
    fn calculate_helmholtz_free_energy_linear_element(&self, nodal_coordinates: &NodalCoordinates<N>) -> Scalar
    {
        self.get_constitutive_model()
        .calculate_helmholtz_free_energy_density(
            &self.calculate_deformation_gradient(nodal_coordinates)
        ) * self.get_integration_weight()
    }
}

pub trait ViscoelasticLinearElement<'a, C, const G: usize, const M: usize, const N: usize, const O: usize>
where
    C: Viscoelastic<'a>,
    Self: LinearElement<'a, C, G, M, N, O>
{
    fn calculate_nodal_forces_linear_element(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> NodalForces<N>
    {
        let first_piola_kirchoff_stress = self.get_constitutive_model()
        .calculate_first_piola_kirchoff_stress(
            &self.calculate_deformation_gradient(nodal_coordinates),
            &self.calculate_deformation_gradient_rate(nodal_coordinates, nodal_velocities)
        );
        self.get_gradient_vectors().iter()
        .map(|gradient_vector|
            (&first_piola_kirchoff_stress * gradient_vector) * self.get_integration_weight()
        ).collect()
    }
    fn calculate_nodal_stiffnesses_linear_element(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> NodalStiffnesses<N>
    {
        let first_piola_kirchoff_tangent_stiffness = self.get_constitutive_model()
        .calculate_first_piola_kirchoff_rate_tangent_stiffness(
            &self.calculate_deformation_gradient(nodal_coordinates),
            &self.calculate_deformation_gradient_rate(nodal_coordinates, nodal_velocities)
        );
        let gradient_vectors = self.get_gradient_vectors();
        gradient_vectors.iter()
        .map(|gradient_vector_a|
            gradient_vectors.iter()
            .map(|gradient_vector_b|
                first_piola_kirchoff_tangent_stiffness
                .contract_second_fourth_indices_with_first_indices_of(
                    gradient_vector_a, gradient_vector_b
                ) * self.get_integration_weight()
            ).collect()
        ).collect()
    }
}

pub trait ElasticHyperviscousLinearElement<'a, C, const G: usize, const M: usize, const N: usize, const O: usize>
where
    C: ElasticHyperviscous<'a>,
    Self: ViscoelasticLinearElement<'a, C, G, M, N, O>
{
    fn calculate_viscous_dissipation_linear_element(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> Scalar
    {
        self.get_constitutive_model()
        .calculate_viscous_dissipation(
            &self.calculate_deformation_gradient(nodal_coordinates),
            &self.calculate_deformation_gradient_rate(nodal_coordinates, nodal_velocities)
        ) * self.get_integration_weight()
    }
    fn calculate_dissipation_potential_linear_element(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> Scalar
    {
        self.get_constitutive_model()
        .calculate_dissipation_potential(
            &self.calculate_deformation_gradient(nodal_coordinates),
            &self.calculate_deformation_gradient_rate(nodal_coordinates, nodal_velocities)
        ) * self.get_integration_weight()
    }
}

pub trait HyperviscoelasticLinearElement<'a, C, const G: usize, const M: usize, const N: usize, const O: usize>
where
    C: Hyperviscoelastic<'a>,
    Self: ElasticHyperviscousLinearElement<'a, C, G, M, N, O>
{
    fn calculate_helmholtz_free_energy_linear_element(&self, nodal_coordinates: &NodalCoordinates<N>) -> Scalar
    {
        self.get_constitutive_model()
        .calculate_helmholtz_free_energy_density(
            &self.calculate_deformation_gradient(nodal_coordinates)
        ) * self.get_integration_weight()
    }
}

macro_rules! linear_element_boilerplate
{
    ($element: ident) =>
    {
        impl<'a, C> ElasticFiniteElement<'a, C, G, N> for $element<C>
        where
            C: Elastic<'a>
        {
            fn calculate_nodal_forces(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalForces<N>
            {
                self.calculate_nodal_forces_linear_element(nodal_coordinates)
            }
            fn calculate_nodal_stiffnesses(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalStiffnesses<N>
            {
                self.calculate_nodal_stiffnesses_linear_element(nodal_coordinates)
            }
        }
        impl<'a, C> ElasticLinearElement<'a, C, G, M, N, O> for $element<C>
        where
            C: Elastic<'a>
        {}
        impl<'a, C> HyperelasticFiniteElement<'a, C, G, N> for $element<C>
        where
            C: Hyperelastic<'a>
        {
            fn calculate_helmholtz_free_energy(&self, nodal_coordinates: &NodalCoordinates<N>) -> Scalar
            {
                self.calculate_helmholtz_free_energy_linear_element(nodal_coordinates)
            }
        }
        impl<'a, C> HyperelasticLinearElement<'a, C, G, M, N, O> for $element<C>
        where
            C: Hyperelastic<'a>
        {}
        impl<'a, C> ViscoelasticFiniteElement<'a, C, G, N> for $element<C>
        where
            C: Viscoelastic<'a>
        {
            fn calculate_nodal_forces(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> NodalForces<N>
            {
                self.calculate_nodal_forces_linear_element(nodal_coordinates, nodal_velocities)
            }
            fn calculate_nodal_stiffnesses(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> NodalStiffnesses<N>
            {
                self.calculate_nodal_stiffnesses_linear_element(nodal_coordinates, nodal_velocities)
            }
        }
        impl<'a, C> ViscoelasticLinearElement<'a, C, G, M, N, O> for $element<C>
        where
            C: Viscoelastic<'a>
        {}
        impl<'a, C> ElasticHyperviscousFiniteElement<'a, C, G, N> for $element<C>
        where
            C: ElasticHyperviscous<'a>
        {
            fn calculate_viscous_dissipation(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> Scalar
            {
                self.calculate_viscous_dissipation_linear_element(nodal_coordinates, nodal_velocities)
            }
            fn calculate_dissipation_potential(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> Scalar
            {
                self.calculate_dissipation_potential_linear_element(nodal_coordinates, nodal_velocities)
            }
        }
        impl<'a, C> ElasticHyperviscousLinearElement<'a, C, G, M, N, O> for $element<C>
        where
            C: ElasticHyperviscous<'a>
        {}
        impl<'a, C> HyperviscoelasticFiniteElement<'a, C, G, N> for $element<C>
        where
            C: Hyperviscoelastic<'a>
        {
            fn calculate_helmholtz_free_energy(&self, nodal_coordinates: &NodalCoordinates<N>) -> Scalar
            {
                self.calculate_helmholtz_free_energy_linear_element(nodal_coordinates)
            }
        }
        impl<'a, C> HyperviscoelasticLinearElement<'a, C, G, M, N, O> for $element<C>
        where
            C: Hyperviscoelastic<'a>
        {}
    }
}
pub(crate) use linear_element_boilerplate;