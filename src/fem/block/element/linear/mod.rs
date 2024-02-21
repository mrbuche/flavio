#[cfg(test)]
mod test;

pub mod tetrahedron;

use super::*;

pub trait LinearFiniteElement<'a, C, const G: usize, const N: usize>
where
    C: Constitutive<'a>,
    Self: FiniteElement<'a, C, G, N>
{
    fn calculate_deformation_gradient(&self, nodal_coordinates: &NodalCoordinates<N>) -> DeformationGradient
    {
        nodal_coordinates.iter()
        .zip(self.get_gradient_vectors().iter())
        .map(|(current_nodal_coordinate, gradient_vector)|
            DeformationGradient::dyad(current_nodal_coordinate, gradient_vector)
        ).sum()
    }
    fn calculate_deformation_gradient_rate(&self, nodal_velocities: &NodalVelocities<N>) -> DeformationGradientRate
    {
        nodal_velocities.iter()
        .zip(self.get_gradient_vectors().iter())
        .map(|(nodal_velocity, gradient_vector)|
            DeformationGradientRate::dyad(nodal_velocity, gradient_vector)
        ).sum()
    }
    fn calculate_gradient_vectors(reference_nodal_coordinates: &ReferenceNodalCoordinates<N>) -> GradientVectors<N>
    {
        let standard_gradient_operator = Self::calculate_standard_gradient_operator();
        (reference_nodal_coordinates * &standard_gradient_operator).inverse_transpose() * standard_gradient_operator
    }
    fn calculate_standard_gradient_operator() -> StandardGradientOperator<N>;
    fn get_gradient_vectors(&self) -> &GradientVectors<N>;
}

pub trait ElasticLinearFiniteElement<'a, C, const G: usize, const N: usize>
where
    C: Elastic<'a>,
    Self: LinearFiniteElement<'a, C, G, N>
{
    fn calculate_nodal_forces_linear_element(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalForces<N>
    {
        let deformation_gradient = self.calculate_deformation_gradient(nodal_coordinates);
        self.get_gradient_vectors().iter()
        .map(|gradient_vector|
            self.get_constitutive_models().iter()
            .zip(self.get_integration_weights().iter())
            .map(|(constitutive_model, integration_weight)|
                constitutive_model.calculate_first_piola_kirchoff_stress(
                    &deformation_gradient
                ) * integration_weight
            ).sum::<FirstPiolaKirchoffStress>() * gradient_vector
        ).collect()
    }
    fn calculate_nodal_stiffnesses_linear_element(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalStiffnesses<N>
    {
        let deformation_gradient = self.calculate_deformation_gradient(nodal_coordinates);
        self.get_gradient_vectors().iter()
        .map(|gradient_vector_a|
            self.get_gradient_vectors().iter()
            .map(|gradient_vector_b|
                self.get_constitutive_models().iter()
                .zip(self.get_integration_weights().iter())
                .map(|(constitutive_model, integration_weight)|
                    constitutive_model.calculate_first_piola_kirchoff_tangent_stiffness(
                        &deformation_gradient
                    ) * integration_weight
                ).sum::<FirstPiolaKirchoffTangentStiffness>()
                .contract_second_fourth_indices_with_first_indices_of(
                    gradient_vector_a, gradient_vector_b
                )
            ).collect()
        ).collect()
    }
}

pub trait HyperelasticLinearFiniteElement<'a, C, const G: usize, const N: usize>
where
    C: Hyperelastic<'a>,
    Self: ElasticLinearFiniteElement<'a, C, G, N>
{
    fn calculate_helmholtz_free_energy_linear_element(&self, nodal_coordinates: &NodalCoordinates<N>) -> Scalar
    {
        let deformation_gradient = self.calculate_deformation_gradient(nodal_coordinates);
        self.get_constitutive_models().iter()
        .zip(self.get_integration_weights().iter())
        .map(|(constitutive_model, integration_weight)|
            constitutive_model.calculate_helmholtz_free_energy_density(
                &deformation_gradient
            ) * integration_weight
        ).sum()
    }
}

pub trait ViscoelasticLinearFiniteElement<'a, C, const G: usize, const N: usize>
where
    C: Viscoelastic<'a>,
    Self: LinearFiniteElement<'a, C, G, N>
{
    fn calculate_nodal_forces_linear_element(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> NodalForces<N>
    {
        let deformation_gradient = self.calculate_deformation_gradient(nodal_coordinates);
        let deformation_gradient_rate = self.calculate_deformation_gradient_rate(nodal_velocities);
        self.get_gradient_vectors().iter()
        .map(|gradient_vector|
            self.get_constitutive_models().iter()
            .zip(self.get_integration_weights().iter())
            .map(|(constitutive_model, integration_weight)|
                constitutive_model.calculate_first_piola_kirchoff_stress(
                    &deformation_gradient, &deformation_gradient_rate
                ) * integration_weight
            ).sum::<FirstPiolaKirchoffStress>() * gradient_vector
        ).collect()
    }
    fn calculate_nodal_stiffnesses_linear_element(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> NodalStiffnesses<N>
    {
        let deformation_gradient = self.calculate_deformation_gradient(nodal_coordinates);
        let deformation_gradient_rate = self.calculate_deformation_gradient_rate(nodal_velocities);
        self.get_gradient_vectors().iter()
        .map(|gradient_vector_a|
            self.get_gradient_vectors().iter()
            .map(|gradient_vector_b|
                self.get_constitutive_models().iter()
                .zip(self.get_integration_weights().iter())
                .map(|(constitutive_model, integration_weight)|
                    constitutive_model.calculate_first_piola_kirchoff_rate_tangent_stiffness(
                        &deformation_gradient, &deformation_gradient_rate
                    ) * integration_weight
                ).sum::<FirstPiolaKirchoffTangentStiffness>()
                .contract_second_fourth_indices_with_first_indices_of(
                    gradient_vector_a, gradient_vector_b
                )
            ).collect()
        ).collect()
    }
}

pub trait HyperviscoelasticLinearFiniteElement<'a, C, const G: usize, const N: usize>
where
    C: Hyperviscoelastic<'a>,
    Self: ViscoelasticLinearFiniteElement<'a, C, G, N>
{
    fn calculate_helmholtz_free_energy_linear_element(&self, nodal_coordinates: &NodalCoordinates<N>) -> Scalar
    {
        let deformation_gradient = self.calculate_deformation_gradient(nodal_coordinates);
        self.get_constitutive_models().iter()
        .zip(self.get_integration_weights().iter())
        .map(|(constitutive_model, integration_weight)|
            constitutive_model.calculate_helmholtz_free_energy_density(
                &deformation_gradient
            ) * integration_weight
        ).sum()
    }
    fn calculate_viscous_dissipation_linear_element(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> Scalar
    {
        let deformation_gradient = self.calculate_deformation_gradient(nodal_coordinates);
        let deformation_gradient_rate = self.calculate_deformation_gradient_rate(nodal_velocities);
        self.get_constitutive_models().iter()
        .zip(self.get_integration_weights().iter())
        .map(|(constitutive_model, integration_weight)|
            constitutive_model.calculate_viscous_dissipation(
                &deformation_gradient, &deformation_gradient_rate
            ) * integration_weight
        ).sum()
    }
    fn calculate_dissipation_potential_linear_element(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> Scalar
    {
        let deformation_gradient = self.calculate_deformation_gradient(nodal_coordinates);
        let deformation_gradient_rate = self.calculate_deformation_gradient_rate(nodal_velocities);
        self.get_constitutive_models().iter()
        .zip(self.get_integration_weights().iter())
        .map(|(constitutive_model, integration_weight)|
            constitutive_model.calculate_dissipation_potential(
                &deformation_gradient, &deformation_gradient_rate
            ) * integration_weight
        ).sum()
    }
}