#[cfg(test)]
mod test;

pub mod tetrahedron;

use super::*;

pub trait LinearFiniteElement<'a, C, const G: usize, const N: usize>
where
    C: ConstitutiveModel<'a>,
    Self: FiniteElement<'a, C, G, N>
{
    fn calculate_deformation_gradient(&self, current_nodal_coordinates: &CurrentNodalCoordinates<N>) -> DeformationGradient
    {
        current_nodal_coordinates.iter()
        .zip(self.get_gradient_vectors().iter())
        .map(|(current_coordinates, gradient_vector)|
            DeformationGradient::dyad(current_coordinates, gradient_vector)
        ).sum()
    }
    fn calculate_gradient_vectors(reference_nodal_coordinates: &ReferenceNodalCoordinates<N>) -> GradientVectors<N>
    {
        let standard_gradient_operator = Self::calculate_standard_gradient_operator();
        (reference_nodal_coordinates * &standard_gradient_operator).inverse_transpose() * standard_gradient_operator
    }
    fn calculate_nodal_forces_linear_element(&self, current_nodal_coordinates: &CurrentNodalCoordinates<N>) -> NodalForces<N>
    {
        let deformation_gradient = self.calculate_deformation_gradient(current_nodal_coordinates);
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
    fn calculate_nodal_stiffnesses_linear_element(&self, current_nodal_coordinates: &CurrentNodalCoordinates<N>) -> NodalStiffnesses<N>
    {
        let deformation_gradient = self.calculate_deformation_gradient(current_nodal_coordinates);
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
    fn calculate_standard_gradient_operator() -> StandardGradientOperator<N>;
    fn get_gradient_vectors(&self) -> &GradientVectors<N>;
}

pub trait HyperelasticLinearFiniteElement<'a, C, const G: usize, const N: usize>
where
    C: ConstitutiveModel<'a> + HyperelasticConstitutiveModel,
    Self: FiniteElement<'a, C, G, N> + LinearFiniteElement<'a, C, G, N>
{
    fn calculate_helmholtz_free_energy_linear_element(&self, current_nodal_coordinates: &CurrentNodalCoordinates<N>) -> Scalar
    {
        let deformation_gradient = self.calculate_deformation_gradient(current_nodal_coordinates);
        self.get_constitutive_models().iter()
        .zip(self.get_integration_weights().iter())
        .map(|(constitutive_model, integration_weight)|
            constitutive_model.calculate_helmholtz_free_energy_density(
                &deformation_gradient
            ) * integration_weight
        ).sum()
    }
}