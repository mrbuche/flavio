#[cfg(test)]
mod test;

pub mod tetrahedron;
pub mod triangle;

use crate::math::TensorRank1List;

use super::*;

type BasisVectors = Vectors<1, 3>;
type ReferenceBasisVectors = Vectors<0, 3>;
type StandardGradientOperator<const M: usize, const N: usize> = TensorRank1List<M, 9, N>;

pub trait LinearElement<'a, C, const G: usize, const M: usize, const N: usize>
where
    C: Constitutive<'a>,
    Self: FiniteElement<'a, C, G, N>
{
    fn calculate_deformation_gradient(&self, nodal_coordinates: &NodalCoordinates<N>) -> DeformationGradient
    {
        nodal_coordinates.iter()
        .zip(self.get_gradient_vectors().iter())
        .map(|(nodal_coordinate, gradient_vector)|
            DeformationGradient::dyad(nodal_coordinate, gradient_vector)
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
    fn calculate_gradient_vectors(reference_nodal_coordinates: &ReferenceNodalCoordinates<N>) -> GradientVectors<N>;
    fn calculate_standard_gradient_operator() -> StandardGradientOperator<M, N>;
    fn get_gradient_vectors(&self) -> &GradientVectors<N>;
}

pub trait LinearSurfaceElement<'a, C, const G: usize, const M: usize, const N: usize>
where
    C: Constitutive<'a>,
    Self: SurfaceElement<'a, C, G, N>
{
    fn calculate_basis_vectors(&self, nodal_coordinates: &NodalCoordinates<N>) -> BasisVectors
    {
        let basis_vectors_surface: Vectors<1, 2> =
            Self::calculate_standard_gradient_operator().iter()
            .zip(nodal_coordinates.iter())
            .map(|(standard_gradient_operator_a, nodal_coordinates_a)|
                standard_gradient_operator_a.iter()
                .map(|standard_gradient_operator_a_m|
                    nodal_coordinates_a.iter()
                    .map(|nodal_coordinates_a_i|
                        nodal_coordinates_a_i*standard_gradient_operator_a_m
                    ).collect()
                ).collect()
            ).sum();
        BasisVectors::new([
            basis_vectors_surface[0].as_array(),
            basis_vectors_surface[1].as_array(),
            basis_vectors_surface[0].cross(&basis_vectors_surface[1]).normalized().as_array()
        ])
    }
    fn calculate_deformation_gradient(&self, nodal_coordinates: &NodalCoordinates<N>) -> DeformationGradient
    {
        self.calculate_basis_vectors(nodal_coordinates).iter()
        .zip(self.get_reference_basis_vectors().iter())
        .map(|(basis_vector, reference_basis_vector)|
            DeformationGradient::dyad(basis_vector, reference_basis_vector)
        ).sum::<DeformationGradient>()
    }
    fn calculate_reference_basis_vectors(reference_nodal_coordinates: &ReferenceNodalCoordinates<N>) -> ReferenceBasisVectors
    {
        let reference_basis_vectors_surface: Vectors<0, 2> =
            Self::calculate_standard_gradient_operator().iter()
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
        ReferenceBasisVectors::new([
            reference_basis_vectors_surface[0].as_array(),
            reference_basis_vectors_surface[1].as_array(),
            reference_basis_vectors_surface[0].cross(&reference_basis_vectors_surface[1]).normalized().as_array()
        ])
    }
    fn calculate_standard_gradient_operator() -> StandardGradientOperator<M, N>;
    fn get_reference_basis_vectors(&self) -> &ReferenceBasisVectors;
    fn get_thickness(&self) -> &Scalar;
}

pub trait ElasticLinearElement<'a, C, const G: usize, const M: usize, const N: usize>
where
    C: Elastic<'a>,
    Self: LinearElement<'a, C, G, M, N>
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
                    gradient_vector_b, gradient_vector_a
                ).transpose()
            ).collect()
        ).collect()
    }
}

pub trait ElasticLinearSurfaceElement<'a, C, const G: usize, const M: usize, const N: usize>
where
    C: Elastic<'a>,
    Self: SurfaceElement<'a, C, G, N>
{
    fn calculate_nodal_forces_linear_surface_element(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalForces<N>
    {
        todo!()
    }
    fn calculate_nodal_stiffnesses_linear_surface_element(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalStiffnesses<N>
    {
        todo!()
    }
}

pub trait HyperelasticLinearElement<'a, C, const G: usize, const M: usize, const N: usize>
where
    C: Hyperelastic<'a>,
    Self: ElasticLinearElement<'a, C, G, M, N>
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

pub trait ViscoelasticLinearElement<'a, C, const G: usize, const M: usize, const N: usize>
where
    C: Viscoelastic<'a>,
    Self: LinearElement<'a, C, G, M, N>
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

pub trait ElasticHyperviscousLinearElement<'a, C, const G: usize, const M: usize, const N: usize>
where
    C: ElasticHyperviscous<'a>,
    Self: ViscoelasticLinearElement<'a, C, G, M, N>
{
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

pub trait HyperviscoelasticLinearElement<'a, C, const G: usize, const M: usize, const N: usize>
where
    C: Hyperviscoelastic<'a>,
    Self: ElasticHyperviscousLinearElement<'a, C, G, M, N>
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