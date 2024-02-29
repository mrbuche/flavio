#[cfg(test)]
mod test;

pub mod tetrahedron;
pub mod triangle;

use crate::math::TensorRank1List;

use super::*;

type Basis = Vectors<1, 2>;
type Normal = Vector<1>;
type NormalRate = Vector<1>;
type ReferenceBasis = Vectors<0, 2>;
type ReferenceNormal = Vector<0>;
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
    fn calculate_deformation_gradient_rate(&self, _: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> DeformationGradientRate
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
    Self: LinearElement<'a, C, G, M, N>
{
    fn calculate_basis_vectors(&self, nodal_coordinates: &NodalCoordinates<N>) -> Basis
    {
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
        ).sum()
    }
    fn calculate_deformation_gradient_linear_surface_element(&self, nodal_coordinates: &NodalCoordinates<N>) -> DeformationGradient
    {
        nodal_coordinates.iter()
        .zip(self.get_gradient_vectors().iter())
        .map(|(nodal_coordinate, gradient_vector)|
            DeformationGradient::dyad(nodal_coordinate, gradient_vector)
        ).sum::<DeformationGradient>() + DeformationGradient::dyad(
            &self.calculate_normal(nodal_coordinates), self.get_reference_normal()
        )
    }
    fn calculate_deformation_gradient_rate_linear_surface_element(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> DeformationGradientRate
    {
        nodal_velocities.iter()
        .zip(self.get_gradient_vectors().iter())
        .map(|(nodal_velocity, gradient_vector)|
            DeformationGradient::dyad(nodal_velocity, gradient_vector)
        ).sum::<DeformationGradient>() + DeformationGradient::dyad(
            &self.calculate_normal_rate(nodal_coordinates, nodal_velocities), self.get_reference_normal()
        )
    }
    fn calculate_gradient_vectors_linear_surface_element(reference_nodal_coordinates: &ReferenceNodalCoordinates<N>) -> GradientVectors<N>
    {
        Self::calculate_standard_gradient_operator().iter()
        .map(|standard_gradient_operator_a|
            standard_gradient_operator_a.iter()
            .zip(Self::calculate_reference_dual_basis_vectors(reference_nodal_coordinates).iter())
            .map(|(standard_gradient_operator_a_m, dual_reference_basis_vector_m)|
                dual_reference_basis_vector_m*standard_gradient_operator_a_m
            ).sum()
        ).collect()
    }
    fn calculate_normal(&self, nodal_coordinates: &NodalCoordinates<N>) -> Normal
    {
        let basis_vectors = self.calculate_basis_vectors(nodal_coordinates);
        basis_vectors[0].cross(&basis_vectors[1]).normalized()
    }
    fn calculate_normal_rate(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> NormalRate
    {
        todo!()
    }
    fn calculate_reference_basis_vectors(reference_nodal_coordinates: &ReferenceNodalCoordinates<N>) -> ReferenceBasis
    {
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
        ).sum()
    }
    fn calculate_reference_normal(reference_dual_basis_vectors: &ReferenceBasis) -> ReferenceNormal
    {
        reference_dual_basis_vectors[0].cross(&reference_dual_basis_vectors[1]).normalized()
    }
    fn calculate_reference_dual_basis_vectors(reference_nodal_coordinates: &ReferenceNodalCoordinates<N>) -> ReferenceBasis
    {
        let reference_basis_vectors_surface = Self::calculate_reference_basis_vectors(reference_nodal_coordinates);
        reference_basis_vectors_surface.iter()
        .map(|reference_basis_vectors_surface_m|
            reference_basis_vectors_surface.iter()
            .map(|reference_basis_vectors_surface_n|
                reference_basis_vectors_surface_m*reference_basis_vectors_surface_n
            ).collect()
        ).collect::<TensorRank2<2, 0, 0>>()
        .inverse()
        .iter()
        .map(|reference_metric_tensor_m|
            reference_metric_tensor_m.iter()
            .zip(reference_basis_vectors_surface.iter())
            .map(|(reference_metric_tensor_mn, reference_basis_vectors_surface_n)|
            reference_basis_vectors_surface_n*reference_metric_tensor_mn
            ).sum()
        ).collect()
    }
    fn get_reference_normal(&self) -> &ReferenceNormal;
}

pub trait ElasticLinearElement<'a, C, const G: usize, const M: usize, const N: usize>
where
    C: Elastic<'a>,
    Self: LinearElement<'a, C, G, M, N>
{
    fn calculate_nodal_forces_linear_element(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalForces<N>
    {
        let first_piola_kirchoff_stress = self.get_constitutive_models()[0]
        .calculate_first_piola_kirchoff_stress(
            &self.calculate_deformation_gradient(nodal_coordinates)
        );
        self.get_gradient_vectors().iter().map(|gradient_vector|
            &first_piola_kirchoff_stress * gradient_vector
        ).collect()
    }
    fn calculate_nodal_stiffnesses_linear_element(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalStiffnesses<N>
    {
        let first_piola_kirchoff_tangent_stiffness = self.get_constitutive_models()[0]
        .calculate_first_piola_kirchoff_tangent_stiffness(
            &self.calculate_deformation_gradient(nodal_coordinates)
        );
        self.get_gradient_vectors().iter()
        .map(|gradient_vector_a|
            self.get_gradient_vectors().iter()
            .map(|gradient_vector_b|
                first_piola_kirchoff_tangent_stiffness
                .contract_second_fourth_indices_with_first_indices_of(
                    gradient_vector_b, gradient_vector_a
                ).transpose()
            ).collect()
        ).collect()
    }
}

pub trait HyperelasticLinearElement<'a, C, const G: usize, const M: usize, const N: usize>
where
    C: Hyperelastic<'a>,
    Self: ElasticLinearElement<'a, C, G, M, N>
{
    fn calculate_helmholtz_free_energy_linear_element(&self, nodal_coordinates: &NodalCoordinates<N>) -> Scalar
    {
        self.get_constitutive_models()[0]
        .calculate_helmholtz_free_energy_density(
            &self.calculate_deformation_gradient(nodal_coordinates)
        )
    }
}

pub trait ViscoelasticLinearElement<'a, C, const G: usize, const M: usize, const N: usize>
where
    C: Viscoelastic<'a>,
    Self: LinearElement<'a, C, G, M, N>
{
    fn calculate_nodal_forces_linear_element(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> NodalForces<N>
    {
        let first_piola_kirchoff_stress = self.get_constitutive_models()[0]
        .calculate_first_piola_kirchoff_stress(
            &self.calculate_deformation_gradient(nodal_coordinates),
            &self.calculate_deformation_gradient_rate(nodal_coordinates, nodal_velocities)
        );
        self.get_gradient_vectors().iter().map(|gradient_vector|
            &first_piola_kirchoff_stress * gradient_vector
        ).collect()
    }
    fn calculate_nodal_stiffnesses_linear_element(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> NodalStiffnesses<N>
    {
        let first_piola_kirchoff_tangent_stiffness = self.get_constitutive_models()[0]
        .calculate_first_piola_kirchoff_rate_tangent_stiffness(
            &self.calculate_deformation_gradient(nodal_coordinates),
            &self.calculate_deformation_gradient_rate(nodal_coordinates, nodal_velocities)
        );
        self.get_gradient_vectors().iter()
        .map(|gradient_vector_a|
            self.get_gradient_vectors().iter()
            .map(|gradient_vector_b|
                first_piola_kirchoff_tangent_stiffness
                .contract_second_fourth_indices_with_first_indices_of(
                    gradient_vector_b, gradient_vector_a
                ).transpose()
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
        self.get_constitutive_models()[0]
        .calculate_viscous_dissipation(
            &self.calculate_deformation_gradient(nodal_coordinates),
            &self.calculate_deformation_gradient_rate(nodal_coordinates, nodal_velocities)
        )
    }
    fn calculate_dissipation_potential_linear_element(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> Scalar
    {
        self.get_constitutive_models()[0]
        .calculate_dissipation_potential(
            &self.calculate_deformation_gradient(nodal_coordinates),
            &self.calculate_deformation_gradient_rate(nodal_coordinates, nodal_velocities)
        )
    }
}

pub trait HyperviscoelasticLinearElement<'a, C, const G: usize, const M: usize, const N: usize>
where
    C: Hyperviscoelastic<'a>,
    Self: ElasticHyperviscousLinearElement<'a, C, G, M, N>
{
    fn calculate_helmholtz_free_energy_linear_element(&self, nodal_coordinates: &NodalCoordinates<N>) -> Scalar
    {
        self.get_constitutive_models()[0]
        .calculate_helmholtz_free_energy_density(
            &self.calculate_deformation_gradient(nodal_coordinates)
        )
    }
}