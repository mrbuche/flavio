#[cfg(test)]
mod test;

pub mod triangle;

use super::*;

pub trait LinearSurfaceElement<'a, C, const G: usize, const M: usize, const N: usize, const O: usize>
where
    C: Constitutive<'a>,
    Self: LinearElement<'a, C, G, M, N, O>
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
        ).sum::<DeformationGradient>() + DeformationGradientRate::dyad(
            &self.calculate_normal_rate(nodal_coordinates, nodal_velocities), self.get_reference_normal()
        )
    }
    fn calculate_gradient_vectors_linear_surface_element(reference_nodal_coordinates: &ReferenceNodalCoordinates<N>) -> GradientVectors<O>
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
        let basis_vectors = self.calculate_basis_vectors(nodal_coordinates);
        let levi_civita_symbol = levi_civita::<1, 1, 1>();
        let normalization = basis_vectors[0].cross(&basis_vectors[1]).norm();
        let normal_vector = basis_vectors[0].cross(&basis_vectors[1])/normalization;
        let standard_gradient_operator = Self::calculate_standard_gradient_operator();
        TensorRank2::<3, 1, 1>::identity().iter()
        .zip(normal_vector.iter())
        .map(|(identity_i, normal_vector_i)|
            nodal_velocities.iter()
            .zip(standard_gradient_operator.iter())
            .map(|(nodal_velocity_a, standard_gradient_operator_a)|
                levi_civita_symbol.iter()
                .zip(nodal_velocity_a.iter())
                .map(|(levi_civita_symbol_m, nodal_velocity_a_m)|
                    levi_civita_symbol_m.iter()
                    .zip(basis_vectors[0].iter().zip(basis_vectors[1].iter()))
                    .map(|(levi_civita_symbol_mn, (basis_vector_0_n, basis_vector_1_n))|
                        levi_civita_symbol_mn.iter()
                        .zip(identity_i.iter().zip(normal_vector.iter()))
                        .map(|(levi_civita_symbol_mno, (identity_io, normal_vector_o))|
                            levi_civita_symbol_mno*(identity_io - normal_vector_i*normal_vector_o)
                        ).sum::<Scalar>()*(standard_gradient_operator_a[0]*basis_vector_1_n
                                         - standard_gradient_operator_a[1]*basis_vector_0_n)
                    ).sum::<Scalar>()*nodal_velocity_a_m
                ).sum::<Scalar>()
            ).sum()
        ).collect::<NormalRate>()/normalization
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
