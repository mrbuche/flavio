#[cfg(test)]
mod test;

pub mod triangle;

use super::*;

pub trait CompositeSurfaceElement<'a, C, const G: usize, const M: usize, const N: usize, const O: usize, const P: usize, const Q: usize>
where
    C: Constitutive<'a>,
    Self: CompositeElement<'a, C, G, M, N, O, P, Q>
{
    fn calculate_bases<const I: usize>(nodal_coordinates: &Coordinates<I, O>) -> Bases<I, P>
    {
        Self::calculate_standard_gradient_operators().iter()
        .map(|calculate_standard_gradient_operator|
            calculate_standard_gradient_operator.iter()
            .zip(nodal_coordinates.iter())
            .map(|(standard_gradient_operator_a, nodal_coordinates_a)|
                standard_gradient_operator_a.iter()
                .map(|standard_gradient_operator_a_m|
                    nodal_coordinates_a.iter()
                    .map(|nodal_coordinates_a_i|
                        nodal_coordinates_a_i * standard_gradient_operator_a_m
                    ).collect()
                ).collect()
            ).sum()
        ).collect()
    }
    fn calculate_deformation_gradients_composite_surface_element(&self, nodal_coordinates: &NodalCoordinates<O>) -> DeformationGradients<G>
    {
        self.get_projected_gradient_vectors().iter()
        .zip(Self::calculate_normals(nodal_coordinates).iter()
        .zip(self.get_reference_normals().iter()))
        .map(|(projected_gradient_vectors, (normal, reference_normal))|
            nodal_coordinates.iter()
            .zip(projected_gradient_vectors.iter())
            .map(|(nodal_coordinate, projected_gradient_vector)|
                DeformationGradient::dyad(nodal_coordinate, projected_gradient_vector)
            ).sum::<DeformationGradient>() + DeformationGradient::dyad(
                normal, reference_normal
            )
        ).collect()
    }
    fn calculate_deformation_gradient_rates_composite_surface_element(&self, nodal_coordinates: &NodalCoordinates<O>, nodal_velocities: &NodalVelocities<O>) -> DeformationGradientRates<G>
    {
        self.get_projected_gradient_vectors().iter()
        .zip(Self::calculate_normal_rates(nodal_coordinates, nodal_velocities).iter()
        .zip(self.get_reference_normals().iter()))
        .map(|(projected_gradient_vectors, (normal_rate, reference_normal))|
            nodal_velocities.iter()
            .zip(projected_gradient_vectors.iter())
            .map(|(nodal_velocity, projected_gradient_vector)|
                DeformationGradientRate::dyad(nodal_velocity, projected_gradient_vector)
            ).sum::<DeformationGradientRate>() + DeformationGradientRate::dyad(
                normal_rate, reference_normal
            )
        ).collect()
    }
    fn calculate_dual_bases<const I: usize>(nodal_coordinates: &Coordinates<I, O>) -> Bases<I, P>
    {
        let bases_vectors = Self::calculate_bases(nodal_coordinates);
        bases_vectors.iter()
        .map(|basis_vectors|
            basis_vectors.iter()
            .map(|basis_vectors_m|
                basis_vectors.iter()
                .map(|basis_vectors_n|
                    basis_vectors_m * basis_vectors_n
                ).collect()
            ).collect::<TensorRank2<M, I, I>>()
            .inverse()
            .iter()
            .map(|metric_tensor_m|
                metric_tensor_m.iter()
                .zip(basis_vectors.iter())
                .map(|(metric_tensor_mn, basis_vectors_n)|
                    basis_vectors_n * metric_tensor_mn
                ).sum()
            ).collect()
        ).collect()
    }
    fn calculate_normals<const I: usize>(nodal_coordinates: &Coordinates<I, O>) -> Normals<I, P>
    {
        Self::calculate_bases(nodal_coordinates).iter()
        .map(|basis_vectors|
            basis_vectors[0].cross(&basis_vectors[1]).normalized()
        ).collect()
    }
    fn calculate_normal_gradients(nodal_coordinates: &Coordinates<1, O>) -> NormalGradientss<G, O>
    {
        todo!()
    }
    fn calculate_normal_rates(nodal_coordinates: &NodalCoordinates<O>, nodal_velocities: &NodalVelocities<O>) -> NormalRates<G>
    {
        todo!()
    }
    fn calculate_normal_tangents(nodal_coordinates: &Coordinates<1, O>) -> NormalTangentss<G, O>
    {
        todo!()
    }
    fn calculate_projected_gradient_vectors_composite_surface_element(reference_nodal_coordinates: &ReferenceNodalCoordinates<O>) -> ProjectedGradientVectors<G, N>
    {
        todo!()
    }
    fn get_reference_normals(&self) -> &ReferenceNormals<P>;
}

macro_rules! composite_surface_element_boilerplate
{
    ($element: ident) =>
    {
        impl<'a, C> CompositeSurfaceElement<'a, C, G, M, N, O, P, Q> for $element<C>
        where
            C: Constitutive<'a>
        {
            fn get_reference_normals(&self) -> &ReferenceNormals<P>
            {
                &self.reference_normals
            }
        }
    }
}
pub(crate) use composite_surface_element_boilerplate;