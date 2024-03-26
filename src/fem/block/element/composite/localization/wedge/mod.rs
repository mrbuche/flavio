#[cfg(test)]
mod test;

use super::*;

const G: usize = 3;
const M: usize = 2;
const N: usize = 12;
const O: usize = 6;
const P: usize = 4;
const Q: usize = 3;

const INTEGRATION_WEIGHT: Scalar = 1.0/6.0;

pub struct Wedge<C>
{
    constitutive_models: [C; G],
    projected_gradient_vectors: ProjectedGradientVectors<G, N>,
    scaled_composite_jacobians: Scalars<G>,
    scaled_reference_normals: ScaledReferenceNormals<G, P>
}

impl<'a, C> FiniteElement<'a, C, G, N> for Wedge<C>
where
    C: Constitutive<'a>
{
    fn new(constitutive_model_parameters: Parameters<'a>, reference_nodal_coordinates: ReferenceNodalCoordinates<N>) -> Self
    {
        let nodal_coordinates_midplane = Self::calculate_midplane(&reference_nodal_coordinates);
        Self
        {
            constitutive_models: std::array::from_fn(|_| <C>::new(constitutive_model_parameters)),
            projected_gradient_vectors: Self::calculate_projected_gradient_vectors(&reference_nodal_coordinates),
            scaled_composite_jacobians: Self::calculate_scaled_composite_jacobian_at_integration_points(&nodal_coordinates_midplane),
            scaled_reference_normals: Self::calculate_scaled_reference_normals(&nodal_coordinates_midplane)
        }
    }
}

impl<'a, C> CompositeElement<'a, C, G, M, N, O, P, Q> for Wedge<C>
where
    C: Constitutive<'a>
{
    fn calculate_deformation_gradients(&self, nodal_coordinates: &NodalCoordinates<N>) -> DeformationGradients<G>
    {
        self.calculate_deformation_gradients_composite_localization_element(nodal_coordinates)
    }
    fn calculate_deformation_gradient_rates(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> DeformationGradientRates<G>
    {
        self.calculate_deformation_gradient_rates_composite_localization_element(nodal_coordinates, nodal_velocities)
    }
    fn calculate_projected_gradient_vectors(reference_nodal_coordinates: &ReferenceNodalCoordinates<N>) -> ProjectedGradientVectors<G, N>
    {
        todo!()
    }
    composite_surface_element_boilerplate_inner!{}
}

impl<'a, C> CompositeLocalizationElement<'a, C, G, M, N, O, P, Q> for Wedge<C>
where
    C: Constitutive<'a>
{
    fn calculate_midplane<const I: usize>(nodal_coordinates: &Coordinates<I, N>) -> Coordinates<I, O>
    {
        nodal_coordinates.iter().skip(3)
        .chain(nodal_coordinates.iter().skip(6).take(3))
        .zip(nodal_coordinates.iter().take(3)
        .chain(nodal_coordinates.iter().skip(3).take(3)))
        .map(|(nodal_coordinates_top, nodal_coordinates_bottom)|
            nodal_coordinates_top.iter()
            .zip(nodal_coordinates_bottom.iter())
            .map(|(nodal_coordinates_top_i, nodal_coordinates_bottom_i)|
                (nodal_coordinates_top_i + nodal_coordinates_bottom_i) * 0.5
            ).collect()
        ).collect()
    }
}
impl<'a, C> ElasticFiniteElement<'a, C, G, N> for Wedge<C>
where
    C: Elastic<'a>
{
    fn calculate_nodal_forces(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalForces<N>
    {
        todo!("Factors of 1/2 and stuff that are not present in surface implementation.")
    }
    fn calculate_nodal_stiffnesses(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalStiffnesses<N>
    {
        todo!()
    }
}
impl<'a, C> ViscoelasticFiniteElement<'a, C, G, N> for Wedge<C>
where
    C: Viscoelastic<'a>
{
    fn calculate_nodal_forces(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> NodalForces<N>
    {
        todo!()
    }
    fn calculate_nodal_stiffnesses(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> NodalStiffnesses<N>
    {
        todo!()
    }
}

composite_surface_or_localization_element_boilerplate!(Wedge);