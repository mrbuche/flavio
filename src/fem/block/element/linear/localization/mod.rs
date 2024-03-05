#[cfg(test)]
mod test;

pub mod wedge;

use super::
{
    *, surface::*
};

pub trait LinearLocalizationElement<'a, C, const G: usize, const M: usize, const N: usize, const O: usize>
where
    C: Constitutive<'a>,
    Self: LinearSurfaceElement<'a, C, G, M, N, O>
{
    fn calculate_deformation_gradient_linear_localization_element(&self, nodal_coordinates: &NodalCoordinates<N>) -> DeformationGradient
    {
        let nodal_coordinates_midplane = Self::calculate_midplane(nodal_coordinates);
        nodal_coordinates.iter()
        .zip(self.get_gradient_vectors().iter())
        .map(|(nodal_coordinate, gradient_vector)|
            DeformationGradient::dyad(nodal_coordinate, gradient_vector)
        ).sum::<DeformationGradient>() + DeformationGradient::dyad(
            &Self::calculate_normal(&nodal_coordinates_midplane),
            self.get_reference_normal()
        )
    }
    fn calculate_deformation_gradient_rate_linear_localization_element(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> DeformationGradientRate
    {
        let nodal_coordinates_midplane = Self::calculate_midplane(nodal_coordinates);
        let nodal_velocities_midplane = Self::calculate_midplane(nodal_velocities);
        nodal_velocities.iter()
        .zip(self.get_gradient_vectors().iter())
        .map(|(nodal_velocity, gradient_vector)|
            DeformationGradientRate::dyad(nodal_velocity, gradient_vector)
        ).sum::<DeformationGradientRate>() + DeformationGradientRate::dyad(
            &Self::calculate_normal_rate(&nodal_coordinates_midplane, &nodal_velocities_midplane),
            self.get_reference_normal()
        )
    }
    fn calculate_midplane<const I: usize>(nodal_coordinates: &Coordinates<I, N>) -> Coordinates<I, O>;
}

macro_rules! linear_localization_element_boilerplate
{
    ($element: ident) =>
    {
        impl<'a, C> ElasticLinearElement<'a, C, G, M, N, O> for Wedge<C>
        where
            C: Elastic<'a>
        {}
        
        impl<'a, C> HyperelasticFiniteElement<'a, C, G, N> for Wedge<C>
        where
            C: Hyperelastic<'a>
        {
            fn calculate_helmholtz_free_energy(&self, nodal_coordinates: &NodalCoordinates<N>) -> Scalar
            {
                self.calculate_helmholtz_free_energy_linear_element(nodal_coordinates)
            }
        }
        impl<'a, C> HyperelasticLinearElement<'a, C, G, M, N, O> for Wedge<C>
        where
            C: Hyperelastic<'a>
        {}
        impl<'a, C> ViscoelasticLinearElement<'a, C, G, M, N, O> for Wedge<C>
        where
            C: Viscoelastic<'a>
        {}
        impl<'a, C> ElasticHyperviscousFiniteElement<'a, C, G, N> for Wedge<C>
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
        impl<'a, C> ElasticHyperviscousLinearElement<'a, C, G, M, N, O> for Wedge<C>
        where
            C: ElasticHyperviscous<'a>
        {}
        
        impl<'a, C> HyperviscoelasticFiniteElement<'a, C, G, N> for Wedge<C>
        where
            C: Hyperviscoelastic<'a>
        {
            fn calculate_helmholtz_free_energy(&self, nodal_coordinates: &NodalCoordinates<N>) -> Scalar
            {
                self.calculate_helmholtz_free_energy_linear_element(nodal_coordinates)
            }
        }
        impl<'a, C> HyperviscoelasticLinearElement<'a, C, G, M, N, O> for Wedge<C>
        where
            C: Hyperviscoelastic<'a>
        {}
    }
}
pub(crate) use linear_localization_element_boilerplate;