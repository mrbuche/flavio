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
        let nodal_coordinates_midplane = Self::calculate_midplane_coordinates(nodal_coordinates);
        nodal_coordinates.iter()
        .zip(self.get_gradient_vectors().iter())
        .map(|(nodal_coordinate, gradient_vector)|
            DeformationGradient::dyad(nodal_coordinate, gradient_vector)
        ).sum::<DeformationGradient>() + DeformationGradient::dyad(
            &Self::calculate_normal(&nodal_coordinates_midplane), self.get_reference_normal()
        )
    }
    fn calculate_deformation_gradient_rate_linear_localization_element(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> DeformationGradientRate
    {
        todo!()
    }
    fn calculate_jump(nodal_coordinates: &NodalCoordinates<N>) -> Jump;
    fn calculate_midplane_coordinates<const I: usize>(nodal_coordinates: &Coordinates<I, N>) -> Coordinates<I, O>;
}