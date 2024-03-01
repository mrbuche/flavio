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
        // nodal_coordinates.iter()
        // .zip(self.get_gradient_vectors().iter())
        // .map(|(nodal_coordinate, gradient_vector)|
        //     DeformationGradient::dyad(nodal_coordinate, gradient_vector)
        // ).sum::<DeformationGradient>() + DeformationGradient::dyad(
        //     &self.calculate_normal(nodal_coordinates), self.get_reference_normal()
        // )
        todo!()
    }
    fn calculate_deformation_gradient_rate_linear_localization_element(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> DeformationGradientRate
    {
        // nodal_velocities.iter()
        // .zip(self.get_gradient_vectors().iter())
        // .map(|(nodal_velocity, gradient_vector)|
        //     DeformationGradient::dyad(nodal_velocity, gradient_vector)
        // ).sum::<DeformationGradient>() + DeformationGradientRate::dyad(
        //     &self.calculate_normal_rate(nodal_coordinates, nodal_velocities), self.get_reference_normal()
        // )
        todo!()
    }
    fn calculate_nodal_jump(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalCoordinates<O>
    {
        todo!()
    }
    fn calculate_nodal_midplane(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalCoordinates<O>
    {
        todo!()
    }
}