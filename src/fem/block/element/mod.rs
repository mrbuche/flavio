pub mod linear;

use super::*;

type CurrentNodalCoordinates<const N: usize> = CurrentCoordinates<N>;
type HelmholtzFreeEnergies<const G: usize> = [Scalar; G];
type IntegrationWeights<const G: usize> = [Scalar; G];
type NodalForces<const N: usize> = Forces<N>;
type NodalStiffnesses<const N: usize> = Stiffnesses<N>;
type ReferenceNodalCoordinates<const N: usize> = ReferenceCoordinates<N>;
type TangentStiffnesses<const G: usize> = FirstPiolaKirchoffTangentStiffnesses<G>;

pub trait FiniteElementTraits<'a, C, const G: usize, const N: usize>
where
    C: ConstitutiveModel<'a>
{
    fn calculate_deformation_gradients(&self, _current_nodal_coordinates: &CurrentNodalCoordinates<N>) -> DeformationGradients<G>;
    fn calculate_first_piola_kirchoff_stresses(&self, deformation_gradients: &DeformationGradients<G>) -> FirstPiolaKirchoffStresses<G>
    {
        todo!();
    }
    fn calculate_helmholtz_free_energy(&self, current_nodal_coordinates: &CurrentNodalCoordinates<N>) -> Scalar
    {
        todo!();
    }
    fn calculate_helmholtz_free_energy_densities(&self, deformation_gradients: &DeformationGradients<G>) -> HelmholtzFreeEnergies<G>
    {
        todo!();
    }
    fn calculate_nodal_forces(&self, current_nodal_coordinates: &CurrentNodalCoordinates<N>) -> NodalForces<N>;
    fn calculate_nodal_stiffnesses(&self, current_nodal_coordinates: &CurrentNodalCoordinates<N>) -> NodalStiffnesses<N>;
    fn calculate_tangent_stiffnesses(&self, deformation_gradients: &DeformationGradients<G>) -> TangentStiffnesses<G>
    {
        todo!()
    }
    fn get_constitutive_model(&self, integration_point: usize) -> &C;
    fn get_integration_weights(&self) -> IntegrationWeights<G>;
    fn new(constitutive_model_parameters: ConstitutiveModelParameters<'a>, reference_nodal_coordinates: ReferenceNodalCoordinates<N>) -> Self;
}