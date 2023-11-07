pub mod linear;

use super::*;

type CurrentNodalCoordinates<const N: usize> = CurrentCoordinates<N>;
type GradientVectors<const N: usize> = Vectors<0, N>;
type HelmholtzFreeEnergies<const G: usize> = Scalars<G>;
type IntegrationWeights<const G: usize> = Scalars<G>;
type NodalForces<const N: usize> = Forces<N>;
type NodalStiffnesses<const N: usize> = Stiffnesses<N>;
type ReferenceNodalCoordinates<const N: usize> = ReferenceCoordinates<N>;
type StandardGradientOperator<const N: usize> = Vectors<9, N>;
type TangentStiffnesses<const G: usize> = FirstPiolaKirchoffTangentStiffnesses<G>;

pub trait FiniteElement<'a, C, const G: usize, const N: usize>
where
    C: ConstitutiveModel<'a>
{
    fn calculate_deformation_gradients(&self, _current_nodal_coordinates: &CurrentNodalCoordinates<N>) -> DeformationGradients<G>;
    fn calculate_helmholtz_free_energy(&self, current_nodal_coordinates: &CurrentNodalCoordinates<N>) -> Scalar
    {
        self.calculate_helmholtz_free_energy_densities(
            &self.calculate_deformation_gradients(current_nodal_coordinates)
        ) * self.get_integration_weights()
    }
    fn calculate_helmholtz_free_energy_densities(&self, deformation_gradients: &DeformationGradients<G>) -> HelmholtzFreeEnergies<G>
    {
        self.get_constitutive_models().iter()
        .enumerate()
        .map(|(integration_point, constitutive_model)|
            constitutive_model
            .calculate_helmholtz_free_energy_density(
                &deformation_gradients[integration_point]
            )
        ).collect()
    }
    fn calculate_nodal_forces(&self, current_nodal_coordinates: &CurrentNodalCoordinates<N>) -> NodalForces<N>;
    fn calculate_nodal_stiffnesses(&self, current_nodal_coordinates: &CurrentNodalCoordinates<N>) -> NodalStiffnesses<N>;
    fn calculate_tangent_stiffnesses(&self, deformation_gradients: &DeformationGradients<G>) -> TangentStiffnesses<G>
    {
        self.get_constitutive_models().iter()
        .enumerate()
        .map(|(integration_point, constitutive_model)|
            constitutive_model
            .calculate_first_piola_kirchoff_tangent_stiffness(
                &deformation_gradients[integration_point]
            )
        ).collect()
    }
    fn get_constitutive_models(&self) -> &[C; G];
    fn get_integration_weights(&self) -> IntegrationWeights<G>;
    fn new(constitutive_model_parameters: ConstitutiveModelParameters<'a>, reference_nodal_coordinates: ReferenceNodalCoordinates<N>) -> Self;
}