pub mod linear;

use super::*;

type CurrentNodalCoordinates<const N: usize> = CurrentCoordinates<N>;
type GradientVectors<const N: usize> = Vectors<0, N>;
type IntegrationWeights<const G: usize> = Scalars<G>;
type NodalForces<const N: usize> = Forces<N>;
type NodalStiffnesses<const N: usize> = Stiffnesses<N>;
type ReferenceNodalCoordinates<const N: usize> = ReferenceCoordinates<N>;
type StandardGradientOperator<const N: usize> = Vectors<9, N>;

pub trait FiniteElement<'a, C, const G: usize, const N: usize>
where
    C: ConstitutiveModel<'a>
{
    fn calculate_deformation_gradients(&self, _current_nodal_coordinates: &CurrentNodalCoordinates<N>) -> DeformationGradients<G>;
    fn calculate_helmholtz_free_energy(&self, current_nodal_coordinates: &CurrentNodalCoordinates<N>) -> Scalar
    {
        self.calculate_deformation_gradients(current_nodal_coordinates).iter()
        .zip(self.get_constitutive_models().iter()
            .zip(self.get_integration_weights().iter())
        ).map(|(deformation_gradient, (constitutive_model, integration_weight))|
            constitutive_model.calculate_helmholtz_free_energy_density(
                deformation_gradient
            ) * integration_weight
        ).sum()
    }
    fn calculate_nodal_forces(&self, current_nodal_coordinates: &CurrentNodalCoordinates<N>) -> NodalForces<N>;
    fn calculate_nodal_stiffnesses(&self, current_nodal_coordinates: &CurrentNodalCoordinates<N>) -> NodalStiffnesses<N>;
    fn get_constitutive_models(&self) -> &[C; G];
    fn get_integration_weights(&self) -> IntegrationWeights<G>;
    fn new(constitutive_model_parameters: ConstitutiveModelParameters<'a>, reference_nodal_coordinates: ReferenceNodalCoordinates<N>) -> Self;
}