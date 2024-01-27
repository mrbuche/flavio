#[cfg(test)]
mod test;

pub mod linear;

use super::*;

type GradientVectors<const N: usize> = Vectors<0, N>;
type IntegrationWeights<const G: usize> = Scalars<G>;
type StandardGradientOperator<const N: usize> = Vectors<9, N>;

pub trait FiniteElement<'a, C, const G: usize, const N: usize, Y>
where
    C: ConstitutiveModel<'a, Y>
{
    fn calculate_nodal_forces(&self, current_nodal_coordinates: &CurrentNodalCoordinates<N>) -> NodalForces<N>;
    fn calculate_nodal_stiffnesses(&self, current_nodal_coordinates: &CurrentNodalCoordinates<N>) -> NodalStiffnesses<N>;
    fn get_constitutive_models(&self) -> &[C; G];
    fn get_integration_weights(&self) -> IntegrationWeights<G>;
    fn new(constitutive_model_parameters: ConstitutiveModelParameters<'a>, reference_nodal_coordinates: ReferenceNodalCoordinates<N>) -> Self;
}

pub trait HyperelasticFiniteElement<'a, C, const G: usize, const N: usize, Y>
where
    C: ConstitutiveModel<'a, Y> + HyperelasticConstitutiveModel<'a, Y>,
    Self: FiniteElement<'a, C, G, N, Y>
{
    fn calculate_helmholtz_free_energy(&self, current_nodal_coordinates: &CurrentNodalCoordinates<N>) -> Scalar;
}