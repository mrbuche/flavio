#[cfg(test)]
mod test;

pub mod linear;

use super::*;

type GradientVectors<const N: usize> = Vectors<0, N>;
type IntegrationWeights<const G: usize> = Scalars<G>;
type StandardGradientOperator<const N: usize> = Vectors<9, N>;

pub trait FiniteElement<'a, C, const G: usize, const N: usize>
where
    C: Constitutive<'a>
{
    fn get_constitutive_models(&self) -> &[C; G];
    fn get_integration_weights(&self) -> IntegrationWeights<G>;
    fn new(constitutive_model_parameters: Parameters<'a>, reference_nodal_coordinates: ReferenceNodalCoordinates<N>) -> Self;
}

pub trait ElasticFiniteElement<'a, C, const G: usize, const N: usize>
where
    C: Elastic<'a>,
    Self: FiniteElement<'a, C, G, N>
{
    fn calculate_nodal_forces(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalForces<N>;
    fn calculate_nodal_stiffnesses(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalStiffnesses<N>;
}

pub trait HyperelasticFiniteElement<'a, C, const G: usize, const N: usize>
where
    C: Hyperelastic<'a>,
    Self: ElasticFiniteElement<'a, C, G, N>
{
    fn calculate_helmholtz_free_energy(&self, nodal_coordinates: &NodalCoordinates<N>) -> Scalar;
}

pub trait ViscoelasticFiniteElement<'a, C, const G: usize, const N: usize>
where
    C: Viscoelastic<'a>,
    Self: FiniteElement<'a, C, G, N>
{
    fn calculate_nodal_forces(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> NodalForces<N>;
    fn calculate_nodal_stiffnesses(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> NodalStiffnesses<N>;
}

pub trait HyperviscoelasticFiniteElement<'a, C, const G: usize, const N: usize>
where
    C: Hyperviscoelastic<'a>,
    Self: ViscoelasticFiniteElement<'a, C, G, N>
{
    fn calculate_helmholtz_free_energy(&self, nodal_coordinates: &NodalCoordinates<N>) -> Scalar;
    fn calculate_viscous_dissipation(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> Scalar;
    fn calculate_dissipation_potential(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> Scalar;
}