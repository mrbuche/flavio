#[cfg(test)]
mod test;

pub mod composite;
pub mod linear;

use super::*;

pub trait FiniteElement<'a, C, const G: usize, const N: usize>
where
    C: Constitutive<'a>,
{
    fn new(
        constitutive_model_parameters: Parameters<'a>,
        reference_nodal_coordinates: ReferenceNodalCoordinates<N>,
    ) -> Self;
}

pub trait SurfaceElement<'a, C, const G: usize, const N: usize>
where
    C: Constitutive<'a>,
{
    fn new(
        constitutive_model_parameters: Parameters<'a>,
        reference_nodal_coordinates: ReferenceNodalCoordinates<N>,
        thickness: &Scalar,
    ) -> Self;
}

pub trait CohesiveElement<'a, C, const G: usize, const N: usize>
where
    C: Cohesive<'a>,
    Self: FiniteElement<'a, C, G, N>,
{
    fn calculate_nodal_forces(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalForces<N>;
    fn calculate_nodal_stiffnesses(
        &self,
        nodal_coordinates: &NodalCoordinates<N>,
    ) -> NodalStiffnesses<N>;
}

pub trait ElasticFiniteElement<'a, C, const G: usize, const N: usize>
where
    C: Elastic<'a>,
{
    fn calculate_nodal_forces(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalForces<N>;
    fn calculate_nodal_stiffnesses(
        &self,
        nodal_coordinates: &NodalCoordinates<N>,
    ) -> NodalStiffnesses<N>;
}

pub trait HyperelasticFiniteElement<'a, C, const G: usize, const N: usize>
where
    C: Hyperelastic<'a>,
    Self: ElasticFiniteElement<'a, C, G, N>,
{
    fn calculate_helmholtz_free_energy(&self, nodal_coordinates: &NodalCoordinates<N>) -> Scalar;
}

pub trait ViscoelasticFiniteElement<'a, C, const G: usize, const N: usize>
where
    C: Viscoelastic<'a>,
{
    fn calculate_nodal_forces(
        &self,
        nodal_coordinates: &NodalCoordinates<N>,
        nodal_velocities: &NodalVelocities<N>,
    ) -> NodalForces<N>;
    fn calculate_nodal_stiffnesses(
        &self,
        nodal_coordinates: &NodalCoordinates<N>,
        nodal_velocities: &NodalVelocities<N>,
    ) -> NodalStiffnesses<N>;
}

pub trait ElasticHyperviscousFiniteElement<'a, C, const G: usize, const N: usize>
where
    C: ElasticHyperviscous<'a>,
    Self: ViscoelasticFiniteElement<'a, C, G, N>,
{
    fn calculate_viscous_dissipation(
        &self,
        nodal_coordinates: &NodalCoordinates<N>,
        nodal_velocities: &NodalVelocities<N>,
    ) -> Scalar;
    fn calculate_dissipation_potential(
        &self,
        nodal_coordinates: &NodalCoordinates<N>,
        nodal_velocities: &NodalVelocities<N>,
    ) -> Scalar;
}

pub trait HyperviscoelasticFiniteElement<'a, C, const G: usize, const N: usize>
where
    C: Hyperviscoelastic<'a>,
    Self: ElasticHyperviscousFiniteElement<'a, C, G, N>,
{
    fn calculate_helmholtz_free_energy(&self, nodal_coordinates: &NodalCoordinates<N>) -> Scalar;
}
