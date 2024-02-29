#[cfg(test)]
mod test;

use super::*;

const G: usize = 1;
const M: usize = 2;
const N: usize = 3;

pub struct Triangle<C>
{
    constitutive_models: [C; G],
    gradient_vectors: GradientVectors<N>,
    reference_normal: ReferenceNormal
}

impl<'a, C> FiniteElement<'a, C, G, N> for Triangle<C>
where
    C: Constitutive<'a>
{
    fn get_constitutive_models(&self) -> &[C; G]
    {
        &self.constitutive_models
    }
    fn get_integration_weights(&self) -> IntegrationWeights<G>
    {
        IntegrationWeights::new([1.0; G])
    }
    fn new(constitutive_model_parameters: Parameters<'a>, reference_nodal_coordinates: ReferenceNodalCoordinates<N>) -> Self
    {
        Self
        {
            constitutive_models: std::array::from_fn(|_| <C>::new(constitutive_model_parameters)),
            gradient_vectors: Self::calculate_gradient_vectors(&reference_nodal_coordinates),
            reference_normal: Self::calculate_reference_normal(&Self::calculate_reference_dual_basis_vectors(&reference_nodal_coordinates))
        }
    }
}

impl<'a, C> LinearElement<'a, C, G, M, N> for Triangle<C>
where
    C: Constitutive<'a>
{
    fn calculate_deformation_gradient(&self, nodal_coordinates: &NodalCoordinates<N>) -> DeformationGradient
    {
        self.calculate_deformation_gradient_linear_surface_element(nodal_coordinates)
    }
    fn calculate_deformation_gradient_rate(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> DeformationGradientRate
    {
        self.calculate_deformation_gradient_rate_linear_surface_element(nodal_coordinates, nodal_velocities)
    }
    fn calculate_gradient_vectors(reference_nodal_coordinates: &ReferenceNodalCoordinates<N>) -> GradientVectors<N>
    {
        Self::calculate_gradient_vectors_linear_surface_element(reference_nodal_coordinates)
    }
    fn calculate_standard_gradient_operator() -> StandardGradientOperator<M, N>
    {
        StandardGradientOperator::new([
            [-1.0, -1.0],
            [ 1.0,  0.0],
            [ 0.0,  1.0]
        ])
    }
    fn get_gradient_vectors(&self) -> &GradientVectors<N>
    {
        &self.gradient_vectors
    }
}

impl<'a, C> LinearSurfaceElement<'a, C, G, M, N> for Triangle<C>
where
    C: Constitutive<'a>
{
    fn get_reference_normal(&self) -> &ReferenceNormal
    {
        &self.reference_normal
    }
}

impl<'a, C> ElasticFiniteElement<'a, C, G, N> for Triangle<C>
where
    C: Elastic<'a>
{
    fn calculate_nodal_forces(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalForces<N>
    {
        self.calculate_nodal_forces_linear_element(nodal_coordinates)
    }
    fn calculate_nodal_stiffnesses(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalStiffnesses<N>
    {
        self.calculate_nodal_stiffnesses_linear_element(nodal_coordinates)
    }
}

impl<'a, C> ElasticLinearElement<'a, C, G, M, N> for Triangle<C>
where
    C: Elastic<'a>
{}

impl<'a, C> HyperelasticFiniteElement<'a, C, G, N> for Triangle<C>
where
    C: Hyperelastic<'a>
{
    fn calculate_helmholtz_free_energy(&self, nodal_coordinates: &NodalCoordinates<N>) -> Scalar
    {
        self.calculate_helmholtz_free_energy_linear_element(nodal_coordinates)
    }
}

impl<'a, C> HyperelasticLinearElement<'a, C, G, M, N> for Triangle<C>
where
    C: Hyperelastic<'a>
{}

impl<'a, C> ViscoelasticFiniteElement<'a, C, G, N> for Triangle<C>
where
    C: Viscoelastic<'a>
{
    fn calculate_nodal_forces(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> NodalForces<N>
    {
        self.calculate_nodal_forces_linear_element(nodal_coordinates, nodal_velocities)
    }
    fn calculate_nodal_stiffnesses(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> NodalStiffnesses<N>
    {
        self.calculate_nodal_stiffnesses_linear_element(nodal_coordinates, nodal_velocities)
    }
}

impl<'a, C> ElasticHyperviscousFiniteElement<'a, C, G, N> for Triangle<C>
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

impl<'a, C> HyperviscoelasticFiniteElement<'a, C, G, N> for Triangle<C>
where
    C: Hyperviscoelastic<'a>
{
    fn calculate_helmholtz_free_energy(&self, nodal_coordinates: &NodalCoordinates<N>) -> Scalar
    {
        self.calculate_helmholtz_free_energy_linear_element(nodal_coordinates)
    }
}

impl<'a, C> ViscoelasticLinearElement<'a, C, G, M, N> for Triangle<C>
where
    C: Viscoelastic<'a>
{}

impl<'a, C> ElasticHyperviscousLinearElement<'a, C, G, M, N> for Triangle<C>
where
    C: ElasticHyperviscous<'a>
{}

impl<'a, C> HyperviscoelasticLinearElement<'a, C, G, M, N> for Triangle<C>
where
    C: Hyperviscoelastic<'a>
{}