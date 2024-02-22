#[cfg(test)]
mod test;

use super::*;

const G: usize = 1;
const N: usize = 4;

pub struct LinearTetrahedron<C>
{
    constitutive_models: [C; G],
    gradient_vectors: GradientVectors<N>
}

impl<'a, C> LinearFiniteElement<'a, C, G, N> for LinearTetrahedron<C>
where
    C: Constitutive<'a>
{
    fn calculate_standard_gradient_operator() -> StandardGradientOperator<N>
    {
        StandardGradientOperator::new([
            [-1.0, -1.0, -1.0],
            [ 1.0,  0.0,  0.0],
            [ 0.0,  1.0,  0.0],
            [ 0.0,  0.0,  1.0]
        ])
    }
    fn get_gradient_vectors(&self) -> &GradientVectors<N>
    {
        &self.gradient_vectors
    }
}

impl<'a, C> FiniteElement<'a, C, G, N> for LinearTetrahedron<C>
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
            gradient_vectors: Self::calculate_gradient_vectors(&reference_nodal_coordinates)
        }
    }
}

impl<'a, C> ElasticFiniteElement<'a, C, G, N> for LinearTetrahedron<C>
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

impl<'a, C> HyperelasticFiniteElement<'a, C, G, N> for LinearTetrahedron<C>
where
    C: Hyperelastic<'a>
{
    fn calculate_helmholtz_free_energy(&self, nodal_coordinates: &NodalCoordinates<N>) -> Scalar
    {
        self.calculate_helmholtz_free_energy_linear_element(nodal_coordinates)
    }
}

impl<'a, C> ElasticLinearFiniteElement<'a, C, G, N> for LinearTetrahedron<C>
where
    C: Elastic<'a>
{}

impl<'a, C> HyperelasticLinearFiniteElement<'a, C, G, N> for LinearTetrahedron<C>
where
    C: Hyperelastic<'a>
{}

impl<'a, C> ViscoelasticFiniteElement<'a, C, G, N> for LinearTetrahedron<C>
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

impl<'a, C> HyperviscoelasticFiniteElement<'a, C, G, N> for LinearTetrahedron<C>
where
    C: Hyperviscoelastic<'a>
{
    fn calculate_helmholtz_free_energy(&self, nodal_coordinates: &NodalCoordinates<N>) -> Scalar
    {
        self.calculate_helmholtz_free_energy_linear_element(nodal_coordinates)
    }
    fn calculate_viscous_dissipation(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> Scalar
    {
        self.calculate_viscous_dissipation_linear_element(nodal_coordinates, nodal_velocities)
    }
    fn calculate_dissipation_potential(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> Scalar
    {
        self.calculate_dissipation_potential_linear_element(nodal_coordinates, nodal_velocities)
    }
}

impl<'a, C> ViscoelasticLinearFiniteElement<'a, C, G, N> for LinearTetrahedron<C>
where
    C: Viscoelastic<'a>
{}

impl<'a, C> HyperviscoelasticLinearFiniteElement<'a, C, G, N> for LinearTetrahedron<C>
where
    C: Hyperviscoelastic<'a>
{}