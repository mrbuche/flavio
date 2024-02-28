#[cfg(test)]
mod test;

use super::*;

const G: usize = 1;
const M: usize = 2;
const N: usize = 3;

pub struct Triangle<'a, C>
{
    constitutive_models: [C; G],
    gradient_vectors: GradientVectors<N>,
    reference_normal: ReferenceNormal,
    thickness: &'a Scalar
}

impl<'a, C> FiniteElement<'a, C, G, N> for Triangle<'a, C>
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
            reference_normal: Self::calculate_reference_normal(&Self::calculate_reference_dual_basis_vectors(&reference_nodal_coordinates)),
            thickness: &1.0
        }
    }
}

impl<'a, C> SurfaceElement<'a> for Triangle<'a, C>
where
    C: Constitutive<'a>
{
    fn get_thickness(&self) -> &Scalar
    {
        self.thickness
    }
    fn set_thickness(&mut self, thickness: &'a Scalar)
    {
        self.thickness = thickness;
    }
}

impl<'a, C> LinearElement<'a, C, G, M, N> for Triangle<'a, C>
where
    C: Constitutive<'a>
{
    fn calculate_deformation_gradient(&self, nodal_coordinates: &NodalCoordinates<N>) -> DeformationGradient
    {
        self.calculate_deformation_gradient_linear_surface_element(nodal_coordinates)
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

impl<'a, C> LinearSurfaceElement<'a, C, G, M, N> for Triangle<'a, C>
where
    C: Constitutive<'a>
{
    fn get_reference_normal(&self) -> &ReferenceNormal
    {
        &self.reference_normal
    }
    fn get_thickness(&self) -> &Scalar
    {
        self.thickness
    }
}

impl<'a, C> ElasticFiniteElement<'a, C, G, N> for Triangle<'a, C>
where
    C: Elastic<'a>
{
    fn calculate_nodal_forces(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalForces<N>
    {
        self.calculate_nodal_forces_linear_surface_element(nodal_coordinates)
    }
    fn calculate_nodal_stiffnesses(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalStiffnesses<N>
    {
        self.calculate_nodal_stiffnesses_linear_surface_element(nodal_coordinates)
    }
}

impl<'a, C> ElasticLinearSurfaceElement<'a, C, G, M, N> for Triangle<'a, C>
where
    C: Elastic<'a>
{}

impl<'a, C> HyperelasticFiniteElement<'a, C, G, N> for Triangle<'a, C>
where
    C: Hyperelastic<'a>
{
    fn calculate_helmholtz_free_energy(&self, nodal_coordinates: &NodalCoordinates<N>) -> Scalar
    {
        self.calculate_helmholtz_free_energy_linear_surface_element(nodal_coordinates)
    }
}

impl<'a, C> HyperelasticLinearSurfaceElement<'a, C, G, M, N> for Triangle<'a, C>
where
    C: Hyperelastic<'a>
{}