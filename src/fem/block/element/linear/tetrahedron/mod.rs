#[cfg(test)]
mod test;

use super::*;

const G: usize = 1;
const N: usize = 4;

pub struct LinearTetrahedron<'a, C>
where
    C: ConstitutiveModel<'a>
{
    constitutive_models: [C; G],
    gradient_vectors: GradientVectors<N>,
    phantom_a: std::marker::PhantomData<*const &'a C>
}

impl<'a, C> LinearFiniteElement<'a, C, G, N> for LinearTetrahedron<'a, C>
where
    C: ConstitutiveModel<'a>
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

impl<'a, C> FiniteElement<'a, C, G, N> for LinearTetrahedron<'a, C>
where
    C: ConstitutiveModel<'a>
{
    fn get_constitutive_models(&self) -> &[C; G]
    {
        &self.constitutive_models
    }
    fn get_integration_weights(&self) -> IntegrationWeights<G>
    {
        IntegrationWeights::new([1.0; G])
    }
    fn new(constitutive_model_parameters: ConstitutiveModelParameters<'a>, reference_nodal_coordinates: ReferenceNodalCoordinates<N>) -> Self
    {
        Self
        {
            constitutive_models: std::array::from_fn(|_| <C>::new(constitutive_model_parameters)),
            gradient_vectors: Self::calculate_gradient_vectors(&reference_nodal_coordinates),
            phantom_a: std::marker::PhantomData
        }
    }
}

impl<'a, C> ElasticFiniteElement<'a, C, G, N> for LinearTetrahedron<'a, C>
where
    C: ElasticConstitutiveModel<'a>
{
    fn calculate_nodal_forces(&self, current_nodal_coordinates: &CurrentNodalCoordinates<N>) -> NodalForces<N>
    {
        self.calculate_nodal_forces_linear_element(current_nodal_coordinates)
    }
    fn calculate_nodal_stiffnesses(&self, current_nodal_coordinates: &CurrentNodalCoordinates<N>) -> NodalStiffnesses<N>
    {
        self.calculate_nodal_stiffnesses_linear_element(current_nodal_coordinates)
    }
}

impl<'a, C> HyperelasticFiniteElement<'a, C, G, N> for LinearTetrahedron<'a, C>
where
    C: HyperelasticConstitutiveModel<'a>
{
    fn calculate_helmholtz_free_energy(&self, current_nodal_coordinates: &CurrentNodalCoordinates<N>) -> Scalar
    {
        self.calculate_helmholtz_free_energy_linear_element(current_nodal_coordinates)
    }
}

impl<'a, C> ElasticLinearFiniteElement<'a, C, G, N> for LinearTetrahedron<'a, C>
where
    C: ElasticConstitutiveModel<'a>
{}

impl<'a, C> HyperelasticLinearFiniteElement<'a, C, G, N> for LinearTetrahedron<'a, C>
where
    C: HyperelasticConstitutiveModel<'a>
{}