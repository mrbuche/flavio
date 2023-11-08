#[cfg(test)]
mod test;

use super::*;
//
// for get_integration_weights():
//
// these are the same for each element in the block
// so should not store them as a field
// this is perhaps not the most efficient way to get them
// but will work for now in terms of reducing memory stored for each element
//
// if integration_points are used in parametric space
// should do something similar

const NUMBER_OF_NODES: usize = 4;

pub struct LinearTetrahedron<'a, C, const G: usize>
where
    C: ConstitutiveModel<'a>
{
    constitutive_models: [C; G],
    gradient_vectors: GradientVectors<NUMBER_OF_NODES>,
    phantom_a: std::marker::PhantomData<*const &'a C>
}

impl<'a, C, const G: usize> LinearFiniteElement<'a, C, G, NUMBER_OF_NODES> for LinearTetrahedron<'a, C, G>
where
    C: ConstitutiveModel<'a>
{
    fn calculate_standard_gradient_operator() -> StandardGradientOperator<NUMBER_OF_NODES>
    {
        StandardGradientOperator::new([
            [-1.0, -1.0, -1.0],
            [ 1.0,  0.0,  0.0],
            [ 0.0,  1.0,  0.0],
            [ 0.0,  0.0,  1.0]
        ])
    }
    fn get_gradient_vectors(&self) -> &GradientVectors<NUMBER_OF_NODES>
    {
        &self.gradient_vectors
    }
}

impl<'a, C, const G: usize> FiniteElement<'a, C, G, NUMBER_OF_NODES> for LinearTetrahedron<'a, C, G>
where
    C: ConstitutiveModel<'a>
{
    fn calculate_deformation_gradients(&self, _current_nodal_coordinates: &CurrentNodalCoordinates<NUMBER_OF_NODES>) -> DeformationGradients<G>
    {
        panic!()
    }
    fn calculate_helmholtz_free_energy(&self, current_nodal_coordinates: &CurrentNodalCoordinates<NUMBER_OF_NODES>) -> Scalar
    {
        self.calculate_helmholtz_free_energy_linear_element(current_nodal_coordinates)
    }
    fn calculate_nodal_forces(&self, current_nodal_coordinates: &CurrentNodalCoordinates<NUMBER_OF_NODES>) -> NodalForces<NUMBER_OF_NODES>
    {
        self.calculate_nodal_forces_linear_element(current_nodal_coordinates)
    }
    fn calculate_nodal_stiffnesses(&self, current_nodal_coordinates: &CurrentNodalCoordinates<NUMBER_OF_NODES>) -> NodalStiffnesses<NUMBER_OF_NODES>
    {
        self.calculate_nodal_stiffnesses_linear_element(current_nodal_coordinates)
    }
    fn get_constitutive_models(&self) -> &[C; G]
    {
        &self.constitutive_models
    }
    fn get_integration_weights(&self) -> IntegrationWeights<G>
    {
        let mut integration_weights =
        {
            if G == 1
            {
                IntegrationWeights::new([1.0; G])
            }
            else if G == 4
            {
                IntegrationWeights::new([0.25; G])
            }
            else if G == 5
            {
                IntegrationWeights::new([0.45; G])
            }
            else if G == 11
            {
                IntegrationWeights::new([0.448/3.0; G])
            }
            else if G == 15
            {
                IntegrationWeights::new([0.06569484936831872; G])
            }
            else
            {
                panic!("Invalid number of integration points.")
            }
        };
        if G == 5
        {
            integration_weights[0] = -0.8;
        }
        else if G == 11
        {
            integration_weights[0] = -0.2368/3.0;
            integration_weights[1] = 0.1372/3.0;
            integration_weights[2] = 0.1372/3.0;
            integration_weights[3] = 0.1372/3.0;
            integration_weights[4] = 0.1372/3.0;
        }
        else if G == 15
        {
            integration_weights[0] = 0.18170206858253513;
            integration_weights[1] = 0.036160714285714296;
            integration_weights[2] = 0.036160714285714296;
            integration_weights[3] = 0.036160714285714296;
            integration_weights[4] = 0.036160714285714296;
            integration_weights[5] = 0.06987149451617385;
            integration_weights[6] = 0.06987149451617385;
            integration_weights[7] = 0.06987149451617385;
            integration_weights[8] = 0.06987149451617385;
        }
        integration_weights
    }
    fn new(constitutive_model_parameters: ConstitutiveModelParameters<'a>, reference_nodal_coordinates: ReferenceNodalCoordinates<NUMBER_OF_NODES>) -> Self
    {
        Self
        {
            constitutive_models: std::array::from_fn(|_| <C>::new(constitutive_model_parameters)),
            gradient_vectors: Self::calculate_gradient_vectors(&reference_nodal_coordinates),
            phantom_a: std::marker::PhantomData
        }
    }
}