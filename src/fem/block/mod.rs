#[cfg(test)]
mod test;

pub mod element;

use super::*;
use self::element::
{
    FiniteElement,
    ElasticFiniteElement,
    HyperelasticFiniteElement
};
use std::array::from_fn;

pub struct Block<const D: usize, const E: usize, F, const G: usize, const N: usize>
{
    connectivity: Connectivity<E, N>,
    current_nodal_coordinates: CurrentNodalCoordinates<D>,
    elements: [F; E]
}

pub struct _ThermalSolidBlock<const D: usize, const E: usize, F, const G: usize, const N: usize>
{
    connectivity: Connectivity<E, N>,
    current_nodal_coordinates: CurrentNodalCoordinates<D>,
    elements: [F; E],
    nodal_temperatures: _NodalTemperatures<D>
}

pub trait FiniteElementBlock<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
where
    C: ConstitutiveModel<'a>,
    F: FiniteElement<'a, C, G, N>
{
    fn calculate_current_nodal_coordinates_element(&self, element_connectivity: &[usize; N]) -> CurrentNodalCoordinates<N>;
    fn get_current_nodal_coordinates(&self) -> &CurrentNodalCoordinates<D>;
    fn get_connectivity(&self) -> &Connectivity<E, N>;
    fn get_elements(&self) -> &[F; E];
    fn new(constitutive_model_parameters: ConstitutiveModelParameters<'a>, connectivity: Connectivity<E, N>, reference_nodal_coordinates: ReferenceNodalCoordinates<D>) -> Self;
    fn set_current_nodal_coordinates(&mut self, current_nodal_coordinates: CurrentNodalCoordinates<D>);
}

pub trait ThermalSolidFiniteElementBlock<'a, C, C1, C2, const D: usize, const E: usize, F, const G: usize, const N: usize>
where
    C: SolidThermal<'a, C1, C2>,
    C1: Solid<'a>,
    C2: Thermal<'a>,
    F: FiniteElement<'a, C, G, N>,
    Self: FiniteElementBlock<'a, C, D, E, F, G, N>
{
    fn get_nodal_temperatures(&self) -> &_NodalTemperatures<D>;
    fn set_nodal_temperatures(&mut self, nodal_temperatures: _NodalTemperatures<D>);
}

pub trait ElasticFiniteElementBlock<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
where
    C: Elastic<'a>,
    F: ElasticFiniteElement<'a, C, G, N>,
    Self: FiniteElementBlock<'a, C, D, E, F, G, N>
{
    fn calculate_nodal_forces(&self) -> NodalForces<D>;
    fn calculate_nodal_stiffnesses(&self) -> NodalStiffnesses<D>;
}

pub trait HyperelasticFiniteElementBlock<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
where
    C: Hyperelastic<'a>,
    F: HyperelasticFiniteElement<'a, C, G, N>,
    Self: ElasticFiniteElementBlock<'a, C, D, E, F, G, N>
{
    fn calculate_helmholtz_free_energy(&self) -> Scalar;
}

impl<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
    FiniteElementBlock<'a, C, D, E, F, G, N>
    for Block<D, E, F, G, N>
where
    C: ConstitutiveModel<'a>,
    F: FiniteElement<'a, C, G, N>
{
    fn calculate_current_nodal_coordinates_element(&self, element_connectivity: &[usize; N]) -> CurrentNodalCoordinates<N>
    {
        let current_nodal_coordinates = self.get_current_nodal_coordinates();
        element_connectivity.iter().map(|node|
            current_nodal_coordinates[*node]
            .iter().copied().collect()
        ).collect()
    }
    fn get_current_nodal_coordinates(&self) -> &CurrentNodalCoordinates<D>
    {
        &self.current_nodal_coordinates
    }
    fn get_connectivity(&self) -> &Connectivity<E, N>
    {
        &self.connectivity
    }
    fn get_elements(&self) -> &[F; E]
    {
        &self.elements
    }
    fn new(constitutive_model_parameters: ConstitutiveModelParameters<'a>, connectivity: Connectivity<E, N>, reference_nodal_coordinates: ReferenceNodalCoordinates<D>) -> Self
    {
        let elements = from_fn(|element|
            <F>::new(
                constitutive_model_parameters,
                connectivity[element].iter().map(|node|
                    reference_nodal_coordinates[*node]
                    .iter().copied().collect()
                ).collect()
            )
        );
        Self
        {
            connectivity,
            current_nodal_coordinates: reference_nodal_coordinates.convert(),
            elements
        }
    }
    fn set_current_nodal_coordinates(&mut self, current_nodal_coordinates: CurrentNodalCoordinates<D>)
    {
        self.current_nodal_coordinates = current_nodal_coordinates;
    }
}

impl<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
    ElasticFiniteElementBlock<'a, C, D, E, F, G, N>
    for Block<D, E, F, G, N>
where
    C: Elastic<'a>,
    F: ElasticFiniteElement<'a, C, G, N>,
    Self: FiniteElementBlock<'a, C, D, E, F, G, N>
{
    fn calculate_nodal_forces(&self) -> NodalForces<D>
    {
        let mut nodal_forces = NodalForces::zero();
        self.get_elements().iter()
        .zip(self.get_connectivity().iter())
        .for_each(|(element, element_connectivity)|
            element.calculate_nodal_forces(
                &self.calculate_current_nodal_coordinates_element(element_connectivity)
            ).iter()
            .zip(element_connectivity.iter())
            .for_each(|(nodal_force, node)|
                nodal_forces[*node] += nodal_force
            )
        );
        nodal_forces
    }
    fn calculate_nodal_stiffnesses(&self) -> NodalStiffnesses<D>
    {
        let mut nodal_stiffnesses = NodalStiffnesses::zero();
        self.get_elements().iter()
        .zip(self.get_connectivity().iter())
        .for_each(|(element, element_connectivity)|
            element.calculate_nodal_stiffnesses(
                &self.calculate_current_nodal_coordinates_element(element_connectivity)
            ).iter()
            .zip(element_connectivity.iter())
            .for_each(|(object, node_a)|
                object.iter()
                .zip(element_connectivity.iter())
                .for_each(|(nodal_stiffness, node_b)|
                    nodal_stiffnesses[*node_a][*node_b] += nodal_stiffness
                )
            )
        );
        nodal_stiffnesses
    }
}

impl<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
    HyperelasticFiniteElementBlock<'a, C, D, E, F, G, N>
    for Block<D, E, F, G, N>
where
    C: Hyperelastic<'a>,
    F: HyperelasticFiniteElement<'a, C, G, N>,
    Self: ElasticFiniteElementBlock<'a, C, D, E, F, G, N>
{
    fn calculate_helmholtz_free_energy(&self) -> Scalar
    {
        self.get_elements().iter()
        .zip(self.get_connectivity().iter())
        .map(|(element, element_connectivity)|
            element.calculate_helmholtz_free_energy(
                &self.calculate_current_nodal_coordinates_element(element_connectivity)
            )
        ).sum()
    }
}