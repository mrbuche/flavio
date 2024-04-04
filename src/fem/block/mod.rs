#[cfg(test)]
mod test;

pub mod element;

use super::*;
use self::element::
{
    FiniteElement,
    SurfaceElement,
    ElasticFiniteElement,
    HyperelasticFiniteElement,
    ViscoelasticFiniteElement,
    ElasticHyperviscousFiniteElement,
    HyperviscoelasticFiniteElement
};
use std::array::from_fn;

pub struct ElasticBlock<const D: usize, const E: usize, F, const G: usize, const N: usize>
{
    connectivity: Connectivity<E, N>,
    elements: [F; E],
    nodal_coordinates: NodalCoordinates<D>
}

pub struct ViscoelasticBlock<const D: usize, const E: usize, F, const G: usize, const N: usize>
{
    connectivity: Connectivity<E, N>,
    elements: [F; E],
    nodal_coordinates: NodalCoordinates<D>,
    nodal_velocities: NodalVelocities<D>
}

pub trait BasicFiniteElementBlock<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
where
    C: Constitutive<'a>
{
    fn calculate_nodal_coordinates_element(&self, element_connectivity: &[usize; N]) -> NodalCoordinates<N>;
    fn get_nodal_coordinates(&self) -> &NodalCoordinates<D>;
    fn get_connectivity(&self) -> &Connectivity<E, N>;
    fn get_elements(&self) -> &[F; E];
    fn set_nodal_coordinates(&mut self, nodal_coordinates: NodalCoordinates<D>);
}

pub trait FiniteElementBlock<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
where
    C: Constitutive<'a>,
    F: FiniteElement<'a, C, G, N>,
    Self: BasicFiniteElementBlock<'a, C, D, E, F, G, N>
{
    fn new(constitutive_model_parameters: Parameters<'a>, connectivity: Connectivity<E, N>, reference_nodal_coordinates: ReferenceNodalCoordinates<D>) -> Self;
}

pub trait SurfaceElementBlock<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
where
    C: Constitutive<'a>,
    F: SurfaceElement<'a, C, G, N>,
    Self: BasicFiniteElementBlock<'a, C, D, E, F, G, N>
{
    fn new(constitutive_model_parameters: Parameters<'a>, connectivity: Connectivity<E, N>, reference_nodal_coordinates: ReferenceNodalCoordinates<D>, thickness: Scalar) -> Self;
}

pub trait ElasticFiniteElementBlock<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
where
    C: Elastic<'a>,
    F: ElasticFiniteElement<'a, C, G, N>
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

pub trait ViscoelasticFiniteElementBlock<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
where
    C: Viscoelastic<'a>,
    F: ViscoelasticFiniteElement<'a, C, G, N>
{
    fn calculate_nodal_forces(&self) -> NodalForces<D>;
    fn calculate_nodal_stiffnesses(&self) -> NodalStiffnesses<D>;
    fn calculate_nodal_velocities_element(&self, element_connectivity: &[usize; N]) -> NodalVelocities<N>;
    fn get_nodal_velocities(&self) -> &NodalVelocities<D>;
    fn set_nodal_velocities(&mut self, nodal_velocities: NodalVelocities<D>);
}

pub trait ElasticHyperviscousFiniteElementBlock<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
where
    C: ElasticHyperviscous<'a>,
    F: ElasticHyperviscousFiniteElement<'a, C, G, N>,
    Self: ViscoelasticFiniteElementBlock<'a, C, D, E, F, G, N>
{
    fn calculate_viscous_dissipation(&self) -> Scalar;
    fn calculate_dissipation_potential(&self) -> Scalar;
}

pub trait HyperviscoelasticFiniteElementBlock<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
where
    C: Hyperviscoelastic<'a>,
    F: HyperviscoelasticFiniteElement<'a, C, G, N>,
    Self: ElasticHyperviscousFiniteElementBlock<'a, C, D, E, F, G, N>
{
    fn calculate_helmholtz_free_energy(&self) -> Scalar;
}

impl<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
    BasicFiniteElementBlock<'a, C, D, E, F, G, N>
    for ElasticBlock<D, E, F, G, N>
where
    C: Elastic<'a>,
    F: ElasticFiniteElement<'a, C, G, N>
{
    fn calculate_nodal_coordinates_element(&self, element_connectivity: &[usize; N]) -> NodalCoordinates<N>
    {
        let nodal_coordinates = self.get_nodal_coordinates();
        element_connectivity.iter().map(|node|
            nodal_coordinates[*node]
            .iter().copied().collect()
        ).collect()
    }
    fn get_nodal_coordinates(&self) -> &NodalCoordinates<D>
    {
        &self.nodal_coordinates
    }
    fn get_connectivity(&self) -> &Connectivity<E, N>
    {
        &self.connectivity
    }
    fn get_elements(&self) -> &[F; E]
    {
        &self.elements
    }
    fn set_nodal_coordinates(&mut self, nodal_coordinates: NodalCoordinates<D>)
    {
        self.nodal_coordinates = nodal_coordinates;
    }
}

impl<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
    FiniteElementBlock<'a, C, D, E, F, G, N>
    for ElasticBlock<D, E, F, G, N>
where
    C: Elastic<'a>,
    F: FiniteElement<'a, C, G, N> + ElasticFiniteElement<'a, C, G, N>
{
    fn new(constitutive_model_parameters: Parameters<'a>, connectivity: Connectivity<E, N>, reference_nodal_coordinates: ReferenceNodalCoordinates<D>) -> Self
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
            elements,
            nodal_coordinates: reference_nodal_coordinates.convert()
        }
    }
}

impl<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
    SurfaceElementBlock<'a, C, D, E, F, G, N>
    for ElasticBlock<D, E, F, G, N>
where
    C: Elastic<'a>,
    F: SurfaceElement<'a, C, G, N> + ElasticFiniteElement<'a, C, G, N>
{
    fn new(constitutive_model_parameters: Parameters<'a>, connectivity: Connectivity<E, N>, reference_nodal_coordinates: ReferenceNodalCoordinates<D>, thickness: Scalar) -> Self
    {
        let elements = from_fn(|element|
            <F>::new(
                constitutive_model_parameters,
                connectivity[element].iter().map(|node|
                    reference_nodal_coordinates[*node]
                    .iter().copied().collect()
                ).collect(),
                &thickness
            )
        );
        Self
        {
            connectivity,
            elements,
            nodal_coordinates: reference_nodal_coordinates.convert()
        }
    }
}

impl<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
    ElasticFiniteElementBlock<'a, C, D, E, F, G, N>
    for ElasticBlock<D, E, F, G, N>
where
    C: Elastic<'a>,
    F: ElasticFiniteElement<'a, C, G, N>,
    Self: BasicFiniteElementBlock<'a, C, D, E, F, G, N>
{
    fn calculate_nodal_forces(&self) -> NodalForces<D>
    {
        let mut nodal_forces = NodalForces::zero();
        self.get_elements().iter()
        .zip(self.get_connectivity().iter())
        .for_each(|(element, element_connectivity)|
            element.calculate_nodal_forces(
                &self.calculate_nodal_coordinates_element(
                    element_connectivity
                )
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
                &self.calculate_nodal_coordinates_element(
                    element_connectivity
                )
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
    for ElasticBlock<D, E, F, G, N>
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
                &self.calculate_nodal_coordinates_element(
                    element_connectivity
                )
            )
        ).sum()
    }
}

impl<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
    BasicFiniteElementBlock<'a, C, D, E, F, G, N>
    for ViscoelasticBlock<D, E, F, G, N>
where
    C: Viscoelastic<'a>,
    F: ViscoelasticFiniteElement<'a, C, G, N>
{
    fn calculate_nodal_coordinates_element(&self, element_connectivity: &[usize; N]) -> NodalCoordinates<N>
    {
        let nodal_coordinates = self.get_nodal_coordinates();
        element_connectivity.iter().map(|node|
            nodal_coordinates[*node]
            .iter().copied().collect()
        ).collect()
    }
    fn get_nodal_coordinates(&self) -> &NodalCoordinates<D>
    {
        &self.nodal_coordinates
    }
    fn get_connectivity(&self) -> &Connectivity<E, N>
    {
        &self.connectivity
    }
    fn get_elements(&self) -> &[F; E]
    {
        &self.elements
    }
    fn set_nodal_coordinates(&mut self, nodal_coordinates: NodalCoordinates<D>)
    {
        self.nodal_coordinates = nodal_coordinates;
    }
}

impl<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
    FiniteElementBlock<'a, C, D, E, F, G, N>
    for ViscoelasticBlock<D, E, F, G, N>
where
    C: Viscoelastic<'a>,
    F: FiniteElement<'a, C, G, N> + ViscoelasticFiniteElement<'a, C, G, N>
{
    fn new(constitutive_model_parameters: Parameters<'a>, connectivity: Connectivity<E, N>, reference_nodal_coordinates: ReferenceNodalCoordinates<D>) -> Self
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
            elements,
            nodal_coordinates: reference_nodal_coordinates.convert(),
            nodal_velocities: NodalVelocities::zero()
        }
    }
}

impl<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
    SurfaceElementBlock<'a, C, D, E, F, G, N>
    for ViscoelasticBlock<D, E, F, G, N>
where
    C: Viscoelastic<'a>,
    F: SurfaceElement<'a, C, G, N> + ViscoelasticFiniteElement<'a, C, G, N>
{
    fn new(constitutive_model_parameters: Parameters<'a>, connectivity: Connectivity<E, N>, reference_nodal_coordinates: ReferenceNodalCoordinates<D>, thickness: Scalar) -> Self
    {
        let elements = from_fn(|element|
            <F>::new(
                constitutive_model_parameters,
                connectivity[element].iter().map(|node|
                    reference_nodal_coordinates[*node]
                    .iter().copied().collect()
                ).collect(),
                &thickness
            )
        );
        Self
        {
            connectivity,
            elements,
            nodal_coordinates: reference_nodal_coordinates.convert(),
            nodal_velocities: NodalVelocities::zero()
        }
    }
}

impl<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
    ViscoelasticFiniteElementBlock<'a, C, D, E, F, G, N>
    for ViscoelasticBlock<D, E, F, G, N>
where
    C: Viscoelastic<'a>,
    F: ViscoelasticFiniteElement<'a, C, G, N>,
    Self: BasicFiniteElementBlock<'a, C, D, E, F, G, N>
{
    fn calculate_nodal_forces(&self) -> NodalForces<D>
    {
        let mut nodal_forces = NodalForces::zero();
        self.get_elements().iter()
        .zip(self.get_connectivity().iter())
        .for_each(|(element, element_connectivity)|
            element.calculate_nodal_forces(
                &self.calculate_nodal_coordinates_element(
                    element_connectivity
                ),
                &self.calculate_nodal_velocities_element(
                    element_connectivity
                )
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
                &self.calculate_nodal_coordinates_element(
                    element_connectivity
                ),
                &self.calculate_nodal_velocities_element(
                    element_connectivity
                )
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
    fn calculate_nodal_velocities_element(&self, element_connectivity: &[usize; N]) -> NodalVelocities<N>
    {
        let nodal_velocities = self.get_nodal_velocities();
        element_connectivity.iter().map(|node|
            nodal_velocities[*node]
            .iter().copied().collect()
        ).collect()
    }
    fn get_nodal_velocities(&self) -> &NodalVelocities<D>
    {
        &self.nodal_velocities
    }
    fn set_nodal_velocities(&mut self, nodal_velocities: NodalVelocities<D>)
    {
        self.nodal_velocities = nodal_velocities;
    }
}

impl<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
    ElasticHyperviscousFiniteElementBlock<'a, C, D, E, F, G, N>
    for ViscoelasticBlock<D, E, F, G, N>
where
    C: ElasticHyperviscous<'a>,
    F: ElasticHyperviscousFiniteElement<'a, C, G, N>,
    Self: ViscoelasticFiniteElementBlock<'a, C, D, E, F, G, N>
{
    fn calculate_viscous_dissipation(&self) -> Scalar
    {
        self.get_elements().iter()
        .zip(self.get_connectivity().iter())
        .map(|(element, element_connectivity)|
            element.calculate_viscous_dissipation(
                &self.calculate_nodal_coordinates_element(
                    element_connectivity
                ),
                &self.calculate_nodal_velocities_element(
                    element_connectivity
                )
            )
        ).sum()
    }
    fn calculate_dissipation_potential(&self) -> Scalar
    {
        self.get_elements().iter()
        .zip(self.get_connectivity().iter())
        .map(|(element, element_connectivity)|
            element.calculate_dissipation_potential(
                &self.calculate_nodal_coordinates_element(
                    element_connectivity
                ),
                &self.calculate_nodal_velocities_element(
                    element_connectivity
                )
            )
        ).sum()
    }
}

impl<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
    HyperviscoelasticFiniteElementBlock<'a, C, D, E, F, G, N>
    for ViscoelasticBlock<D, E, F, G, N>
where
    C: Hyperviscoelastic<'a>,
    F: HyperviscoelasticFiniteElement<'a, C, G, N>,
    Self: ElasticHyperviscousFiniteElementBlock<'a, C, D, E, F, G, N>
{
    fn calculate_helmholtz_free_energy(&self) -> Scalar
    {
        self.get_elements().iter()
        .zip(self.get_connectivity().iter())
        .map(|(element, element_connectivity)|
            element.calculate_helmholtz_free_energy(
                &self.calculate_nodal_coordinates_element(
                    element_connectivity
                )
            )
        ).sum()
    }
}