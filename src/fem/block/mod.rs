#[cfg(test)]
mod test;

pub mod element;

use crate::
{
    ABS_TOL,
    REL_TOL,
    constitutive::
    {
        ConstitutiveModel,
        ConstitutiveModelParameters,
        hyperelastic::HyperelasticConstitutiveModel
    },
    math::
    {
        ContractSecondFourthIndicesWithFirstIndicesOf,
        Convert,
        TensorRank0ListTrait,
        TensorRank1ListTrait,
        TensorRank2Trait,
        TensorRank2List2DTrait
    },
    mechanics::
    {
        CurrentCoordinates,
        DeformationGradient,
        FirstPiolaKirchoffStress,
        FirstPiolaKirchoffTangentStiffness,
        Forces,
        ReferenceCoordinates,
        Scalar,
        Scalars,
        Stiffnesses,
        Vectors
    }
};
use self::element::FiniteElement;
use std::array::from_fn;

type Connectivity<const E: usize, const N: usize> = [[usize; N]; E];
type CurrentNodalCoordinates<const D: usize> = CurrentCoordinates<D>;
type NodalForces<const D: usize> = Forces<D>;
type NodalStiffnesses<const D: usize> = Stiffnesses<D>;
type ReferenceNodalCoordinates<const D: usize> = ReferenceCoordinates<D>;

pub struct FiniteElementBlock<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
where
    C: ConstitutiveModel<'a> + HyperelasticConstitutiveModel,
    F: FiniteElement<'a, C, G, N>
{
    connectivity: Connectivity<E, N>,
    current_nodal_coordinates: CurrentNodalCoordinates<D>,
    elements: [F; E],
    phantom_a: std::marker::PhantomData<*const &'a C>
}

pub trait FiniteElementBlockTrait<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
where
    C: ConstitutiveModel<'a> + HyperelasticConstitutiveModel,
    F: FiniteElement<'a, C, G, N>
{
    fn calculate_current_nodal_coordinates_element(&self, element_connectivity: &[usize; N]) -> CurrentNodalCoordinates<N>;
    fn calculate_helmholtz_free_energy(&self) -> Scalar;
    fn calculate_nodal_forces(&self) -> NodalForces<D>;
    fn calculate_nodal_stiffnesses(&self) -> NodalStiffnesses<D>;
    fn get_current_nodal_coordinates(&self) -> &CurrentNodalCoordinates<D>;
    fn get_connectivity(&self) -> &Connectivity<E, N>;
    fn get_elements(&self) -> &[F; E];
    fn new(constitutive_model_parameters: ConstitutiveModelParameters<'a>, connectivity: Connectivity<E, N>, reference_nodal_coordinates: ReferenceNodalCoordinates<D>) -> Self;
    fn set_current_nodal_coordinates(&mut self, current_nodal_coordinates: CurrentNodalCoordinates<D>);
    fn solve_using_gradient_descent(&mut self);
    fn solve_using_newtons_method(&mut self);
}

impl<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
    FiniteElementBlockTrait<'a, C, D, E, F, G, N>
    for FiniteElementBlock<'a, C, D, E, F, G, N>
where
    C: ConstitutiveModel<'a> + HyperelasticConstitutiveModel,
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
            elements,
            phantom_a: std::marker::PhantomData
        }
    }
    fn set_current_nodal_coordinates(&mut self, current_nodal_coordinates: CurrentNodalCoordinates<D>)
    {
        self.current_nodal_coordinates = current_nodal_coordinates;
    }
    fn solve_using_gradient_descent(&mut self)
    {
        let mut step_size = REL_TOL;
        let mut coordinates = self.get_current_nodal_coordinates() * 1.0;
        let mut nodal_forces = self.calculate_nodal_forces();
        let mut residual = nodal_forces.norm();
        let mut coordinates_old: CurrentNodalCoordinates<D>;
        let mut nodal_forces_old: NodalForces<D>;
        let mut nodal_forces_diff: NodalForces<D>;
        while residual > ABS_TOL
        {
            coordinates_old = &coordinates * 1.0;
            nodal_forces_old = &nodal_forces * 1.0;
            self.set_current_nodal_coordinates(
                nodal_forces * (-step_size) + coordinates
            );
            coordinates = self.get_current_nodal_coordinates() * 1.0;
            nodal_forces = self.calculate_nodal_forces();
            nodal_forces_diff = nodal_forces_old - &nodal_forces;
            step_size = (
                (coordinates_old - &coordinates).dot(&nodal_forces_diff)
            ).abs() / nodal_forces_diff.dot_self();
            residual = nodal_forces.norm();
        }
    }
    fn solve_using_newtons_method(&mut self)
    {
        let mut nodal_stiffnesses = self.calculate_nodal_stiffnesses();
        // nodal_stiffnesses.inverse();
        todo!()
    }
}