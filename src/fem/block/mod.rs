#[cfg(test)]
mod test;

pub mod element;

use self::element::{
    ElasticFiniteElement, ElasticHyperviscousFiniteElement, FiniteElement,
    HyperelasticFiniteElement, HyperviscoelasticFiniteElement, SurfaceElement,
    ViscoelasticFiniteElement,
};
use super::*;
use std::array::from_fn;

pub struct ElasticBlock<const D: usize, const E: usize, F, const G: usize, const N: usize> {
    connectivity: Connectivity<E, N>,
    elements: [F; E],
    nodal_coordinates: NodalCoordinates<D>,
}

pub struct ViscoelasticBlock<const D: usize, const E: usize, F, const G: usize, const N: usize> {
    connectivity: Connectivity<E, N>,
    elements: [F; E],
    nodal_coordinates: NodalCoordinates<D>,
    nodal_velocities: NodalVelocities<D>,
}

pub trait BasicFiniteElementBlock<
    'a,
    C,
    const D: usize,
    const E: usize,
    F,
    const G: usize,
    const N: usize,
> where
    C: Constitutive<'a>,
{
    fn calculate_nodal_coordinates_element(
        &self,
        element_connectivity: &[usize; N],
    ) -> NodalCoordinates<N>;
    fn get_nodal_coordinates(&self) -> &NodalCoordinates<D>;
    fn get_connectivity(&self) -> &Connectivity<E, N>;
    fn get_elements(&self) -> &[F; E];
    fn set_nodal_coordinates(&mut self, nodal_coordinates: NodalCoordinates<D>);
}

pub trait FiniteElementBlock<
    'a,
    C,
    const D: usize,
    const E: usize,
    F,
    const G: usize,
    const N: usize,
> where
    C: Constitutive<'a>,
    F: FiniteElement<'a, C, G, N>,
    Self: BasicFiniteElementBlock<'a, C, D, E, F, G, N>,
{
    fn new(
        constitutive_model_parameters: Parameters<'a>,
        connectivity: Connectivity<E, N>,
        reference_nodal_coordinates: ReferenceNodalCoordinates<D>,
    ) -> Self;
}

pub trait SurfaceElementBlock<
    'a,
    C,
    const D: usize,
    const E: usize,
    F,
    const G: usize,
    const N: usize,
> where
    C: Constitutive<'a>,
    F: SurfaceElement<'a, C, G, N>,
    Self: BasicFiniteElementBlock<'a, C, D, E, F, G, N>,
{
    fn new(
        constitutive_model_parameters: Parameters<'a>,
        connectivity: Connectivity<E, N>,
        reference_nodal_coordinates: ReferenceNodalCoordinates<D>,
        thickness: Scalar,
    ) -> Self;
}

pub trait ElasticFiniteElementBlock<
    'a,
    C,
    const D: usize,
    const E: usize,
    F,
    const G: usize,
    const N: usize,
> where
    C: Elastic<'a>,
    F: ElasticFiniteElement<'a, C, G, N>,
{
    fn calculate_nodal_forces(&self) -> Result<NodalForces<D>, ConstitutiveError>;
    fn calculate_nodal_stiffnesses(&self) -> Result<NodalStiffnesses<D>, ConstitutiveError>;
    fn solve(&mut self, fixed_nodes: &[usize]) -> Result<(), ConstitutiveError>;
}

pub trait HyperelasticFiniteElementBlock<
    'a,
    C,
    const D: usize,
    const E: usize,
    F,
    const G: usize,
    const N: usize,
> where
    C: Hyperelastic<'a>,
    F: HyperelasticFiniteElement<'a, C, G, N>,
    Self: ElasticFiniteElementBlock<'a, C, D, E, F, G, N>,
{
    fn calculate_helmholtz_free_energy(&self) -> Result<Scalar, ConstitutiveError>;
}

pub trait ViscoelasticFiniteElementBlock<
    'a,
    C,
    const D: usize,
    const E: usize,
    F,
    const G: usize,
    const N: usize,
> where
    C: Viscoelastic<'a>,
    F: ViscoelasticFiniteElement<'a, C, G, N>,
{
    fn calculate_nodal_forces(&self) -> Result<NodalForces<D>, ConstitutiveError>;
    fn calculate_nodal_stiffnesses(&self) -> Result<NodalStiffnesses<D>, ConstitutiveError>;
    fn calculate_nodal_velocities_element(
        &self,
        element_connectivity: &[usize; N],
    ) -> NodalVelocities<N>;
    fn get_nodal_velocities(&self) -> &NodalVelocities<D>;
    fn set_nodal_velocities(&mut self, nodal_velocities: NodalVelocities<D>);
}

pub trait ElasticHyperviscousFiniteElementBlock<
    'a,
    C,
    const D: usize,
    const E: usize,
    F,
    const G: usize,
    const N: usize,
> where
    C: ElasticHyperviscous<'a>,
    F: ElasticHyperviscousFiniteElement<'a, C, G, N>,
    Self: ViscoelasticFiniteElementBlock<'a, C, D, E, F, G, N>,
{
    fn calculate_viscous_dissipation(&self) -> Result<Scalar, ConstitutiveError>;
    fn calculate_dissipation_potential(&self) -> Result<Scalar, ConstitutiveError>;
}

pub trait HyperviscoelasticFiniteElementBlock<
    'a,
    C,
    const D: usize,
    const E: usize,
    F,
    const G: usize,
    const N: usize,
> where
    C: Hyperviscoelastic<'a>,
    F: HyperviscoelasticFiniteElement<'a, C, G, N>,
    Self: ElasticHyperviscousFiniteElementBlock<'a, C, D, E, F, G, N>,
{
    fn calculate_helmholtz_free_energy(&self) -> Result<Scalar, ConstitutiveError>;
}

impl<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
    BasicFiniteElementBlock<'a, C, D, E, F, G, N> for ElasticBlock<D, E, F, G, N>
where
    C: Elastic<'a>,
    F: ElasticFiniteElement<'a, C, G, N>,
{
    fn calculate_nodal_coordinates_element(
        &self,
        element_connectivity: &[usize; N],
    ) -> NodalCoordinates<N> {
        let nodal_coordinates = self.get_nodal_coordinates();
        element_connectivity
            .iter()
            .map(|node| nodal_coordinates[*node].copy())
            .collect()
    }
    fn get_nodal_coordinates(&self) -> &NodalCoordinates<D> {
        &self.nodal_coordinates
    }
    fn get_connectivity(&self) -> &Connectivity<E, N> {
        &self.connectivity
    }
    fn get_elements(&self) -> &[F; E] {
        &self.elements
    }
    fn set_nodal_coordinates(&mut self, nodal_coordinates: NodalCoordinates<D>) {
        self.nodal_coordinates = nodal_coordinates;
    }
}

impl<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
    FiniteElementBlock<'a, C, D, E, F, G, N> for ElasticBlock<D, E, F, G, N>
where
    C: Elastic<'a>,
    F: FiniteElement<'a, C, G, N> + ElasticFiniteElement<'a, C, G, N>,
{
    fn new(
        constitutive_model_parameters: Parameters<'a>,
        connectivity: Connectivity<E, N>,
        reference_nodal_coordinates: ReferenceNodalCoordinates<D>,
    ) -> Self {
        let elements = from_fn(|element| {
            <F>::new(
                constitutive_model_parameters,
                connectivity[element]
                    .iter()
                    .map(|node| reference_nodal_coordinates[*node].copy())
                    .collect(),
            )
        });
        Self {
            connectivity,
            elements,
            nodal_coordinates: reference_nodal_coordinates.into(),
        }
    }
}

impl<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
    SurfaceElementBlock<'a, C, D, E, F, G, N> for ElasticBlock<D, E, F, G, N>
where
    C: Elastic<'a>,
    F: SurfaceElement<'a, C, G, N> + ElasticFiniteElement<'a, C, G, N>,
{
    fn new(
        constitutive_model_parameters: Parameters<'a>,
        connectivity: Connectivity<E, N>,
        reference_nodal_coordinates: ReferenceNodalCoordinates<D>,
        thickness: Scalar,
    ) -> Self {
        let elements = from_fn(|element| {
            <F>::new(
                constitutive_model_parameters,
                connectivity[element]
                    .iter()
                    .map(|node| reference_nodal_coordinates[*node].copy())
                    .collect(),
                &thickness,
            )
        });
        Self {
            connectivity,
            elements,
            nodal_coordinates: reference_nodal_coordinates.into(),
        }
    }
}

impl<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
    ElasticFiniteElementBlock<'a, C, D, E, F, G, N> for ElasticBlock<D, E, F, G, N>
where
    C: Elastic<'a>,
    F: ElasticFiniteElement<'a, C, G, N>,
    Self: BasicFiniteElementBlock<'a, C, D, E, F, G, N>,
{
    fn calculate_nodal_forces(&self) -> Result<NodalForces<D>, ConstitutiveError> {
        let mut nodal_forces = NodalForces::zero();
        self.get_elements()
            .iter()
            .zip(self.get_connectivity().iter())
            .try_for_each(|(element, element_connectivity)| {
                element
                    .calculate_nodal_forces(
                        &self.calculate_nodal_coordinates_element(element_connectivity),
                    )?
                    .iter()
                    .zip(element_connectivity.iter())
                    .for_each(|(nodal_force, node)| nodal_forces[*node] += nodal_force);
                Ok(())
            })?;
        Ok(nodal_forces)
    }
    fn calculate_nodal_stiffnesses(&self) -> Result<NodalStiffnesses<D>, ConstitutiveError> {
        let mut nodal_stiffnesses = NodalStiffnesses::zero();
        self.get_elements()
            .iter()
            .zip(self.get_connectivity().iter())
            .try_for_each(|(element, element_connectivity)| {
                element
                    .calculate_nodal_stiffnesses(
                        &self.calculate_nodal_coordinates_element(element_connectivity),
                    )?
                    .iter()
                    .zip(element_connectivity.iter())
                    .for_each(|(object, node_a)| {
                        object.iter().zip(element_connectivity.iter()).for_each(
                            |(nodal_stiffness, node_b)| {
                                nodal_stiffnesses[*node_a][*node_b] += nodal_stiffness
                            },
                        )
                    });
                Ok(())
            })?;
        Ok(nodal_stiffnesses)
    }
    //
    // What about constrained optimization instead?
    // Would that work nicely, and be more robust to large step deformations?
    //
    fn solve(&mut self, fixed_nodes: &[usize]) -> Result<(), ConstitutiveError> {
        //
        // can you just make the calculate() methods take a reference rather than reference self?
        // as these solvers show, you in many cases want to compute states that are not the current one
        // and then you can use the templated solvers
        //
        // prescribed BCs
        // &[usize] for node ids
        // &[Option<f64>; 3] for values (None = free DOF)
        //
        let mut guess = self.get_nodal_coordinates().copy();
        let mut guess_0 = NodalCoordinates::zero();
        let mut residual;
        let mut residual_0 = NodalForces::zero();
        let mut residual_difference;
        let mut residual_norm: Scalar;
        let mut step_size: Scalar;
        for step in 0..1000 {
            residual = self.calculate_nodal_forces()?;
            fixed_nodes
                .iter()
                .for_each(|node| residual[*node].iter_mut().for_each(|entry| *entry = 0.0));
            residual_norm = residual
                .iter()
                .map(|entry| entry.norm_squared())
                .sum::<Scalar>()
                .sqrt();
            // if residual_norm < 1e-6 {
            if residual
                .iter()
                .filter(|entry| entry.norm() >= crate::ABS_TOL)
                .count()
                == 0
            {
                return Ok(());
            } else {
                residual_difference = residual_0 - &residual;
                //
                // how to choose short (below, dx*dg/dg*dg) or long (dx*dx/dx*dg) steps?
                //
                step_size = residual_difference.dot(&(guess_0 - &guess))
                    / residual_difference
                        .iter()
                        .map(|entry| entry.norm_squared())
                        .sum::<Scalar>();
                guess_0 = guess.copy();
                residual_0 = residual.copy();
                println!("{:?}", (step, step_size, residual_norm));
                guess = self.get_nodal_coordinates() - residual * step_size;
                self.set_nodal_coordinates(guess.copy());
            }
        }
        panic!("Maximum steps reached.")
    }
}

impl<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
    HyperelasticFiniteElementBlock<'a, C, D, E, F, G, N> for ElasticBlock<D, E, F, G, N>
where
    C: Hyperelastic<'a>,
    F: HyperelasticFiniteElement<'a, C, G, N>,
    Self: ElasticFiniteElementBlock<'a, C, D, E, F, G, N>,
{
    fn calculate_helmholtz_free_energy(&self) -> Result<Scalar, ConstitutiveError> {
        self.get_elements()
            .iter()
            .zip(self.get_connectivity().iter())
            .map(|(element, element_connectivity)| {
                element.calculate_helmholtz_free_energy(
                    &self.calculate_nodal_coordinates_element(element_connectivity),
                )
            })
            .sum()
    }
}

impl<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
    BasicFiniteElementBlock<'a, C, D, E, F, G, N> for ViscoelasticBlock<D, E, F, G, N>
where
    C: Viscoelastic<'a>,
    F: ViscoelasticFiniteElement<'a, C, G, N>,
{
    fn calculate_nodal_coordinates_element(
        &self,
        element_connectivity: &[usize; N],
    ) -> NodalCoordinates<N> {
        let nodal_coordinates = self.get_nodal_coordinates();
        element_connectivity
            .iter()
            .map(|node| nodal_coordinates[*node].copy())
            .collect()
    }
    fn get_nodal_coordinates(&self) -> &NodalCoordinates<D> {
        &self.nodal_coordinates
    }
    fn get_connectivity(&self) -> &Connectivity<E, N> {
        &self.connectivity
    }
    fn get_elements(&self) -> &[F; E] {
        &self.elements
    }
    fn set_nodal_coordinates(&mut self, nodal_coordinates: NodalCoordinates<D>) {
        self.nodal_coordinates = nodal_coordinates;
    }
}

impl<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
    FiniteElementBlock<'a, C, D, E, F, G, N> for ViscoelasticBlock<D, E, F, G, N>
where
    C: Viscoelastic<'a>,
    F: FiniteElement<'a, C, G, N> + ViscoelasticFiniteElement<'a, C, G, N>,
{
    fn new(
        constitutive_model_parameters: Parameters<'a>,
        connectivity: Connectivity<E, N>,
        reference_nodal_coordinates: ReferenceNodalCoordinates<D>,
    ) -> Self {
        let elements = from_fn(|element| {
            <F>::new(
                constitutive_model_parameters,
                connectivity[element]
                    .iter()
                    .map(|node| reference_nodal_coordinates[*node].copy())
                    .collect(),
            )
        });
        Self {
            connectivity,
            elements,
            nodal_coordinates: reference_nodal_coordinates.into(),
            nodal_velocities: NodalVelocities::zero(),
        }
    }
}

impl<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
    SurfaceElementBlock<'a, C, D, E, F, G, N> for ViscoelasticBlock<D, E, F, G, N>
where
    C: Viscoelastic<'a>,
    F: SurfaceElement<'a, C, G, N> + ViscoelasticFiniteElement<'a, C, G, N>,
{
    fn new(
        constitutive_model_parameters: Parameters<'a>,
        connectivity: Connectivity<E, N>,
        reference_nodal_coordinates: ReferenceNodalCoordinates<D>,
        thickness: Scalar,
    ) -> Self {
        let elements = from_fn(|element| {
            <F>::new(
                constitutive_model_parameters,
                connectivity[element]
                    .iter()
                    .map(|node| reference_nodal_coordinates[*node].copy())
                    .collect(),
                &thickness,
            )
        });
        Self {
            connectivity,
            elements,
            nodal_coordinates: reference_nodal_coordinates.into(),
            nodal_velocities: NodalVelocities::zero(),
        }
    }
}

impl<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
    ViscoelasticFiniteElementBlock<'a, C, D, E, F, G, N> for ViscoelasticBlock<D, E, F, G, N>
where
    C: Viscoelastic<'a>,
    F: ViscoelasticFiniteElement<'a, C, G, N>,
    Self: BasicFiniteElementBlock<'a, C, D, E, F, G, N>,
{
    fn calculate_nodal_forces(&self) -> Result<NodalForces<D>, ConstitutiveError> {
        let mut nodal_forces = NodalForces::zero();
        self.get_elements()
            .iter()
            .zip(self.get_connectivity().iter())
            .try_for_each(|(element, element_connectivity)| {
                element
                    .calculate_nodal_forces(
                        &self.calculate_nodal_coordinates_element(element_connectivity),
                        &self.calculate_nodal_velocities_element(element_connectivity),
                    )?
                    .iter()
                    .zip(element_connectivity.iter())
                    .for_each(|(nodal_force, node)| nodal_forces[*node] += nodal_force);
                Ok(())
            })?;
        Ok(nodal_forces)
    }
    fn calculate_nodal_stiffnesses(&self) -> Result<NodalStiffnesses<D>, ConstitutiveError> {
        let mut nodal_stiffnesses = NodalStiffnesses::zero();
        self.get_elements()
            .iter()
            .zip(self.get_connectivity().iter())
            .try_for_each(|(element, element_connectivity)| {
                element
                    .calculate_nodal_stiffnesses(
                        &self.calculate_nodal_coordinates_element(element_connectivity),
                        &self.calculate_nodal_velocities_element(element_connectivity),
                    )?
                    .iter()
                    .zip(element_connectivity.iter())
                    .for_each(|(object, node_a)| {
                        object.iter().zip(element_connectivity.iter()).for_each(
                            |(nodal_stiffness, node_b)| {
                                nodal_stiffnesses[*node_a][*node_b] += nodal_stiffness
                            },
                        )
                    });
                Ok(())
            })?;
        Ok(nodal_stiffnesses)
    }
    fn calculate_nodal_velocities_element(
        &self,
        element_connectivity: &[usize; N],
    ) -> NodalVelocities<N> {
        let nodal_velocities = self.get_nodal_velocities();
        element_connectivity
            .iter()
            .map(|node| nodal_velocities[*node].copy())
            .collect()
    }
    fn get_nodal_velocities(&self) -> &NodalVelocities<D> {
        &self.nodal_velocities
    }
    fn set_nodal_velocities(&mut self, nodal_velocities: NodalVelocities<D>) {
        self.nodal_velocities = nodal_velocities;
    }
}

impl<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
    ElasticHyperviscousFiniteElementBlock<'a, C, D, E, F, G, N> for ViscoelasticBlock<D, E, F, G, N>
where
    C: ElasticHyperviscous<'a>,
    F: ElasticHyperviscousFiniteElement<'a, C, G, N>,
    Self: ViscoelasticFiniteElementBlock<'a, C, D, E, F, G, N>,
{
    fn calculate_viscous_dissipation(&self) -> Result<Scalar, ConstitutiveError> {
        self.get_elements()
            .iter()
            .zip(self.get_connectivity().iter())
            .map(|(element, element_connectivity)| {
                element.calculate_viscous_dissipation(
                    &self.calculate_nodal_coordinates_element(element_connectivity),
                    &self.calculate_nodal_velocities_element(element_connectivity),
                )
            })
            .sum()
    }
    fn calculate_dissipation_potential(&self) -> Result<Scalar, ConstitutiveError> {
        self.get_elements()
            .iter()
            .zip(self.get_connectivity().iter())
            .map(|(element, element_connectivity)| {
                element.calculate_dissipation_potential(
                    &self.calculate_nodal_coordinates_element(element_connectivity),
                    &self.calculate_nodal_velocities_element(element_connectivity),
                )
            })
            .sum()
    }
}

impl<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
    HyperviscoelasticFiniteElementBlock<'a, C, D, E, F, G, N> for ViscoelasticBlock<D, E, F, G, N>
where
    C: Hyperviscoelastic<'a>,
    F: HyperviscoelasticFiniteElement<'a, C, G, N>,
    Self: ElasticHyperviscousFiniteElementBlock<'a, C, D, E, F, G, N>,
{
    fn calculate_helmholtz_free_energy(&self) -> Result<Scalar, ConstitutiveError> {
        self.get_elements()
            .iter()
            .zip(self.get_connectivity().iter())
            .map(|(element, element_connectivity)| {
                element.calculate_helmholtz_free_energy(
                    &self.calculate_nodal_coordinates_element(element_connectivity),
                )
            })
            .sum()
    }
}
