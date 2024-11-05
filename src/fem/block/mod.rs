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
}

pub struct ViscoelasticBlock<const D: usize, const E: usize, F, const G: usize, const N: usize> {
    connectivity: Connectivity<E, N>,
    elements: [F; E],
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
        nodal_coordinates: &NodalCoordinates<D>,
    ) -> NodalCoordinates<N>;
    fn get_connectivity(&self) -> &Connectivity<E, N>;
    fn get_elements(&self) -> &[F; E];
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
    fn calculate_nodal_forces(
        &self,
        nodal_coordinates: &NodalCoordinates<D>,
    ) -> Result<NodalForces<D>, ConstitutiveError>;
    fn calculate_nodal_stiffnesses(
        &self,
        nodal_coordinates: &NodalCoordinates<D>,
    ) -> Result<NodalStiffnesses<D>, ConstitutiveError>;
    fn solve(
        &mut self,
        fixed_nodes: &[usize],
        initial_coordinates: &ReferenceNodalCoordinates<D>,
    ) -> Result<(), ConstitutiveError>;
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
    fn calculate_helmholtz_free_energy(
        &self,
        nodal_coordinates: &NodalCoordinates<D>,
    ) -> Result<Scalar, ConstitutiveError>;
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
    fn calculate_nodal_forces(
        &self,
        nodal_coordinates: &NodalCoordinates<D>,
        nodal_velocities: &NodalVelocities<D>,
    ) -> Result<NodalForces<D>, ConstitutiveError>;
    fn calculate_nodal_stiffnesses(
        &self,
        nodal_coordinates: &NodalCoordinates<D>,
        nodal_velocities: &NodalVelocities<D>,
    ) -> Result<NodalStiffnesses<D>, ConstitutiveError>;
    fn calculate_nodal_velocities_element(
        &self,
        element_connectivity: &[usize; N],
        nodal_velocities: &NodalVelocities<D>,
    ) -> NodalVelocities<N>;
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
    fn calculate_viscous_dissipation(
        &self,
        nodal_coordinates: &NodalCoordinates<D>,
        nodal_velocities: &NodalVelocities<D>,
    ) -> Result<Scalar, ConstitutiveError>;
    fn calculate_dissipation_potential(
        &self,
        nodal_coordinates: &NodalCoordinates<D>,
        nodal_velocities: &NodalVelocities<D>,
    ) -> Result<Scalar, ConstitutiveError>;
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
    fn calculate_helmholtz_free_energy(
        &self,
        nodal_coordinates: &NodalCoordinates<D>,
    ) -> Result<Scalar, ConstitutiveError>;
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
        nodal_coordinates: &NodalCoordinates<D>,
    ) -> NodalCoordinates<N> {
        element_connectivity
            .iter()
            .map(|node| nodal_coordinates[*node].copy())
            .collect()
    }
    fn get_connectivity(&self) -> &Connectivity<E, N> {
        &self.connectivity
    }
    fn get_elements(&self) -> &[F; E] {
        &self.elements
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
    fn calculate_nodal_forces(
        &self,
        nodal_coordinates: &NodalCoordinates<D>,
    ) -> Result<NodalForces<D>, ConstitutiveError> {
        let mut nodal_forces = NodalForces::zero();
        self.get_elements()
            .iter()
            .zip(self.get_connectivity().iter())
            .try_for_each(|(element, element_connectivity)| {
                element
                    .calculate_nodal_forces(&self.calculate_nodal_coordinates_element(
                        element_connectivity,
                        nodal_coordinates,
                    ))?
                    .iter()
                    .zip(element_connectivity.iter())
                    .for_each(|(nodal_force, node)| nodal_forces[*node] += nodal_force);
                Ok(())
            })?;
        Ok(nodal_forces)
    }
    fn calculate_nodal_stiffnesses(
        &self,
        nodal_coordinates: &NodalCoordinates<D>,
    ) -> Result<NodalStiffnesses<D>, ConstitutiveError> {
        let mut nodal_stiffnesses = NodalStiffnesses::zero();
        self.get_elements()
            .iter()
            .zip(self.get_connectivity().iter())
            .try_for_each(|(element, element_connectivity)| {
                element
                    .calculate_nodal_stiffnesses(&self.calculate_nodal_coordinates_element(
                        element_connectivity,
                        nodal_coordinates,
                    ))?
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
    fn solve(
        &mut self,
        fixed_nodes: &[usize],
        initial_coordinates: &ReferenceNodalCoordinates<D>,
    ) -> Result<(), ConstitutiveError> {
        //
        // prescribed BCs
        // &[usize] for node ids
        // &[Option<f64>; 3] for values (None = free DOF)
        //
        let mut guess = initial_coordinates.copy().into();
        let mut guess_0 = NodalCoordinates::zero();
        let mut residual;
        let mut residual_0 = NodalForces::zero();
        let mut residual_difference;
        let mut residual_norm: Scalar;
        let mut step_size: Scalar;
        for step in 0..1000 {
            residual = self.calculate_nodal_forces(&guess)?;
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
                guess -= residual * step_size;
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
    fn calculate_helmholtz_free_energy(
        &self,
        nodal_coordinates: &NodalCoordinates<D>,
    ) -> Result<Scalar, ConstitutiveError> {
        self.get_elements()
            .iter()
            .zip(self.get_connectivity().iter())
            .map(|(element, element_connectivity)| {
                element.calculate_helmholtz_free_energy(
                    &self.calculate_nodal_coordinates_element(
                        element_connectivity,
                        nodal_coordinates,
                    ),
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
        nodal_coordinates: &NodalCoordinates<D>,
    ) -> NodalCoordinates<N> {
        element_connectivity
            .iter()
            .map(|node| nodal_coordinates[*node].copy())
            .collect()
    }
    fn get_connectivity(&self) -> &Connectivity<E, N> {
        &self.connectivity
    }
    fn get_elements(&self) -> &[F; E] {
        &self.elements
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
    fn calculate_nodal_forces(
        &self,
        nodal_coordinates: &NodalCoordinates<D>,
        nodal_velocities: &NodalVelocities<D>,
    ) -> Result<NodalForces<D>, ConstitutiveError> {
        let mut nodal_forces = NodalForces::zero();
        self.get_elements()
            .iter()
            .zip(self.get_connectivity().iter())
            .try_for_each(|(element, element_connectivity)| {
                element
                    .calculate_nodal_forces(
                        &self.calculate_nodal_coordinates_element(
                            element_connectivity,
                            nodal_coordinates,
                        ),
                        &self.calculate_nodal_velocities_element(
                            element_connectivity,
                            nodal_velocities,
                        ),
                    )?
                    .iter()
                    .zip(element_connectivity.iter())
                    .for_each(|(nodal_force, node)| nodal_forces[*node] += nodal_force);
                Ok(())
            })?;
        Ok(nodal_forces)
    }
    fn calculate_nodal_stiffnesses(
        &self,
        nodal_coordinates: &NodalCoordinates<D>,
        nodal_velocities: &NodalVelocities<D>,
    ) -> Result<NodalStiffnesses<D>, ConstitutiveError> {
        let mut nodal_stiffnesses = NodalStiffnesses::zero();
        self.get_elements()
            .iter()
            .zip(self.get_connectivity().iter())
            .try_for_each(|(element, element_connectivity)| {
                element
                    .calculate_nodal_stiffnesses(
                        &self.calculate_nodal_coordinates_element(
                            element_connectivity,
                            nodal_coordinates,
                        ),
                        &self.calculate_nodal_velocities_element(
                            element_connectivity,
                            nodal_velocities,
                        ),
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
        nodal_velocities: &NodalVelocities<D>,
    ) -> NodalVelocities<N> {
        element_connectivity
            .iter()
            .map(|node| nodal_velocities[*node].copy())
            .collect()
    }
}

impl<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
    ElasticHyperviscousFiniteElementBlock<'a, C, D, E, F, G, N> for ViscoelasticBlock<D, E, F, G, N>
where
    C: ElasticHyperviscous<'a>,
    F: ElasticHyperviscousFiniteElement<'a, C, G, N>,
    Self: ViscoelasticFiniteElementBlock<'a, C, D, E, F, G, N>,
{
    fn calculate_viscous_dissipation(
        &self,
        nodal_coordinates: &NodalCoordinates<D>,
        nodal_velocities: &NodalVelocities<D>,
    ) -> Result<Scalar, ConstitutiveError> {
        self.get_elements()
            .iter()
            .zip(self.get_connectivity().iter())
            .map(|(element, element_connectivity)| {
                element.calculate_viscous_dissipation(
                    &self.calculate_nodal_coordinates_element(
                        element_connectivity,
                        nodal_coordinates,
                    ),
                    &self
                        .calculate_nodal_velocities_element(element_connectivity, nodal_velocities),
                )
            })
            .sum()
    }
    fn calculate_dissipation_potential(
        &self,
        nodal_coordinates: &NodalCoordinates<D>,
        nodal_velocities: &NodalVelocities<D>,
    ) -> Result<Scalar, ConstitutiveError> {
        self.get_elements()
            .iter()
            .zip(self.get_connectivity().iter())
            .map(|(element, element_connectivity)| {
                element.calculate_dissipation_potential(
                    &self.calculate_nodal_coordinates_element(
                        element_connectivity,
                        nodal_coordinates,
                    ),
                    &self
                        .calculate_nodal_velocities_element(element_connectivity, nodal_velocities),
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
    fn calculate_helmholtz_free_energy(
        &self,
        nodal_coordinates: &NodalCoordinates<D>,
    ) -> Result<Scalar, ConstitutiveError> {
        self.get_elements()
            .iter()
            .zip(self.get_connectivity().iter())
            .map(|(element, element_connectivity)| {
                element.calculate_helmholtz_free_energy(
                    &self.calculate_nodal_coordinates_element(
                        element_connectivity,
                        nodal_coordinates,
                    ),
                )
            })
            .sum()
    }
}
