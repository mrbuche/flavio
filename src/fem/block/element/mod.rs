#[cfg(test)]
mod test;

pub mod composite;
pub mod linear;

use super::*;
use std::fmt;

#[derive(Debug)]
pub enum FiniteElementError {
    Custom(String, DeformationGradient, String),
    InvalidJacobianElastic(Scalar, DeformationGradient, String),
    InvalidJacobianThermoelastic(Scalar, DeformationGradient, Scalar, String),
}

impl From<ConstitutiveError> for FiniteElementError {
    fn from(constitutive_error: ConstitutiveError) -> Self {
        match constitutive_error {
            ConstitutiveError::Custom(message, deformation_gradient, constitutive_model) => {
                Self::Custom(message, deformation_gradient, constitutive_model)
            }
            ConstitutiveError::InvalidJacobianElastic(
                jacobian,
                deformation_gradient,
                constitutive_model,
            ) => Self::InvalidJacobianElastic(jacobian, deformation_gradient, constitutive_model),
            ConstitutiveError::InvalidJacobianThermoelastic(
                jacobian,
                deformation_gradient,
                temperature,
                constitutive_model,
            ) => Self::InvalidJacobianThermoelastic(
                jacobian,
                deformation_gradient,
                temperature,
                constitutive_model,
            ),
        }
    }
}

impl fmt::Display for FiniteElementError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let error = match self {
            Self::Custom(message, deformation_gradient, constitutive_model) => format!(
                "{}\n\
                 From deformation gradient: {}.\n\
                 In constitutive model: {}.",
                message, deformation_gradient, constitutive_model
            ),
            Self::InvalidJacobianElastic(jacobian, deformation_gradient, constitutive_model) => {
                format!(
                    "Invalid Jacobian: {:.6e}.\n\
                 From deformation gradient: {}.\n\
                 In constitutive model: {}.",
                    jacobian, deformation_gradient, constitutive_model
                )
            }
            Self::InvalidJacobianThermoelastic(
                jacobian,
                deformation_gradient,
                temperature,
                constitutive_model,
            ) => format!(
                "Invalid Jacobian: {:.6e}.\n\
                 From deformation gradient: {}.\n\
                 For temperature: {:.6e}.\n\
                 In constitutive model: {}.",
                jacobian, deformation_gradient, temperature, constitutive_model
            ),
        };
        write!(
            f,
            "\x1b[91m{}\n\x1b[0;2;31m{}\x1b[0m\n",
            error,
            get_error_message()
        )
    }
}

pub trait FiniteElement<'a, C, const G: usize, const N: usize>
where
    C: Constitutive<'a>,
{
    fn new(
        constitutive_model_parameters: Parameters<'a>,
        reference_nodal_coordinates: ReferenceNodalCoordinates<N>,
    ) -> Self;
}

pub trait SurfaceElement<'a, C, const G: usize, const N: usize>
where
    C: Constitutive<'a>,
{
    fn new(
        constitutive_model_parameters: Parameters<'a>,
        reference_nodal_coordinates: ReferenceNodalCoordinates<N>,
        thickness: &Scalar,
    ) -> Self;
}

pub trait CohesiveElement<'a, C, const G: usize, const N: usize>
where
    C: Cohesive<'a>,
    Self: FiniteElement<'a, C, G, N>,
{
    fn calculate_nodal_forces(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalForces<N>;
    fn calculate_nodal_stiffnesses(
        &self,
        nodal_coordinates: &NodalCoordinates<N>,
    ) -> NodalStiffnesses<N>;
}

pub trait ElasticFiniteElement<'a, C, const G: usize, const N: usize>
where
    C: Elastic<'a>,
{
    fn calculate_nodal_forces(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalForces<N>;
    fn calculate_nodal_stiffnesses(
        &self,
        nodal_coordinates: &NodalCoordinates<N>,
    ) -> NodalStiffnesses<N>;
}

pub trait HyperelasticFiniteElement<'a, C, const G: usize, const N: usize>
where
    C: Hyperelastic<'a>,
    Self: ElasticFiniteElement<'a, C, G, N>,
{
    fn calculate_helmholtz_free_energy(
        &self,
        nodal_coordinates: &NodalCoordinates<N>,
    ) -> Result<Scalar, FiniteElementError>;
}

pub trait ViscoelasticFiniteElement<'a, C, const G: usize, const N: usize>
where
    C: Viscoelastic<'a>,
{
    fn calculate_nodal_forces(
        &self,
        nodal_coordinates: &NodalCoordinates<N>,
        nodal_velocities: &NodalVelocities<N>,
    ) -> NodalForces<N>;
    fn calculate_nodal_stiffnesses(
        &self,
        nodal_coordinates: &NodalCoordinates<N>,
        nodal_velocities: &NodalVelocities<N>,
    ) -> NodalStiffnesses<N>;
}

pub trait ElasticHyperviscousFiniteElement<'a, C, const G: usize, const N: usize>
where
    C: ElasticHyperviscous<'a>,
    Self: ViscoelasticFiniteElement<'a, C, G, N>,
{
    fn calculate_viscous_dissipation(
        &self,
        nodal_coordinates: &NodalCoordinates<N>,
        nodal_velocities: &NodalVelocities<N>,
    ) -> Scalar;
    fn calculate_dissipation_potential(
        &self,
        nodal_coordinates: &NodalCoordinates<N>,
        nodal_velocities: &NodalVelocities<N>,
    ) -> Scalar;
}

pub trait HyperviscoelasticFiniteElement<'a, C, const G: usize, const N: usize>
where
    C: Hyperviscoelastic<'a>,
    Self: ElasticHyperviscousFiniteElement<'a, C, G, N>,
{
    fn calculate_helmholtz_free_energy(
        &self,
        nodal_coordinates: &NodalCoordinates<N>,
    ) -> Result<Scalar, FiniteElementError>;
}
