#[cfg(test)]
mod test;

use super::*;

const G: usize = 1;
const M: usize = 2;
const N: usize = 6;
const O: usize = 3;

pub struct Wedge<C>
{
    constitutive_models: [C; G],
    gradient_vectors: GradientVectors<N>,
    reference_normal: ReferenceNormal
}

impl<'a, C> FiniteElement<'a, C, G, N> for Wedge<C>
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
            reference_normal: Self::calculate_normal(&Self::calculate_midplane_coordinates(&reference_nodal_coordinates))
        }
    }
}

impl<'a, C> LinearElement<'a, C, G, M, N, O> for Wedge<C>
where
    C: Constitutive<'a>
{
    fn calculate_deformation_gradient(&self, nodal_coordinates: &NodalCoordinates<N>) -> DeformationGradient
    {
        self.calculate_deformation_gradient_linear_localization_element(nodal_coordinates)
    }
    fn calculate_deformation_gradient_rate(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> DeformationGradientRate
    {
        self.calculate_deformation_gradient_rate_linear_localization_element(nodal_coordinates, nodal_velocities)
    }
    fn calculate_gradient_vectors(reference_nodal_coordinates: &ReferenceNodalCoordinates<N>) -> GradientVectors<N>
    {
        let reference_nodal_coordinates_midplane = Self::calculate_midplane_coordinates(&reference_nodal_coordinates);
        let reference_dual_basis_vectors = Self::calculate_dual_basis(&reference_nodal_coordinates_midplane);
        let reference_normal = Self::calculate_normal(&reference_nodal_coordinates_midplane);
        let gradient_vectors_midplane = Self::calculate_standard_gradient_operator().iter()
        .map(|standard_gradient_operator_a|
            standard_gradient_operator_a.iter()
            .zip(reference_dual_basis_vectors.iter())
            .map(|(standard_gradient_operator_a_m, dual_reference_basis_vector_m)|
                dual_reference_basis_vector_m*standard_gradient_operator_a_m
            ).sum()
        ).collect::<GradientVectors<O>>();
        let mut gradient_vectors = GradientVectors::zero();
        (0..2).for_each(|bottom_or_top|
            (0..O).zip(gradient_vectors_midplane.iter())
            .for_each(|(a, gradient_vectors_midplane_a)|
                gradient_vectors[a + O*bottom_or_top] = gradient_vectors_midplane_a * 0.5
                    + &reference_normal * (2.0*(bottom_or_top as f64) - 1.0) / 3.0
            )
        );
        gradient_vectors
    }
    fn calculate_standard_gradient_operator() -> StandardGradientOperator<M, O>
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

impl<'a, C> LinearSurfaceElement<'a, C, G, M, N, O> for Wedge<C>
where
    C: Constitutive<'a>
{
    fn get_reference_normal(&self) -> &ReferenceNormal
    {
        &self.reference_normal
    }
}

impl<'a, C> LinearLocalizationElement<'a, C, G, M, N, O> for Wedge<C>
where
    C: Constitutive<'a>
{
    fn calculate_jump(nodal_coordinates: &NodalCoordinates<N>) -> Jump
    {
        (0..O).map(|a|
            nodal_coordinates[a + O].iter()
            .zip(nodal_coordinates[a].iter())
            .map(|(nodal_coordinates_ap3, nodal_coordinates_a)|
                (nodal_coordinates_ap3 - nodal_coordinates_a)
            ).collect()
        ).sum::<Jump>() / 3.0
    }
    fn calculate_midplane_coordinates<const I: usize>(nodal_coordinates: &Coordinates<I, N>) -> Coordinates<I, O>
    {
        (0..O).map(|a|
            nodal_coordinates[a + O].iter()
            .zip(nodal_coordinates[a].iter())
            .map(|(nodal_coordinates_ap3, nodal_coordinates_a)|
                (nodal_coordinates_ap3 + nodal_coordinates_a) * 0.5
            ).collect()
        ).collect()
    }
}

impl<'a, C> ElasticFiniteElement<'a, C, G, N> for Wedge<C>
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

impl<'a, C> ElasticLinearElement<'a, C, G, M, N, O> for Wedge<C>
where
    C: Elastic<'a>
{}

impl<'a, C> HyperelasticFiniteElement<'a, C, G, N> for Wedge<C>
where
    C: Hyperelastic<'a>
{
    fn calculate_helmholtz_free_energy(&self, nodal_coordinates: &NodalCoordinates<N>) -> Scalar
    {
        self.calculate_helmholtz_free_energy_linear_element(nodal_coordinates)
    }
}

impl<'a, C> HyperelasticLinearElement<'a, C, G, M, N, O> for Wedge<C>
where
    C: Hyperelastic<'a>
{}

impl<'a, C> ViscoelasticFiniteElement<'a, C, G, N> for Wedge<C>
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

impl<'a, C> ViscoelasticLinearElement<'a, C, G, M, N, O> for Wedge<C>
where
    C: Viscoelastic<'a>
{}

impl<'a, C> ElasticHyperviscousFiniteElement<'a, C, G, N> for Wedge<C>
where
    C: ElasticHyperviscous<'a>
{
    fn calculate_viscous_dissipation(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> Scalar
    {
        self.calculate_viscous_dissipation_linear_element(nodal_coordinates, nodal_velocities)
    }
    fn calculate_dissipation_potential(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> Scalar
    {
        self.calculate_dissipation_potential_linear_element(nodal_coordinates, nodal_velocities)
    }
}

impl<'a, C> ElasticHyperviscousLinearElement<'a, C, G, M, N, O> for Wedge<C>
where
    C: ElasticHyperviscous<'a>
{}

impl<'a, C> HyperviscoelasticFiniteElement<'a, C, G, N> for Wedge<C>
where
    C: Hyperviscoelastic<'a>
{
    fn calculate_helmholtz_free_energy(&self, nodal_coordinates: &NodalCoordinates<N>) -> Scalar
    {
        self.calculate_helmholtz_free_energy_linear_element(nodal_coordinates)
    }
}

impl<'a, C> HyperviscoelasticLinearElement<'a, C, G, M, N, O> for Wedge<C>
where
    C: Hyperviscoelastic<'a>
{}