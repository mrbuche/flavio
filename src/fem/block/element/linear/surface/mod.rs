#[cfg(test)]
pub mod test;

pub mod triangle;

use super::*;

pub trait LinearSurfaceElement<'a, C, const G: usize, const M: usize, const N: usize, const O: usize>
where
    C: Constitutive<'a>,
    Self: LinearElement<'a, C, G, M, N, O>
{
    fn calculate_basis<const I: usize>(nodal_coordinates: &Coordinates<I, O>) -> Basis<I>
    {
        Self::calculate_standard_gradient_operator().iter()
        .zip(nodal_coordinates.iter())
        .map(|(standard_gradient_operator_a, nodal_coordinates_a)|
            standard_gradient_operator_a.iter()
            .map(|standard_gradient_operator_a_m|
                nodal_coordinates_a * standard_gradient_operator_a_m
            ).collect()
        ).sum()
    }
    fn calculate_deformation_gradient_linear_surface_element(&self, nodal_coordinates: &NodalCoordinates<O>) -> DeformationGradient
    {
        nodal_coordinates.iter()
        .zip(self.get_gradient_vectors().iter())
        .map(|(nodal_coordinate, gradient_vector)|
            DeformationGradient::dyad(nodal_coordinate, gradient_vector)
        ).sum::<DeformationGradient>() + DeformationGradient::dyad(
            &Self::calculate_normal(nodal_coordinates),
            self.get_reference_normal()
        )
    }
    fn calculate_deformation_gradient_rate_linear_surface_element(&self, nodal_coordinates: &NodalCoordinates<O>, nodal_velocities: &NodalVelocities<O>) -> DeformationGradientRate
    {
        nodal_velocities.iter()
        .zip(self.get_gradient_vectors().iter())
        .map(|(nodal_velocity, gradient_vector)|
            DeformationGradientRate::dyad(nodal_velocity, gradient_vector)
        ).sum::<DeformationGradientRate>() + DeformationGradientRate::dyad(
            &Self::calculate_normal_rate(nodal_coordinates, nodal_velocities),
            self.get_reference_normal()
        )
    }
    fn calculate_dual_basis<const I: usize>(nodal_coordinates: &Coordinates<I, O>) -> Basis<I>
    {
        let basis_vectors = Self::calculate_basis(nodal_coordinates);
        basis_vectors.iter()
        .map(|basis_vectors_m|
            basis_vectors.iter()
            .map(|basis_vectors_n|
                basis_vectors_m * basis_vectors_n
            ).collect()
        ).collect::<TensorRank2<M, I, I>>()
        .inverse()
        .iter()
        .map(|metric_tensor_m|
            metric_tensor_m.iter()
            .zip(basis_vectors.iter())
            .map(|(metric_tensor_mn, basis_vectors_n)|
                basis_vectors_n * metric_tensor_mn
            ).sum()
        ).collect()
    }
    fn calculate_gradient_vectors_linear_surface_element(reference_nodal_coordinates: &ReferenceNodalCoordinates<O>) -> GradientVectors<N>
    {
        let reference_dual_basis_vectors = Self::calculate_dual_basis(reference_nodal_coordinates);
        Self::calculate_standard_gradient_operator().iter()
        .map(|standard_gradient_operator_a|
            standard_gradient_operator_a.iter()
            .zip(reference_dual_basis_vectors.iter())
            .map(|(standard_gradient_operator_a_m, dual_reference_basis_vector_m)|
                dual_reference_basis_vector_m * standard_gradient_operator_a_m
            ).sum()
        ).collect()
    }
    fn calculate_normal(nodal_coordinates: &NodalCoordinates<O>) -> Normal
    {
        let basis_vectors = Self::calculate_basis(nodal_coordinates);
        basis_vectors[0].cross(&basis_vectors[1]).normalized()
    }
    fn calculate_normal_gradients(nodal_coordinates: &Coordinates<1, O>) -> NormalGradients<O>
    {
        let basis_vectors = Self::calculate_basis(nodal_coordinates);
        let identity = TensorRank2::<3, 1, 1>::identity();
        let levi_civita_symbol = levi_civita::<1, 1, 1>();
        let normalization = basis_vectors[0].cross(&basis_vectors[1]).norm();
        let normal_vector = basis_vectors[0].cross(&basis_vectors[1])/normalization;
        Self::calculate_standard_gradient_operator().iter()
        .map(|standard_gradient_operator_a|
            levi_civita_symbol.iter()
            .map(|levi_civita_symbol_m|
                identity.iter()
                .zip(normal_vector.iter())
                .map(|(identity_i, normal_vector_i)|
                    levi_civita_symbol_m.iter()
                    .zip(basis_vectors[0].iter()
                    .zip(basis_vectors[1].iter()))
                    .map(|(levi_civita_symbol_mn, (basis_vector_0_n, basis_vector_1_n))|
                        levi_civita_symbol_mn.iter()
                        .zip(identity_i.iter()
                        .zip(normal_vector.iter()))
                        .map(|(levi_civita_symbol_mno, (identity_io, normal_vector_o))|
                            levi_civita_symbol_mno * (identity_io - normal_vector_i * normal_vector_o)
                        ).sum::<Scalar>() * (
                            standard_gradient_operator_a[0] * basis_vector_1_n
                          - standard_gradient_operator_a[1] * basis_vector_0_n
                        )
                    ).sum::<Scalar>() / normalization
                ).collect()
            ).collect()
        ).collect()
    }
    fn calculate_normal_rate(nodal_coordinates: &NodalCoordinates<O>, nodal_velocities: &NodalVelocities<O>) -> NormalRate
    {
        let basis_vectors = Self::calculate_basis(nodal_coordinates);
        let levi_civita_symbol = levi_civita::<1, 1, 1>();
        let normalization = basis_vectors[0].cross(&basis_vectors[1]).norm();
        let normal_vector = basis_vectors[0].cross(&basis_vectors[1])/normalization;
        let standard_gradient_operator = Self::calculate_standard_gradient_operator();
        TensorRank2::<3, 1, 1>::identity().iter()
        .zip(normal_vector.iter())
        .map(|(identity_i, normal_vector_i)|
            nodal_velocities.iter()
            .zip(standard_gradient_operator.iter())
            .map(|(nodal_velocity_a, standard_gradient_operator_a)|
                levi_civita_symbol.iter()
                .zip(nodal_velocity_a.iter())
                .map(|(levi_civita_symbol_m, nodal_velocity_a_m)|
                    levi_civita_symbol_m.iter()
                    .zip(basis_vectors[0].iter()
                    .zip(basis_vectors[1].iter()))
                    .map(|(levi_civita_symbol_mn, (basis_vector_0_n, basis_vector_1_n))|
                        levi_civita_symbol_mn.iter()
                        .zip(identity_i.iter()
                        .zip(normal_vector.iter()))
                        .map(|(levi_civita_symbol_mno, (identity_io, normal_vector_o))|
                            levi_civita_symbol_mno * (identity_io - normal_vector_i * normal_vector_o)
                        ).sum::<Scalar>() * (
                            standard_gradient_operator_a[0] * basis_vector_1_n
                          - standard_gradient_operator_a[1] * basis_vector_0_n
                        )
                    ).sum::<Scalar>() * nodal_velocity_a_m
                ).sum::<Scalar>()
            ).sum::<Scalar>() / normalization
        ).collect()
    }
    fn calculate_normal_tangents(nodal_coordinates: &Coordinates<1, O>) -> NormalTangents<O>
    {
        let basis_vectors = Self::calculate_basis(nodal_coordinates);
        let identity = TensorRank2::<3, 1, 1>::identity();
        let levi_civita_symbol = levi_civita::<1, 1, 1>();
        let normalization = basis_vectors[0].cross(&basis_vectors[1]).norm();
        let normal_gradients = Self::calculate_normal_gradients(nodal_coordinates);
        let normal_vector = basis_vectors[0].cross(&basis_vectors[1])/normalization;
        let standard_gradient_operator = Self::calculate_standard_gradient_operator();
        normal_gradients.iter()
        .zip(standard_gradient_operator.iter())
        .map(|(normal_gradient_a, standard_gradient_operator_a)|
            normal_gradients.iter()
            .zip(standard_gradient_operator.iter())
            .map(|(normal_gradient_b, standard_gradient_operator_b)|
                normal_gradient_a.iter()
                .zip(levi_civita_symbol.iter())
                .map(|(normal_gradient_a_m, levi_civita_symbol_m)|
                    normal_gradient_b.iter()
                    .zip(levi_civita_symbol.iter()
                    .zip(levi_civita_symbol_m.iter()))
                    .map(|(normal_gradient_b_n, (levi_civita_symbol_n, levi_civita_symbol_mn))|
                        normal_gradient_a_m.iter()
                        .zip(normal_gradient_b_n.iter()
                        .zip(normal_vector.iter()
                        .zip(identity.iter())))
                        .map(|(normal_gradient_a_m_i, (normal_gradient_b_n_i, (normal_vector_i, identity_i)))|
                            (levi_civita_symbol_m.iter()
                            .zip(levi_civita_symbol_n.iter()
                            .zip(basis_vectors[0].iter()
                            .zip(basis_vectors[1].iter())))
                            .map(|(levi_civita_symbol_mr, (levi_civita_symbol_nr, (basis_vector_0_r, basis_vector_1_r)))|
                                levi_civita_symbol_mr.iter()
                                .zip(normal_vector.iter()
                                .zip(normal_gradient_b_n.iter()))
                                .map(|(levi_civita_symbol_mrs, (normal_vector_s, normal_gradient_b_n_s))|
                                    levi_civita_symbol_mrs * (
                                        normal_gradient_b_n_i * normal_vector_s
                                      + normal_gradient_b_n_s * normal_vector_i
                                    )
                                ).sum::<Scalar>() * (
                                    standard_gradient_operator_a[1] * basis_vector_0_r
                                  - standard_gradient_operator_a[0] * basis_vector_1_r
                                ) +
                                levi_civita_symbol_nr.iter()
                                .zip(normal_vector.iter())
                                .map(|(levi_civita_symbol_nrs, normal_vector_s)|
                                    levi_civita_symbol_nrs * normal_vector_s * normal_gradient_a_m_i
                                ).sum::<Scalar>() * (
                                    standard_gradient_operator_b[1] * basis_vector_0_r
                                  - standard_gradient_operator_b[0] * basis_vector_1_r
                                )
                            ).sum::<Scalar>() +
                            levi_civita_symbol_mn * (
                                identity_i - &normal_vector * normal_vector_i
                            ) * (
                                standard_gradient_operator_a[0] * standard_gradient_operator_b[1]
                              - standard_gradient_operator_a[1] * standard_gradient_operator_b[0]
                            )) / normalization
                        ).collect()
                    ).collect()
                ).collect()
            ).collect()
        ).collect()
    }
    fn calculate_reference_normal(reference_nodal_coordinates: &ReferenceNodalCoordinates<O>) -> ReferenceNormal
    {
        let dual_basis_vectors = Self::calculate_dual_basis(reference_nodal_coordinates);
        dual_basis_vectors[0].cross(&dual_basis_vectors[1]).normalized()
    }
    fn get_reference_normal(&self) -> &ReferenceNormal;
}

macro_rules! linear_surface_element_boilerplate
{
    ($element: ident) =>
    {
        impl<'a, C> ElasticFiniteElement<'a, C, G, N> for $element<C>
        where
            C: Elastic<'a>
        {
            fn calculate_nodal_forces(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalForces<N>
            {
                self.calculate_nodal_forces_linear_element(nodal_coordinates)
            }
            fn calculate_nodal_stiffnesses(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalStiffnesses<N>
            {
                let first_piola_kirchoff_tangent_stiffness = self.get_constitutive_model()
                .calculate_first_piola_kirchoff_tangent_stiffness(
                    &self.calculate_deformation_gradient(nodal_coordinates)
                );
                let gradient_vectors = self.get_gradient_vectors();
                let identity = TensorRank2::<3, 1, 1>::identity();
                let normal_gradients = Self::calculate_normal_gradients(nodal_coordinates);
                let reference_normal = self.get_reference_normal();
                gradient_vectors.iter()
                .map(|gradient_vector_a|
                    gradient_vectors.iter()
                    .zip(normal_gradients.iter())
                    .map(|(gradient_vector_b, normal_gradient_b)|
                        first_piola_kirchoff_tangent_stiffness.iter()
                        .map(|first_piola_kirchoff_tangent_stiffness_m|
                            identity.iter()
                            .zip(normal_gradient_b.iter())
                            .map(|(identity_n, normal_gradient_b_n)|
                                first_piola_kirchoff_tangent_stiffness_m.iter()
                                .zip(gradient_vector_a.iter())
                                .map(|(first_piola_kirchoff_tangent_stiffness_mj, gradient_vector_a_j)|
                                    first_piola_kirchoff_tangent_stiffness_mj.iter()
                                    .zip(identity_n.iter()
                                    .zip(normal_gradient_b_n.iter()))
                                    .map(|(first_piola_kirchoff_tangent_stiffness_mjk, (identity_nk, normal_gradient_b_n_k))|
                                        first_piola_kirchoff_tangent_stiffness_mjk.iter()
                                        .zip(gradient_vector_b.iter()
                                        .zip(reference_normal.iter()))
                                        .map(|(first_piola_kirchoff_tangent_stiffness_mjkl, (gradient_vector_b_l, reference_normal_l))|
                                            first_piola_kirchoff_tangent_stiffness_mjkl * gradient_vector_a_j * (
                                                identity_nk * gradient_vector_b_l + normal_gradient_b_n_k * reference_normal_l
                                            ) * self.get_integration_weight()
                                        ).sum::<Scalar>()
                                    ).sum::<Scalar>()
                                ).sum::<Scalar>()
                            ).collect()
                        ).collect()
                    ).collect()
                ).collect()
            }
        }
        impl<'a, C> ViscoelasticFiniteElement<'a, C, G, N> for $element<C>
        where
            C: Viscoelastic<'a>
        {
            fn calculate_nodal_forces(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> NodalForces<N>
            {
                self.calculate_nodal_forces_linear_element(nodal_coordinates, nodal_velocities)
            }
            fn calculate_nodal_stiffnesses(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> NodalStiffnesses<N>
            {
                let first_piola_kirchoff_rate_tangent_stiffness = self.get_constitutive_model()
                .calculate_first_piola_kirchoff_rate_tangent_stiffness(
                    &self.calculate_deformation_gradient(nodal_coordinates),
                    &self.calculate_deformation_gradient_rate(nodal_coordinates, nodal_velocities)
                );
                let gradient_vectors = self.get_gradient_vectors();
                let identity = TensorRank2::<3, 1, 1>::identity();
                let normal_gradients = Self::calculate_normal_gradients(nodal_coordinates);
                let reference_normal = self.get_reference_normal();
                gradient_vectors.iter()
                .map(|gradient_vector_a|
                    gradient_vectors.iter()
                    .zip(normal_gradients.iter())
                    .map(|(gradient_vector_b, normal_gradient_b)|
                        first_piola_kirchoff_rate_tangent_stiffness.iter()
                        .map(|first_piola_kirchoff_rate_tangent_stiffness_m|
                            identity.iter()
                            .zip(normal_gradient_b.iter())
                            .map(|(identity_n, normal_gradient_b_n)|
                                first_piola_kirchoff_rate_tangent_stiffness_m.iter()
                                .zip(gradient_vector_a.iter())
                                .map(|(first_piola_kirchoff_rate_tangent_stiffness_mj, gradient_vector_a_j)|
                                    first_piola_kirchoff_rate_tangent_stiffness_mj.iter()
                                    .zip(identity_n.iter()
                                    .zip(normal_gradient_b_n.iter()))
                                    .map(|(first_piola_kirchoff_rate_tangent_stiffness_mjk, (identity_nk, normal_gradient_b_n_k))|
                                        first_piola_kirchoff_rate_tangent_stiffness_mjk.iter()
                                        .zip(gradient_vector_b.iter()
                                        .zip(reference_normal.iter()))
                                        .map(|(first_piola_kirchoff_rate_tangent_stiffness_mjkl, (gradient_vector_b_l, reference_normal_l))|
                                            first_piola_kirchoff_rate_tangent_stiffness_mjkl * gradient_vector_a_j * (
                                                identity_nk * gradient_vector_b_l + normal_gradient_b_n_k * reference_normal_l
                                            ) * self.get_integration_weight()
                                        ).sum::<Scalar>()
                                    ).sum::<Scalar>()
                                ).sum::<Scalar>()
                            ).collect()
                        ).collect()
                    ).collect()
                ).collect()
            }
        }
        linear_surface_or_localization_element_boilerplate!($element);
    }
}
pub(crate) use linear_surface_element_boilerplate;

macro_rules! linear_surface_or_localization_element_boilerplate
{
    ($element: ident) =>
    {
        impl<'a, C> LinearSurfaceElement<'a, C, G, M, N, O> for $element<C>
        where
            C: Constitutive<'a>
        {
            fn get_reference_normal(&self) -> &ReferenceNormal
            {
                &self.reference_normal
            }
        }
        impl<'a, C> ElasticLinearElement<'a, C, G, M, N, O> for $element<C>
        where
            C: Elastic<'a>
        {}
        impl<'a, C> HyperelasticFiniteElement<'a, C, G, N> for $element<C>
        where
            C: Hyperelastic<'a>
        {
            fn calculate_helmholtz_free_energy(&self, nodal_coordinates: &NodalCoordinates<N>) -> Scalar
            {
                self.calculate_helmholtz_free_energy_linear_element(nodal_coordinates)
            }
        }
        impl<'a, C> HyperelasticLinearElement<'a, C, G, M, N, O> for $element<C>
        where
            C: Hyperelastic<'a>
        {}
        impl<'a, C> ViscoelasticLinearElement<'a, C, G, M, N, O> for $element<C>
        where
            C: Viscoelastic<'a>
        {}
        impl<'a, C> ElasticHyperviscousFiniteElement<'a, C, G, N> for $element<C>
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
        impl<'a, C> ElasticHyperviscousLinearElement<'a, C, G, M, N, O> for $element<C>
        where
            C: ElasticHyperviscous<'a>
        {}
        
        impl<'a, C> HyperviscoelasticFiniteElement<'a, C, G, N> for $element<C>
        where
            C: Hyperviscoelastic<'a>
        {
            fn calculate_helmholtz_free_energy(&self, nodal_coordinates: &NodalCoordinates<N>) -> Scalar
            {
                self.calculate_helmholtz_free_energy_linear_element(nodal_coordinates)
            }
        }
        impl<'a, C> HyperviscoelasticLinearElement<'a, C, G, M, N, O> for $element<C>
        where
            C: Hyperviscoelastic<'a>
        {}
    }
}
pub(crate) use linear_surface_or_localization_element_boilerplate;

macro_rules! linear_surface_element_boilerplate_inner
{
    () =>
    {
        fn calculate_standard_gradient_operator() -> StandardGradientOperator<M, O>
        {
            StandardGradientOperator::new([
                [-1.0, -1.0],
                [ 1.0,  0.0],
                [ 0.0,  1.0]
            ])
        }
        fn get_constitutive_model(&self) -> &C
        {
            &self.constitutive_model
        }
        fn get_gradient_vectors(&self) -> &GradientVectors<N>
        {
            &self.gradient_vectors
        }
        fn get_integration_weight(&self) -> &Scalar
        {
            &INTEGRATION_WEIGHT
        }
    }
}
pub(crate) use linear_surface_element_boilerplate_inner;