#[cfg(test)]
mod test;

pub mod triangle;

use super::*;

pub trait CompositeSurfaceElement<'a, C, const G: usize, const M: usize, const N: usize, const O: usize, const P: usize, const Q: usize>
where
    C: Constitutive<'a>,
    Self: CompositeElement<'a, C, G, M, N, O, P, Q>
{
    fn calculate_bases<const I: usize>(nodal_coordinates: &Coordinates<I, O>) -> Bases<I, P>
    {
        Self::calculate_standard_gradient_operators().iter()
        .map(|standard_gradient_operator|
            standard_gradient_operator.iter()
            .zip(nodal_coordinates.iter())
            .map(|(standard_gradient_operator_a, nodal_coordinate_a)|
                standard_gradient_operator_a.iter()
                .map(|standard_gradient_operator_a_m|
                    nodal_coordinate_a * standard_gradient_operator_a_m
                ).collect()
            ).sum()
        ).collect()
    }
    fn calculate_deformation_gradients_composite_surface_element(&self, nodal_coordinates: &NodalCoordinates<O>) -> DeformationGradients<G>
    {
        let normals = Self::calculate_normals(nodal_coordinates);
        self.get_projected_gradient_vectors().iter()
        .zip(self.get_scaled_reference_normals().iter())
        .map(|(projected_gradient_vectors, scaled_reference_normals)|
            nodal_coordinates.iter()
            .zip(projected_gradient_vectors.iter())
            .map(|(nodal_coordinate, projected_gradient_vector)|
                DeformationGradient::dyad(nodal_coordinate, projected_gradient_vector)
            ).sum::<DeformationGradient>() +
            scaled_reference_normals.iter()
            .zip(normals.iter())
            .map(|(scaled_reference_normal, normal)|
                DeformationGradient::dyad(normal, scaled_reference_normal)
            ).sum::<DeformationGradient>()
        ).collect()
    }
    fn calculate_deformation_gradient_rates_composite_surface_element(&self, nodal_coordinates: &NodalCoordinates<O>, nodal_velocities: &NodalVelocities<O>) -> DeformationGradientRates<G>
    {
        let normal_rates = Self::calculate_normal_rates(nodal_coordinates, nodal_velocities);
        self.get_projected_gradient_vectors().iter()
        .zip(self.get_scaled_reference_normals().iter())
        .map(|(projected_gradient_vectors, scaled_reference_normals)|
            nodal_velocities.iter()
            .zip(projected_gradient_vectors.iter())
            .map(|(nodal_velocity, projected_gradient_vector)|
                DeformationGradientRate::dyad(nodal_velocity, projected_gradient_vector)
            ).sum::<DeformationGradientRate>() +
            scaled_reference_normals.iter()
            .zip(normal_rates.iter())
            .map(|(scaled_reference_normal, normal_rate)|
                DeformationGradientRate::dyad(normal_rate, scaled_reference_normal)
            ).sum::<DeformationGradientRate>()
        ).collect()
    }
    fn calculate_dual_bases<const I: usize>(nodal_coordinates: &Coordinates<I, O>) -> Bases<I, P>
    {
        Self::calculate_bases(nodal_coordinates).iter()
        .map(|basis_vectors|
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
        ).collect()
    }
    fn calculate_normals(nodal_coordinates: &NodalCoordinates<O>) -> Normals<P>
    {
        Self::calculate_bases(nodal_coordinates).iter()
        .map(|basis_vectors|
            basis_vectors[0].cross(&basis_vectors[1]).normalized()
        ).collect()
    }
    fn calculate_normal_gradients(nodal_coordinates: &Coordinates<1, O>) -> NormalGradientss<P, O>
    {
        let identity = TensorRank2::<3, 1, 1>::identity();
        let levi_civita_symbol = levi_civita::<1, 1, 1>();
        let mut normalization: Scalar = 0.0;
        let mut normal_vector = Normal::new([0.0, 0.0, 0.0]);
        Self::calculate_standard_gradient_operators().iter()
        .zip(Self::calculate_bases(nodal_coordinates).iter())
        .map(|(standard_gradient_operator, basis_vectors)|{
            normalization = basis_vectors[0].cross(&basis_vectors[1]).norm();
            normal_vector = basis_vectors[0].cross(&basis_vectors[1])/normalization;
            standard_gradient_operator.iter()
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
        }).collect()
    }
    fn calculate_normal_rates(nodal_coordinates: &NodalCoordinates<O>, nodal_velocities: &NodalVelocities<O>) -> NormalRates<P>
    {
        let identity = TensorRank2::<3, 1, 1>::identity();
        let levi_civita_symbol = levi_civita::<1, 1, 1>();
        let mut normalization: Scalar = 0.0;
        let mut normal_vector = Normal::new([0.0, 0.0, 0.0]);
        Self::calculate_standard_gradient_operators().iter()
        .zip(Self::calculate_bases(nodal_coordinates).iter())
        .map(|(standard_gradient_operator, basis_vectors)|{
            normalization = basis_vectors[0].cross(&basis_vectors[1]).norm();
            normal_vector = basis_vectors[0].cross(&basis_vectors[1])/normalization;
            identity.iter()
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
        }).collect()
    }
    fn calculate_normal_tangents(nodal_coordinates: &Coordinates<1, O>) -> NormalTangentss<P, O>
    {
        let identity = TensorRank2::<3, 1, 1>::identity();
        let levi_civita_symbol = levi_civita::<1, 1, 1>();
        let mut normalization: Scalar = 0.0;
        let mut normal_vector = Normal::new([0.0, 0.0, 0.0]);
        Self::calculate_standard_gradient_operators().iter()
        .zip(Self::calculate_bases(nodal_coordinates).iter()
        .zip(Self::calculate_normal_gradients(nodal_coordinates).iter()))
        .map(|(standard_gradient_operator, (basis_vectors, normal_gradients))|{
            normalization = basis_vectors[0].cross(&basis_vectors[1]).norm();
            normal_vector = basis_vectors[0].cross(&basis_vectors[1])/normalization;
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
        }).collect()
    }
    fn calculate_reference_normals(reference_nodal_coordinates: &ReferenceNodalCoordinates<O>) -> ReferenceNormals<P>
    {
        Self::calculate_dual_bases(reference_nodal_coordinates).iter()
        .map(|dual_basis_vectors|
            dual_basis_vectors[0].cross(&dual_basis_vectors[1]).normalized()
        ).collect()
    }
    fn calculate_scaled_reference_normals(reference_nodal_coordinates: &ReferenceNodalCoordinates<O>) -> ScaledReferenceNormals<G, P>
    {
        let inverse_normalized_projection_matrix = Self::calculate_inverse_normalized_projection_matrix();
        let reference_normals = Self::calculate_reference_normals(reference_nodal_coordinates);
        let shape_function_integrals = Self::calculate_shape_function_integrals();
        Self::calculate_shape_functions_at_integration_points().iter()
        .map(|shape_function|
            reference_normals.iter()
            .zip(shape_function_integrals.iter())
            .map(|(reference_normal, shape_function_integral)|
                reference_normal * (shape_function * (&inverse_normalized_projection_matrix * shape_function_integral))
            ).collect()
        ).collect()
    }
    fn get_scaled_reference_normals(&self) -> &ScaledReferenceNormals<G, P>;
}

macro_rules! composite_surface_element_boilerplate
{
    ($element: ident) =>
    {
        impl<'a, C> CompositeSurfaceElement<'a, C, G, M, N, O, P, Q> for $element<C>
        where
            C: Constitutive<'a>
        {
            fn get_scaled_reference_normals(&self) -> &ScaledReferenceNormals<G, P>
            {
                &self.scaled_reference_normals
            }
        }
        impl<'a, C> ElasticFiniteElement<'a, C, G, N> for $element<C>
        where
            C: Elastic<'a>
        {
            fn calculate_nodal_forces(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalForces<N>
            {
                self.calculate_nodal_forces_composite_element(nodal_coordinates)
            }
            fn calculate_nodal_stiffnesses(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalStiffnesses<N>
            {
                todo!()
            }
        }
        impl<'a, C> ElasticCompositeElement<'a, C, G, M, N, O, P, Q> for $element<C>
        where
            C: Elastic<'a>
        {}
        impl<'a, C> HyperelasticFiniteElement<'a, C, G, N> for $element<C>
        where
            C: Hyperelastic<'a>
        {
            fn calculate_helmholtz_free_energy(&self, nodal_coordinates: &NodalCoordinates<N>) -> Scalar
            {
                self.calculate_helmholtz_free_energy_composite_element(nodal_coordinates)
            }
        }
        impl<'a, C> HyperelasticCompositeElement<'a, C, G, M, N, O, P, Q> for $element<C>
        where
            C: Hyperelastic<'a>
        {}
        impl<'a, C> ViscoelasticFiniteElement<'a, C, G, N> for $element<C>
        where
            C: Viscoelastic<'a>
        {
            fn calculate_nodal_forces(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> NodalForces<N>
            {
                self.calculate_nodal_forces_composite_element(nodal_coordinates, nodal_velocities)
            }
            fn calculate_nodal_stiffnesses(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> NodalStiffnesses<N>
            {
                todo!()
            }
        }
        impl<'a, C> ViscoelasticCompositeElement<'a, C, G, M, N, O, P, Q> for $element<C>
        where
            C: Viscoelastic<'a>
        {}
        impl<'a, C> ElasticHyperviscousFiniteElement<'a, C, G, N> for $element<C>
        where
            C: ElasticHyperviscous<'a>
        {
            fn calculate_viscous_dissipation(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> Scalar
            {
                self.calculate_viscous_dissipation_composite_element(nodal_coordinates, nodal_velocities)
            }
            fn calculate_dissipation_potential(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> Scalar
            {
                self.calculate_dissipation_potential_composite_element(nodal_coordinates, nodal_velocities)
            }
        }
        impl<'a, C> ElasticHyperviscousCompositeElement<'a, C, G, M, N, O, P, Q> for $element<C>
        where
            C: ElasticHyperviscous<'a>
        {}
        impl<'a, C> HyperviscoelasticFiniteElement<'a, C, G, N> for $element<C>
        where
            C: Hyperviscoelastic<'a>
        {
            fn calculate_helmholtz_free_energy(&self, nodal_coordinates: &NodalCoordinates<N>) -> Scalar
            {
                self.calculate_helmholtz_free_energy_composite_element(nodal_coordinates)
            }
        }
        impl<'a, C> HyperviscoelasticCompositeElement<'a, C, G, M, N, O, P, Q> for $element<C>
        where
            C: Hyperviscoelastic<'a>
        {}
    }
}
pub(crate) use composite_surface_element_boilerplate;