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
    fn calculate_objects(&self, normal_gradients: &NormalGradientss<P, O>) -> TensorRank3List2D<3, 1, 1, 0, O, G>
    {
        self.get_scaled_reference_normals().iter()
        .map(|scaled_reference_normals|
            normal_gradients.iter()
            .zip(scaled_reference_normals.iter())
            .map(|(normal_gradient, scaled_reference_normal)|
                normal_gradient.iter()
                .map(|normal_gradient_a|
                    normal_gradient_a.iter()
                    .map(|normal_gradient_a_m|
                        TensorRank2::dyad(normal_gradient_a_m, scaled_reference_normal)
                    ).collect::<TensorRank3<3, 1, 1, 0>>()
                ).collect::<TensorRank3List<3, 1, 1, 0, O>>()
            ).sum::<TensorRank3List<3, 1, 1, 0, O>>()
        ).collect()
    }
    fn calculate_projected_gradient_vectors_composite_surface_element(reference_nodal_coordinates: &ReferenceNodalCoordinates<O>) -> ProjectedGradientVectors<G, N>
    {
        let reference_dual_bases = Self::calculate_dual_bases(reference_nodal_coordinates);
        let reference_jacobians = Self::calculate_reference_jacobians(reference_nodal_coordinates);
        let inverse_projection_matrix = Self::calculate_inverse_projection_matrix(&reference_jacobians);
        Self::calculate_shape_functions_at_integration_points().iter()
        .map(|shape_functions_at_integration_point|
            Self::calculate_standard_gradient_operators_transposed().iter()
            .map(|standard_gradient_operators_a|
                Self::calculate_shape_function_integrals().iter()
                .zip(standard_gradient_operators_a.iter()
                .zip(reference_dual_bases.iter()
                .zip(reference_jacobians.iter())))
                .map(|(shape_function_integral, (standard_gradient_operator, (reference_dual_basis_vectors, reference_jacobian)))|
                    reference_dual_basis_vectors.iter()
                    .zip(standard_gradient_operator.iter())
                    .map(|(reference_dual_basis_vector, standard_gradient_operator_mu)|
                        reference_dual_basis_vector * standard_gradient_operator_mu
                    ).sum::<Vector<0>>() * reference_jacobian * (
                        shape_functions_at_integration_point * (&inverse_projection_matrix * shape_function_integral)
                    )
                ).sum()
            ).collect()
        ).collect()
    }
    fn calculate_reference_normals(reference_nodal_coordinates: &ReferenceNodalCoordinates<O>) -> ReferenceNormals<P>
    {
        Self::calculate_dual_bases(reference_nodal_coordinates).iter()
        .map(|reference_dual_basis_vectors|
            reference_dual_basis_vectors[0].cross(&reference_dual_basis_vectors[1]).normalized()
        ).collect()
    }
    fn calculate_scaled_reference_normals(reference_nodal_coordinates: &ReferenceNodalCoordinates<O>) -> ScaledReferenceNormals<G, P>
    {
        let reference_jacobians = Self::calculate_reference_jacobians(reference_nodal_coordinates);
        let inverse_projection_matrix = Self::calculate_inverse_projection_matrix(&reference_jacobians);
        let reference_normals = Self::calculate_reference_normals(reference_nodal_coordinates);
        let shape_function_integrals = Self::calculate_shape_function_integrals();
        Self::calculate_shape_functions_at_integration_points().iter()
        .map(|shape_function|
            reference_normals.iter()
            .zip(shape_function_integrals.iter()
            .zip(reference_jacobians.iter()))
            .map(|(reference_normal, (shape_function_integral, reference_jacobian))|
                reference_normal * ((shape_function * (&inverse_projection_matrix * shape_function_integral)) * reference_jacobian)
            ).collect()
        ).collect()
    }
    fn get_scaled_reference_normals(&self) -> &ScaledReferenceNormals<G, P>;
}

macro_rules! composite_surface_element_boilerplate
{
    ($element: ident) =>
    {
        impl<'a, C> ElasticFiniteElement<'a, C, G, N> for $element<C>
        where
            C: Elastic<'a>
        {
            fn calculate_nodal_forces(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalForces<N>
            {
                let identity = TensorRank2::<3, 1, 1>::identity();
                self.get_constitutive_models().iter()
                .zip(self.calculate_deformation_gradients(nodal_coordinates).iter())
                .map(|(constitutive_model, deformation_gradient)|
                    constitutive_model.calculate_first_piola_kirchoff_stress(deformation_gradient)
                ).collect::<FirstPiolaKirchoffStresses<G>>().iter()
                .zip(self.get_projected_gradient_vectors().iter()
                .zip(self.get_scaled_composite_jacobians().iter()
                .zip(self.calculate_objects(&Self::calculate_normal_gradients(nodal_coordinates)).iter())))
                .map(|(first_piola_kirchoff_stress, (projected_gradient_vectors, (scaled_composite_jacobian, objects)))|
                    projected_gradient_vectors.iter()
                    .zip(objects.iter())
                    .map(|(projected_gradient_vector, object)|
                        identity.iter()
                        .zip(object.iter())
                        .map(|(identity_m, object_m)|
                            first_piola_kirchoff_stress.iter()
                            .zip(identity_m.iter()
                            .zip(object_m.iter()))
                            .map(|(first_piola_kirchoff_stress_i, (identity_mi, object_mi))|
                                first_piola_kirchoff_stress_i.iter()
                                .zip(projected_gradient_vector.iter()
                                .zip(object_mi.iter()))
                                .map(|(first_piola_kirchoff_stress_ij, (projected_gradient_vector_j, object_mij))|
                                    first_piola_kirchoff_stress_ij * (
                                        identity_mi * projected_gradient_vector_j + object_mij
                                    ) * scaled_composite_jacobian
                                ).sum::<Scalar>()
                            ).sum::<Scalar>()
                        ).collect()
                    ).collect()
                ).sum()
            }
            fn calculate_nodal_stiffnesses(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalStiffnesses<N>
            {
                let deformation_gradients = self.calculate_deformation_gradients(nodal_coordinates);
                let identity = TensorRank2::<3, 1, 1>::identity();
                let normal_tangentss = Self::calculate_normal_tangents(&nodal_coordinates);
                let objectss = self.calculate_objects(&Self::calculate_normal_gradients(&nodal_coordinates));
                let mut scaled_traction = Vector::zero();
                self.get_constitutive_models().iter()
                .zip(deformation_gradients.iter())
                .map(|(constitutive_model, deformation_gradient)|
                    constitutive_model.calculate_first_piola_kirchoff_stress(deformation_gradient)
                ).collect::<FirstPiolaKirchoffStresses<G>>().iter()
                .zip(self.get_constitutive_models().iter()
                .zip(deformation_gradients.iter())
                .map(|(constitutive_model, deformation_gradient)|
                    constitutive_model.calculate_first_piola_kirchoff_tangent_stiffness(&deformation_gradient)
                ).collect::<FirstPiolaKirchoffTangentStiffnesses<G>>().iter()
                .zip(self.get_projected_gradient_vectors().iter()
                .zip(self.get_scaled_composite_jacobians().iter()
                .zip(self.get_scaled_reference_normals().iter()
                .zip(objectss.iter())))))
                .map(|(first_piola_kirchoff_stress, (first_piola_kirchoff_tangent_stiffness, (projected_gradient_vectors, (scaled_composite_jacobian, (scaled_reference_normals, objects)))))|
                    projected_gradient_vectors.iter()
                    .zip(objects.iter())
                    .map(|(projected_gradient_vector_a, object_a)|
                        projected_gradient_vectors.iter()
                        .zip(objects.iter())
                        .map(|(projected_gradient_vector_b, object_b)|
                            identity.iter()
                            .zip(object_a.iter())
                            .map(|(identity_m, object_a_m)|
                                identity.iter()
                                .zip(object_b.iter())
                                .map(|(identity_n, object_b_n)|
                                    first_piola_kirchoff_tangent_stiffness.iter()
                                    .zip(identity_m.iter()
                                    .zip(object_a_m.iter()))
                                    .map(|(first_piola_kirchoff_tangent_stiffness_i, (identity_mi, object_a_mi))|
                                        first_piola_kirchoff_tangent_stiffness_i.iter()
                                        .zip(projected_gradient_vector_a.iter()
                                        .zip(object_a_mi.iter()))
                                        .map(|(first_piola_kirchoff_tangent_stiffness_ij, (projected_gradient_vector_a_j, object_a_mij))|
                                            first_piola_kirchoff_tangent_stiffness_ij.iter()
                                            .zip(identity_n.iter()
                                            .zip(object_b_n.iter()))
                                            .map(|(first_piola_kirchoff_tangent_stiffness_ijk, (identity_nk, object_b_nk))|
                                                first_piola_kirchoff_tangent_stiffness_ijk.iter()
                                                .zip(projected_gradient_vector_b.iter()
                                                .zip(object_b_nk.iter()))
                                                .map(|(first_piola_kirchoff_tangent_stiffness_ijkl, (projected_gradient_vector_b_l, object_b_nkl))|
                                                    first_piola_kirchoff_tangent_stiffness_ijkl * (
                                                        identity_mi * projected_gradient_vector_a_j + object_a_mij
                                                    ) * (
                                                        identity_nk * projected_gradient_vector_b_l + object_b_nkl
                                                    ) * scaled_composite_jacobian
                                                ).sum::<Scalar>()
                                            ).sum::<Scalar>()
                                        ).sum::<Scalar>()
                                    ).sum::<Scalar>()
                                ).collect()
                            ).collect()
                        ).collect()
                    ).collect::<NodalStiffnesses<N>>() +
                    normal_tangentss.iter()
                    .zip(scaled_reference_normals.iter())
                    .map(|(normal_tangents, scaled_reference_normal)|{
                        scaled_traction = (first_piola_kirchoff_stress * scaled_reference_normal) * scaled_composite_jacobian;
                        normal_tangents.iter()
                        .map(|normal_tangent_a|
                            normal_tangent_a.iter()
                            .map(|normal_tangent_ab|
                                normal_tangent_ab.iter()
                                .map(|normal_tangent_ab_m|
                                    normal_tangent_ab_m.iter()
                                    .map(|normal_tangent_ab_mn|
                                        normal_tangent_ab_mn * &scaled_traction
                                    ).collect()
                                ).collect()
                            ).collect()
                        ).collect::<NodalStiffnesses<N>>()
                    }).sum::<NodalStiffnesses<N>>()
                ).sum()
            }
        }
        impl<'a, C> ViscoelasticFiniteElement<'a, C, G, N> for $element<C>
        where
            C: Viscoelastic<'a>
        {
            fn calculate_nodal_forces(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> NodalForces<N>
            {
                let identity = TensorRank2::<3, 1, 1>::identity();
                let normal_gradients = Self::calculate_normal_gradients(nodal_coordinates);
                self.get_constitutive_models().iter()
                .zip(self.calculate_deformation_gradients(nodal_coordinates).iter()
                .zip(self.calculate_deformation_gradient_rates(nodal_coordinates, nodal_velocities).iter()))
                .map(|(constitutive_model, (deformation_gradient, deformation_gradient_rate))|
                    constitutive_model.calculate_first_piola_kirchoff_stress(deformation_gradient, deformation_gradient_rate)
                ).collect::<FirstPiolaKirchoffStresses<G>>().iter()
                .zip(self.get_projected_gradient_vectors().iter()
                .zip(self.get_scaled_composite_jacobians().iter()
                .zip(self.calculate_objects(&normal_gradients).iter())))
                .map(|(first_piola_kirchoff_stress, (projected_gradient_vectors, (scaled_composite_jacobian, objects)))|
                    projected_gradient_vectors.iter()
                    .zip(objects.iter())
                    .map(|(projected_gradient_vector, object)|
                        identity.iter()
                        .zip(object.iter())
                        .map(|(identity_m, object_m)|
                            first_piola_kirchoff_stress.iter()
                            .zip(identity_m.iter()
                            .zip(object_m.iter()))
                            .map(|(first_piola_kirchoff_stress_i, (identity_mi, object_mi))|
                                first_piola_kirchoff_stress_i.iter()
                                .zip(projected_gradient_vector.iter()
                                .zip(object_mi.iter()))
                                .map(|(first_piola_kirchoff_stress_ij, (projected_gradient_vector_j, object_mij))|
                                    first_piola_kirchoff_stress_ij * (
                                        identity_mi * projected_gradient_vector_j + object_mij
                                    ) * scaled_composite_jacobian
                                ).sum::<Scalar>()
                            ).sum::<Scalar>()
                        ).collect()
                    ).collect()
                ).sum()
            }
            fn calculate_nodal_stiffnesses(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> NodalStiffnesses<N>
            {
                let identity = TensorRank2::<3, 1, 1>::identity();
                let normal_gradients = Self::calculate_normal_gradients(&nodal_coordinates);
                let objectss = self.calculate_objects(&normal_gradients);
                self.get_constitutive_models().iter()
                .zip(self.calculate_deformation_gradients(nodal_coordinates).iter()
                .zip(self.calculate_deformation_gradient_rates(nodal_coordinates, nodal_velocities).iter()))
                .map(|(constitutive_model, (deformation_gradient, deformation_gradient_rate))|
                    constitutive_model.calculate_first_piola_kirchoff_rate_tangent_stiffness(deformation_gradient, deformation_gradient_rate)
                ).collect::<FirstPiolaKirchoffRateTangentStiffnesses<G>>().iter()
                .zip(self.get_projected_gradient_vectors().iter()
                .zip(self.get_scaled_composite_jacobians().iter()
                .zip(objectss.iter())))
                .map(|(first_piola_kirchoff_rate_tangent_stiffness, (projected_gradient_vectors, (scaled_composite_jacobian, objects)))|
                    projected_gradient_vectors.iter()
                    .zip(objects.iter())
                    .map(|(projected_gradient_vector_a, object_a)|
                        projected_gradient_vectors.iter()
                        .zip(objects.iter())
                        .map(|(projected_gradient_vector_b, object_b)|
                            identity.iter()
                            .zip(object_a.iter())
                            .map(|(identity_m, object_a_m)|
                                identity.iter()
                                .zip(object_b.iter())
                                .map(|(identity_n, object_b_n)|
                                    first_piola_kirchoff_rate_tangent_stiffness.iter()
                                    .zip(identity_m.iter()
                                    .zip(object_a_m.iter()))
                                    .map(|(first_piola_kirchoff_rate_tangent_stiffness_i, (identity_mi, object_a_mi))|
                                        first_piola_kirchoff_rate_tangent_stiffness_i.iter()
                                        .zip(projected_gradient_vector_a.iter()
                                        .zip(object_a_mi.iter()))
                                        .map(|(first_piola_kirchoff_rate_tangent_stiffness_ij, (projected_gradient_vector_a_j, object_a_mij))|
                                            first_piola_kirchoff_rate_tangent_stiffness_ij.iter()
                                            .zip(identity_n.iter()
                                            .zip(object_b_n.iter()))
                                            .map(|(first_piola_kirchoff_rate_tangent_stiffness_ijk, (identity_nk, object_b_nk))|
                                                first_piola_kirchoff_rate_tangent_stiffness_ijk.iter()
                                                .zip(projected_gradient_vector_b.iter()
                                                .zip(object_b_nk.iter()))
                                                .map(|(first_piola_kirchoff_rate_tangent_stiffness_ijkl, (projected_gradient_vector_b_l, object_b_nkl))|
                                                    first_piola_kirchoff_rate_tangent_stiffness_ijkl * (
                                                        identity_mi * projected_gradient_vector_a_j + object_a_mij
                                                    ) * (
                                                        identity_nk * projected_gradient_vector_b_l + object_b_nkl
                                                    ) * scaled_composite_jacobian
                                                ).sum::<Scalar>()
                                            ).sum::<Scalar>()
                                        ).sum::<Scalar>()
                                    ).sum::<Scalar>()
                                ).collect()
                            ).collect()
                        ).collect()
                    ).collect::<NodalStiffnesses<N>>()
                ).sum()
            }
        }
        composite_surface_or_localization_element_boilerplate!($element);
    }
}
pub(crate) use composite_surface_element_boilerplate;

macro_rules! composite_surface_or_localization_element_boilerplate
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
pub(crate) use composite_surface_or_localization_element_boilerplate;

macro_rules! composite_surface_element_boilerplate_inner
{
    () =>
    {
        fn calculate_inverse_normalized_projection_matrix() -> NormalizedProjectionMatrix<Q>
        {
            let diag: Scalar = 3.0/64.0;
            let off: Scalar = -1.0/64.0;
            NormalizedProjectionMatrix::new([
                [diag,  off,  off],
                [ off, diag,  off],
                [ off,  off, diag]
            ])
        }
        fn calculate_reference_jacobians(reference_nodal_coordinates: &ReferenceNodalCoordinates<O>) -> Scalars<P>
        {
            Self::calculate_bases(reference_nodal_coordinates).iter()
            .map(|reference_basis_vectors|
                reference_basis_vectors[0].cross(&reference_basis_vectors[1]).norm()
            ).collect()
        }
        fn calculate_shape_function_integrals() -> ShapeFunctionIntegrals<P, Q>
        {
            ShapeFunctionIntegrals::new([
                [32.0,  8.0,  8.0],
                [ 8.0, 32.0,  8.0],
                [ 8.0,  8.0, 32.0],
                [16.0, 16.0, 16.0]
            ])
        }
        fn calculate_shape_function_integrals_products() -> ShapeFunctionIntegralsProducts<P, Q>
        {
            ShapeFunctionIntegralsProducts::new([[
                [22.0,  5.0,  5.0],
                [ 5.0,  2.0,  1.0],
                [ 5.0,  1.0,  2.0]
            ], [
                [ 2.0,  5.0,  1.0],
                [ 5.0, 22.0,  5.0],
                [ 1.0,  5.0,  2.0]
            ], [
                [ 2.0,  1.0,  5.0],
                [ 1.0,  2.0,  5.0],
                [ 5.0,  5.0, 22.0]
            ], [
                [ 6.0,  5.0,  5.0],
                [ 5.0,  6.0,  5.0],
                [ 5.0,  5.0,  6.0]
            ]])
        }
        fn calculate_shape_functions_at_integration_points() -> ShapeFunctionsAtIntegrationPoints<G, Q>
        {
            let diag: Scalar = 0.666_666_666_666_666_6;
            let off: Scalar = 0.166_666_666_666_666_7;
            ShapeFunctionsAtIntegrationPoints::new([
                [diag,  off,  off],
                [ off, diag,  off],
                [ off,  off, diag]
            ])
        }
        fn calculate_standard_gradient_operators() -> StandardGradientOperators<M, O, P>
        {
            StandardGradientOperators::new([[
                [ 2.0,  0.0],
                [ 0.0,  0.0],
                [ 0.0,  0.0],
                [ 0.0,  2.0],
                [ 0.0,  0.0],
                [-2.0, -2.0]
            ], [
                [ 0.0,  0.0],
                [ 0.0,  2.0],
                [ 0.0,  0.0],
                [ 2.0,  0.0],
                [-2.0, -2.0],
                [ 0.0,  0.0]
            ], [
                [ 0.0,  0.0],
                [ 0.0,  0.0],
                [-2.0, -2.0],
                [ 0.0,  0.0],
                [ 0.0,  2.0],
                [ 2.0,  0.0]
            ], [
                [ 0.0,  0.0],
                [ 0.0,  0.0],
                [ 0.0,  0.0],
                [ 2.0,  2.0],
                [-2.0,  0.0],
                [ 0.0, -2.0]
            ]])
        }
        fn calculate_standard_gradient_operators_transposed() -> StandardGradientOperatorsTransposed<M, O, P>
        {
            let standard_gradient_operators = Self::calculate_standard_gradient_operators();
            let mut standard_gradient_operators_transposed = StandardGradientOperatorsTransposed::zero();
            standard_gradient_operators_transposed.iter_mut().enumerate()
            .for_each(|(n, standard_gradient_operators_transposed_n)|
                standard_gradient_operators_transposed_n.iter_mut()
                .zip(standard_gradient_operators.iter())
                .for_each(|(standard_gradient_operators_transposed_n_e, standard_gradient_operators_e)|
                    standard_gradient_operators_transposed_n_e.iter_mut()
                    .zip(standard_gradient_operators_e[n].iter())
                    .for_each(|(standard_gradient_operators_transposed_n_e_i, standard_gradient_operators_e_n_i)|
                        *standard_gradient_operators_transposed_n_e_i = *standard_gradient_operators_e_n_i
                    )
                )
            );
            standard_gradient_operators_transposed
        }
        fn get_constitutive_models(&self) -> &[C; G]
        {
            &self.constitutive_models
        }
        fn get_integration_weight() -> Scalar
        {
            INTEGRATION_WEIGHT
        }
        fn get_projected_gradient_vectors(&self) -> &ProjectedGradientVectors<G, N>
        {
            &self.projected_gradient_vectors
        }
        fn get_scaled_composite_jacobians(&self) -> &Scalars<G>
        {
            &self.scaled_composite_jacobians
        }
    }
}
pub(crate) use composite_surface_element_boilerplate_inner;