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
                nodal_coordinates_a.iter()
                .map(|nodal_coordinates_a_i|
                    nodal_coordinates_a_i*standard_gradient_operator_a_m
                ).collect()
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
                basis_vectors_m*basis_vectors_n
            ).collect()
        ).collect::<TensorRank2<M, I, I>>()
        .inverse()
        .iter()
        .map(|metric_tensor_m|
            metric_tensor_m.iter()
            .zip(basis_vectors.iter())
            .map(|(metric_tensor_mn, basis_vectors_n)|
                basis_vectors_n*metric_tensor_mn
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
                dual_reference_basis_vector_m*standard_gradient_operator_a_m
            ).sum()
        ).collect()
    }
    fn calculate_normal<const I: usize>(nodal_coordinates: &Coordinates<I, O>) -> Normal<I>
    {
        let basis_vectors = Self::calculate_basis(nodal_coordinates);
        basis_vectors[0].cross(&basis_vectors[1]).normalized()
    }
    fn calculate_normal_gradients(nodal_coordinates: &Coordinates<1, O>) -> NormalGradients<O>
    {
        let basis_vectors = Self::calculate_basis(nodal_coordinates);
        let levi_civita_symbol = levi_civita::<1, 1, 1>();
        let normalization = basis_vectors[0].cross(&basis_vectors[1]).norm();
        let normal_vector = basis_vectors[0].cross(&basis_vectors[1])/normalization;
        let identity = TensorRank2::<3, 1, 1>::identity();
        Self::calculate_standard_gradient_operator().iter()
        .map(|standard_gradient_operator_a|
            levi_civita_symbol.iter()
            .map(|levi_civita_symbol_m|
                identity.iter()
                .zip(normal_vector.iter())
                .map(|(identity_i, normal_vector_i)|
                    levi_civita_symbol_m.iter()
                    .zip(basis_vectors[0].iter().zip(basis_vectors[1].iter()))
                    .map(|(levi_civita_symbol_mn, (basis_vector_0_n, basis_vector_1_n))|
                        levi_civita_symbol_mn.iter()
                        .zip(identity_i.iter().zip(normal_vector.iter()))
                        .map(|(levi_civita_symbol_mno, (identity_io, normal_vector_o))|
                            levi_civita_symbol_mno*(identity_io - normal_vector_i*normal_vector_o)
                        ).sum::<Scalar>()*(standard_gradient_operator_a[0]*basis_vector_1_n
                                         - standard_gradient_operator_a[1]*basis_vector_0_n)
                    ).sum::<Scalar>()/normalization
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
                    .zip(basis_vectors[0].iter().zip(basis_vectors[1].iter()))
                    .map(|(levi_civita_symbol_mn, (basis_vector_0_n, basis_vector_1_n))|
                        levi_civita_symbol_mn.iter()
                        .zip(identity_i.iter().zip(normal_vector.iter()))
                        .map(|(levi_civita_symbol_mno, (identity_io, normal_vector_o))|
                            levi_civita_symbol_mno*(identity_io - normal_vector_i*normal_vector_o)
                        ).sum::<Scalar>()*(standard_gradient_operator_a[0]*basis_vector_1_n
                                         - standard_gradient_operator_a[1]*basis_vector_0_n)
                    ).sum::<Scalar>()*nodal_velocity_a_m
                ).sum::<Scalar>()
            ).sum::<Scalar>()/normalization
        ).collect()
    }
    // fn dummy_convert<const I: usize>(nodal_coordinates: &Coordinates<I, N>) -> Coordinates<I, O>
    // {
    //     nodal_coordinates.iter()
    //     .map(|nodal_coordinate|
    //         nodal_coordinate.iter().copied().collect()
    //     ).collect()
    // }
    fn get_reference_normal(&self) -> &ReferenceNormal;
}

pub trait ElasticLinearSurfaceElement<'a, C, const G: usize, const M: usize, const N: usize, const O: usize>
where
    C: Elastic<'a>,
    Self: LinearSurfaceElement<'a, C, G, M, N, O>
{
    fn calculate_nodal_forces_linear_surface_element(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalForces<N>
    {
        let first_piola_kirchoff_stress = self.get_constitutive_models()[0]
        .calculate_first_piola_kirchoff_stress(
            &self.calculate_deformation_gradient(nodal_coordinates)
        );
        self.get_gradient_vectors().iter().map(|gradient_vector|
            &first_piola_kirchoff_stress * gradient_vector
        ).collect()
        // let first_piola_kirchoff_stress = self.get_constitutive_models()[0]
        // .calculate_first_piola_kirchoff_stress(
        //     &self.calculate_deformation_gradient(nodal_coordinates)
        // );
        // let traction = first_piola_kirchoff_stress * self.get_reference_normal();
        // self.get_gradient_vectors().iter()
        // .zip(Self::calculate_normal_gradients(
        //     &Self::dummy_convert(nodal_coordinates)
        // ).iter())
        // .map(|(gradient_vector_a, normal_gradient_a)|
        //     &first_piola_kirchoff_stress * gradient_vector_a +
        //     normal_gradient_a * (&first_piola_kirchoff_stress * self.get_reference_normal())
        // ).collect()
    }
    fn calculate_nodal_stiffnesses_linear_surface_element(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalStiffnesses<N>
    {
        let first_piola_kirchoff_tangent_stiffness = self.get_constitutive_models()[0]
        .calculate_first_piola_kirchoff_tangent_stiffness(
            &self.calculate_deformation_gradient(nodal_coordinates)
        );
        self.get_gradient_vectors().iter()
        .map(|gradient_vector_a|
            self.get_gradient_vectors().iter()
            .map(|gradient_vector_b|
                first_piola_kirchoff_tangent_stiffness
                .contract_second_fourth_indices_with_first_indices_of(
                    gradient_vector_b, gradient_vector_a
                ).transpose()
            ).collect()
        ).collect()
    }
}

macro_rules! linear_surface_element_boilerplate
{
    ($element: ident) =>
    {
        crate::fem::block::element::linear::linear_element_boilerplate!($element);
        impl<'a, C> LinearSurfaceElement<'a, C, G, M, N, O> for $element<C>
        where
            C: Constitutive<'a>
        {
            fn get_reference_normal(&self) -> &ReferenceNormal
            {
                &self.reference_normal
            }
        }
        impl<'a, C> ElasticLinearSurfaceElement<'a, C, G, M, N, O> for $element<C>
        where
            C: Elastic<'a>
        {}
    }
}
pub(crate) use linear_surface_element_boilerplate;