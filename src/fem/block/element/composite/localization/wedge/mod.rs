#[cfg(test)]
mod test;

use super::*;

const G: usize = 3;
const M: usize = 2;
const N: usize = 12;
const O: usize = 6;
const P: usize = 4;
const Q: usize = 3;

const INTEGRATION_WEIGHT: Scalar = 1.0/6.0;

pub struct Wedge<C>
{
    constitutive_models: [C; G],
    projected_gradient_vectors: ProjectedGradientVectors<G, N>,
    scaled_composite_jacobians: Scalars<G>,
    scaled_reference_normals: ScaledReferenceNormals<G, P>
}

impl<'a, C> FiniteElement<'a, C, G, N> for Wedge<C>
where
    C: Constitutive<'a>
{
    fn new(constitutive_model_parameters: Parameters<'a>, reference_nodal_coordinates: ReferenceNodalCoordinates<N>) -> Self
    {
        let nodal_coordinates_midplane = Self::calculate_midplane(&reference_nodal_coordinates);
        Self
        {
            constitutive_models: std::array::from_fn(|_| <C>::new(constitutive_model_parameters)),
            projected_gradient_vectors: Self::calculate_projected_gradient_vectors(&reference_nodal_coordinates),
            scaled_composite_jacobians: Self::calculate_scaled_composite_jacobian_at_integration_points(&nodal_coordinates_midplane),
            scaled_reference_normals: Self::calculate_scaled_reference_normals(&nodal_coordinates_midplane)
        }
    }
}

impl<'a, C> CompositeElement<'a, C, G, M, N, O, P, Q> for Wedge<C>
where
    C: Constitutive<'a>
{
    fn calculate_deformation_gradients(&self, nodal_coordinates: &NodalCoordinates<N>) -> DeformationGradients<G>
    {
        self.calculate_deformation_gradients_composite_localization_element(nodal_coordinates)
    }
    fn calculate_deformation_gradient_rates(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> DeformationGradientRates<G>
    {
        self.calculate_deformation_gradient_rates_composite_localization_element(nodal_coordinates, nodal_velocities)
    }
    fn calculate_projected_gradient_vectors(reference_nodal_coordinates: &ReferenceNodalCoordinates<N>) -> ProjectedGradientVectors<G, N>
    {
        let reference_nodal_coordinates_midplane = Self::calculate_midplane(reference_nodal_coordinates);
        let reference_dual_bases = Self::calculate_dual_bases(&reference_nodal_coordinates_midplane);
        let reference_jacobians = Self::calculate_reference_jacobians(&reference_nodal_coordinates_midplane);
        let reference_normals = Self::calculate_reference_normals(&reference_nodal_coordinates_midplane);
        let inverse_projection_matrix = Self::calculate_inverse_projection_matrix(&reference_jacobians);
        let projected_gradient_vectors_midplane =
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
        ).collect::<ProjectedGradientVectors<G, N>>();
        let other_scaled_reference_normals =
        Self::calculate_shape_functions_at_integration_points().iter()
        .map(|shape_function|
            Self::calculate_mixed_shape_function_integrals_products().iter()
            .map(|mixed_shape_function_integrals_products|
                reference_normals.iter()
                .zip(reference_jacobians.iter()
                .zip(mixed_shape_function_integrals_products.iter()))
                .map(|(reference_normal, (reference_jacobian, mixed_shape_function_integrals_product))|
                    reference_normal * ((shape_function * (&inverse_projection_matrix * mixed_shape_function_integrals_product)) * reference_jacobian)
                ).sum()
            ).collect()
        ).collect::<ScaledReferenceNormals<G, O>>();
        let mut projected_gradient_vectors = ProjectedGradientVectors::zero();
        projected_gradient_vectors.iter_mut()
        .zip(projected_gradient_vectors_midplane.iter()
        .zip(other_scaled_reference_normals.iter()))
        .for_each(|(projected_gradient_vectors_g, (projected_gradient_vectors_midplane_g, other_scaled_reference_normals_g))|{
            projected_gradient_vectors_g.iter_mut().skip(3).take(3)
            .zip(projected_gradient_vectors_midplane_g.iter().take(3)
            .zip(other_scaled_reference_normals_g.iter().take(3)))
            .for_each(|(projected_gradient_vector_g_a, (projected_gradient_vector_midplane_g_a, other_scaled_reference_normal_g_a))|
                *projected_gradient_vector_g_a = projected_gradient_vector_midplane_g_a * 0.5 + other_scaled_reference_normal_g_a
            );
            projected_gradient_vectors_g.iter_mut().skip(9)
            .zip(projected_gradient_vectors_midplane_g.iter().skip(3)
            .zip(other_scaled_reference_normals_g.iter().skip(3)))
            .for_each(|(projected_gradient_vector_g_a, (projected_gradient_vector_midplane_g_a, other_scaled_reference_normal_g_a))|
                *projected_gradient_vector_g_a = projected_gradient_vector_midplane_g_a * 0.5 + other_scaled_reference_normal_g_a
            );
            projected_gradient_vectors_g.iter_mut().take(3)
            .zip(projected_gradient_vectors_midplane_g.iter().take(3)
            .zip(other_scaled_reference_normals_g.iter().take(3)))
            .for_each(|(projected_gradient_vector_g_a, (projected_gradient_vector_midplane_g_a, other_scaled_reference_normal_g_a))|
                *projected_gradient_vector_g_a = projected_gradient_vector_midplane_g_a * 0.5 - other_scaled_reference_normal_g_a
            );
            projected_gradient_vectors_g.iter_mut().skip(6).take(3)
            .zip(projected_gradient_vectors_midplane_g.iter().skip(3)
            .zip(other_scaled_reference_normals_g.iter().skip(3)))
            .for_each(|(projected_gradient_vector_g_a, (projected_gradient_vector_midplane_g_a, other_scaled_reference_normal_g_a))|
                *projected_gradient_vector_g_a = projected_gradient_vector_midplane_g_a * 0.5 - other_scaled_reference_normal_g_a
            );
        });
        projected_gradient_vectors
    }
    composite_surface_element_boilerplate_inner!{}
}

impl<'a, C> CompositeLocalizationElement<'a, C, G, M, N, O, P, Q> for Wedge<C>
where
    C: Constitutive<'a>
{
    fn calculate_midplane<const I: usize>(nodal_coordinates: &Coordinates<I, N>) -> Coordinates<I, O>
    {
        nodal_coordinates.iter().skip(3).take(3)
        .chain(nodal_coordinates.iter().skip(9))
        .zip(nodal_coordinates.iter().take(3)
        .chain(nodal_coordinates.iter().skip(6).take(3)))
        .map(|(nodal_coordinates_top, nodal_coordinates_bottom)|
            nodal_coordinates_top.iter()
            .zip(nodal_coordinates_bottom.iter())
            .map(|(nodal_coordinates_top_i, nodal_coordinates_bottom_i)|
                (nodal_coordinates_top_i + nodal_coordinates_bottom_i) * 0.5
            ).collect()
        ).collect()
    }
    fn calculate_mixed_shape_function_integrals_products() -> ShapeFunctionIntegralsProductsMixed<O, P>
    {
        ShapeFunctionIntegralsProductsMixed::new([[
            [12.0,  2.0,  2.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0]
            ], [
            [ 0.0,  0.0,  0.0],
            [ 2.0, 12.0,  2.0],
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0]
            ], [
            [ 0.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 2.0,  2.0, 12.0],
            [ 0.0,  0.0,  0.0]
            ], [
            [ 0.0,  0.0,  0.0],
            [ 4.0, 10.0,  2.0],
            [ 0.0,  0.0,  0.0],
            [ 6.0,  6.0,  4.0]
            ], [
            [10.0,  4.0,  2.0],
            [ 2.0, 10.0,  4.0],
            [ 2.0,  4.0, 10.0],
            [ 4.0,  6.0,  6.0]
            ], [
            [10.0,  2.0,  4.0],
            [ 0.0,  0.0,  0.0],
            [ 4.0,  2.0, 10.0],
            [ 6.0,  4.0,  6.0]
        ]])
    }
}
impl<'a, C> ElasticFiniteElement<'a, C, G, N> for Wedge<C>
where
    C: Elastic<'a>
{
    fn calculate_nodal_forces(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalForces<N>
    {
        let identity = TensorRank2::<3, 1, 1>::identity();
        let normal_gradients = Self::calculate_normal_gradients(
            &Self::calculate_midplane(nodal_coordinates)
        );
        let objects = self.calculate_objects(&normal_gradients);
        self.get_constitutive_models().iter()
        .zip(self.calculate_deformation_gradients(nodal_coordinates).iter())
        .map(|(constitutive_model, deformation_gradient)|
            constitutive_model.calculate_first_piola_kirchoff_stress(deformation_gradient)
        ).collect::<FirstPiolaKirchoffStresses<G>>().iter()
        .zip(self.get_projected_gradient_vectors().iter()
        .zip(self.get_scaled_composite_jacobians().iter()
        .zip(objects.iter())))
        .map(|(first_piola_kirchoff_stress, (projected_gradient_vectors, (scaled_composite_jacobian, objects)))|
            projected_gradient_vectors.iter()
            .zip(objects.iter().take(3)
            .chain(objects.iter().take(3))
            .chain(objects.iter().skip(3).take(3))
            .chain(objects.iter().skip(3).take(3))
            )
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
                                identity_mi * projected_gradient_vector_j + object_mij * 0.5
                            ) * scaled_composite_jacobian
                        ).sum::<Scalar>()
                    ).sum::<Scalar>()
                ).collect()
            ).collect()
        ).sum()
        // todo!("Factors of 1/2 and stuff that are not present in surface implementation.
        //        Remember why? It's from the normal gradients on the midplane.
        //        So maybe multiply the object by 0.5 in the loop.
        //        And don't forget to chain the normal gradients, since you calculate them on the midplane!
        //        And make sure those chains are consistent with projected gradients vectors chains.")
        //
        // need to chain objects<O> consistent with node list and projected gradient vectors
    }
    fn calculate_nodal_stiffnesses(&self, nodal_coordinates: &NodalCoordinates<N>) -> NodalStiffnesses<N>
    {
        todo!("Factors of 1/2 and stuff that are not present in surface implementation.")
    }
}
impl<'a, C> ViscoelasticFiniteElement<'a, C, G, N> for Wedge<C>
where
    C: Viscoelastic<'a>
{
    fn calculate_nodal_forces(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> NodalForces<N>
    {
        todo!("Factors of 1/2 and stuff that are not present in surface implementation.")
    }
    fn calculate_nodal_stiffnesses(&self, nodal_coordinates: &NodalCoordinates<N>, nodal_velocities: &NodalVelocities<N>) -> NodalStiffnesses<N>
    {
        todo!("Factors of 1/2 and stuff that are not present in surface implementation.")
    }
}

composite_surface_or_localization_element_boilerplate!(Wedge);