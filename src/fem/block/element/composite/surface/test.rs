macro_rules! test_composite_surface_element
{
    ($element: ident) =>
    {
        crate::fem::block::element::test::setup_for_surface_elements!($element);
        crate::fem::block::element::composite::test::test_composite_element_inner!($element);
        crate::fem::block::element::composite::surface::test::test_composite_surface_element_inner!($element);
    }
}
pub(crate) use test_composite_surface_element;

macro_rules! test_composite_surface_element_inner
{
    ($element: ident) =>
    {
        mod composite_surface_element
        {
            use crate::
            {
                fem::block::element::composite::surface::test::test_composite_surface_element_with_constitutive_model,
                test::assert_eq_within_tols
            };
            use super::*;
            mod elastic
            {
                use crate::constitutive::solid::elastic::AlmansiHamel;
                use super::*;
                mod almansi_hamel
                {
                    use super::*;
                    test_composite_surface_element_with_constitutive_model!($element, AlmansiHamel);
                }
            }
            mod hyperelastic
            {
                use crate::constitutive::solid::hyperelastic::
                    {
                    ArrudaBoyce,
                    Fung,
                    Gent,
                    MooneyRivlin,
                    NeoHookean,
                    SaintVenantKirchoff,
                    Yeoh
                };
                use super::*;
                mod arruda_boyce
                {
                    use super::*;
                    test_composite_surface_element_with_constitutive_model!($element, ArrudaBoyce);
                }
                mod fung
                {
                    use super::*;
                    test_composite_surface_element_with_constitutive_model!($element, Fung);
                }
                mod gent
                {
                    use super::*;
                    test_composite_surface_element_with_constitutive_model!($element, Gent);
                }
                mod mooney_rivlin
                {
                    use super::*;
                    test_composite_surface_element_with_constitutive_model!($element, MooneyRivlin);
                }
                mod neo_hookean
                {
                    use super::*;
                    test_composite_surface_element_with_constitutive_model!($element, NeoHookean);
                }
                mod saint_venant_kirchoff
                {
                    use super::*;
                    test_composite_surface_element_with_constitutive_model!($element, SaintVenantKirchoff);
                }
                mod yeoh
                {
                    use super::*;
                    test_composite_surface_element_with_constitutive_model!($element, Yeoh);
                }
            }
            mod elastic_hyperviscous
            {
                use crate::constitutive::solid::elastic_hyperviscous::AlmansiHamel;
                use super::*;
                mod almansi_hamel
                {
                    use super::*;
                    test_composite_surface_element_with_constitutive_model!($element, AlmansiHamel);
                }
            }
            mod hyperviscoelastic
            {
                use crate::constitutive::solid::hyperviscoelastic::SaintVenantKirchoff;
                use super::*;
                mod saint_venant_kirchoff
                {
                    use super::*;
                    test_composite_surface_element_with_constitutive_model!($element, SaintVenantKirchoff);
                }
            }
        }
    }
}
pub(crate) use test_composite_surface_element_inner;

macro_rules! setup_for_test_composite_element_with_constitutive_model
{
    ($element: ident, $constitutive_model: ident) =>
    {
        #[test]
        fn size()
        {
            assert_eq!(
                std::mem::size_of::<$element::<$constitutive_model>>(),
                std::mem::size_of::<[$constitutive_model; G]>()
                + std::mem::size_of::<ProjectedGradientVectors<G, N>>()
                + std::mem::size_of::<Scalars<G>>()
                + std::mem::size_of::<ScaledReferenceNormals<G, P>>()
            )
        }
    }
}
pub(crate) use setup_for_test_composite_element_with_constitutive_model;

macro_rules! setup_for_test_composite_surface_element_with_constitutive_model
{
    ($element: ident, $constitutive_model: ident) =>
    {
        fn get_bases(is_deformed: bool, is_transformed: bool) -> Bases<1, P>
        {
            if is_transformed
            {
                if is_deformed
                {
                    $element::<$constitutive_model>::calculate_bases(
                        &(get_rotation_current_configuration() * get_coordinates())
                    )
                }
                else
                {
                    $element::<$constitutive_model>::calculate_bases(
                        &(get_rotation_reference_configuration() * get_reference_coordinates()).convert()
                    )
                }
            }
            else
            {
                if is_deformed
                {
                    $element::<$constitutive_model>::calculate_bases(
                        &get_coordinates()
                    )
                }
                else
                {
                    $element::<$constitutive_model>::calculate_bases(
                        &get_reference_coordinates().convert()
                    )
                }
            }
        }
        fn get_dual_bases(is_deformed: bool, is_transformed: bool) -> Bases<1, P>
        {
            if is_transformed
            {
                if is_deformed
                {
                    $element::<$constitutive_model>::calculate_dual_bases(
                        &(get_rotation_current_configuration() * get_coordinates())
                    )
                }
                else
                {
                    $element::<$constitutive_model>::calculate_dual_bases(
                        &(get_rotation_reference_configuration() * get_reference_coordinates()).convert()
                    )
                }
            }
            else
            {
                if is_deformed
                {
                    $element::<$constitutive_model>::calculate_dual_bases(
                        &get_coordinates()
                    )
                }
                else
                {
                    $element::<$constitutive_model>::calculate_dual_bases(
                        &get_reference_coordinates().convert()
                    )
                }
            }
        }
        fn get_normals(is_deformed: bool, is_transformed: bool) -> Normals<P>
        {
            if is_transformed
            {
                if is_deformed
                {
                    $element::<$constitutive_model>::calculate_normals(
                        &(get_rotation_current_configuration() * get_coordinates())
                    )
                }
                else
                {
                    $element::<$constitutive_model>::calculate_normals(
                        &(get_rotation_reference_configuration() * get_reference_coordinates()).convert()
                    )
                }
            }
            else
            {
                if is_deformed
                {
                    $element::<$constitutive_model>::calculate_normals(
                        &get_coordinates()
                    )
                }
                else
                {
                    $element::<$constitutive_model>::calculate_normals(
                        &get_reference_coordinates().convert()
                    )
                }
            }
        }
        fn get_normal_gradients(is_deformed: bool, is_transformed: bool) -> NormalGradientss<P, O>
        {
            if is_transformed
            {
                if is_deformed
                {
                    $element::<$constitutive_model>::calculate_normal_gradients(
                        &(get_rotation_current_configuration() * get_coordinates())
                    )
                }
                else
                {
                    $element::<$constitutive_model>::calculate_normal_gradients(
                        &(get_rotation_reference_configuration() * get_reference_coordinates()).convert()
                    )
                }
            }
            else
            {
                if is_deformed
                {
                    $element::<$constitutive_model>::calculate_normal_gradients(
                        &get_coordinates()
                    )
                }
                else
                {
                    $element::<$constitutive_model>::calculate_normal_gradients(
                        &get_reference_coordinates().convert()
                    )
                }
            }
        }
        fn get_normal_gradients_from_finite_difference(is_deformed: bool) -> NormalGradientss<P, O>
        {
            let mut finite_difference = 0.0;
            (0..P).map(|p|
                (0..O).map(|a|
                    (0..3).map(|m|
                        (0..3).map(|i|{
                            let mut nodal_coordinates = 
                            if is_deformed
                            {
                                get_coordinates()
                            }
                            else
                            {
                                get_reference_coordinates().convert()
                            };
                            nodal_coordinates[a][m] += 0.5 * EPSILON;
                            finite_difference = $element::<$constitutive_model>::calculate_normals(
                                &nodal_coordinates
                            )[p][i];
                            nodal_coordinates[a][m] -= EPSILON;
                            finite_difference -= $element::<$constitutive_model>::calculate_normals(
                                &nodal_coordinates
                            )[p][i];
                            finite_difference/EPSILON
                        }).collect()
                    ).collect()
                ).collect()
            ).collect()
        }
        fn get_normal_rates(is_deformed: bool, is_transformed: bool) -> NormalRates<P>
        {
            if is_transformed
            {
                if is_deformed
                {
                    $element::<$constitutive_model>::calculate_normal_rates(
                        &(get_rotation_current_configuration() * get_coordinates()),
                        &(get_rotation_current_configuration() * get_velocities() + get_rotation_rate_current_configuration() * get_coordinates())
                    )
                }
                else
                {
                    $element::<$constitutive_model>::calculate_normal_rates(
                        &(get_rotation_reference_configuration() * get_reference_coordinates()).convert(),
                        &NodalVelocities::zero()
                    )
                }
            }
            else
            {
                if is_deformed
                {
                    $element::<$constitutive_model>::calculate_normal_rates(
                        &get_coordinates(),
                        &get_velocities()
                    )
                }
                else
                {
                    $element::<$constitutive_model>::calculate_normal_rates(
                        &get_reference_coordinates().convert(),
                        &NodalVelocities::zero()
                    )
                }
            }
        }
        fn get_normal_rates_from_finite_difference(is_deformed: bool) -> NormalRates<P>
        {
            let mut finite_difference = 0.0;
            (0..P).map(|p|
                (0..3).map(|i|
                    get_velocities().iter().enumerate()
                    .map(|(a, velocity_a)|
                        velocity_a.iter().enumerate()
                        .map(|(k, velocity_a_k)|{
                            let mut nodal_coordinates = 
                            if is_deformed
                            {
                                get_coordinates()
                            }
                            else
                            {
                                get_reference_coordinates().convert()
                            };
                            nodal_coordinates[a][k] += 0.5 * EPSILON;
                            finite_difference = $element::<$constitutive_model>::calculate_normals(
                                &nodal_coordinates
                            )[p][i];
                            nodal_coordinates[a][k] -= EPSILON;
                            finite_difference -= $element::<$constitutive_model>::calculate_normals(
                                &nodal_coordinates
                            )[p][i];
                            finite_difference/EPSILON * velocity_a_k
                        }).sum::<Scalar>()
                    ).sum()
                ).collect()
            ).collect()
        }
        fn get_normal_tangents(is_deformed: bool, is_transformed: bool) -> NormalTangentss<P, O>
        {
            if is_transformed
            {
                if is_deformed
                {
                    $element::<$constitutive_model>::calculate_normal_tangents(
                        &(get_rotation_current_configuration() * get_coordinates())
                    )
                }
                else
                {
                    $element::<$constitutive_model>::calculate_normal_tangents(
                        &(get_rotation_reference_configuration() * get_reference_coordinates()).convert()
                    )
                }
            }
            else
            {
                if is_deformed
                {
                    $element::<$constitutive_model>::calculate_normal_tangents(
                        &get_coordinates()
                    )
                }
                else
                {
                    $element::<$constitutive_model>::calculate_normal_tangents(
                        &get_reference_coordinates().convert()
                    )
                }
            }
        }
        fn get_normal_tangents_from_finite_difference(is_deformed: bool) -> NormalTangentss<P, O>
        {
            let mut finite_difference = 0.0;
            (0..P).map(|p|
                (0..O).map(|a|
                    (0..O).map(|b|
                        (0..3).map(|m|
                            (0..3).map(|n|
                                (0..3).map(|i|{
                                    let mut nodal_coordinates = 
                                    if is_deformed
                                    {
                                        get_coordinates()
                                    }
                                    else
                                    {
                                        get_reference_coordinates().convert()
                                    };
                                    nodal_coordinates[b][n] += 0.5 * EPSILON;
                                    finite_difference = $element::<$constitutive_model>::calculate_normal_gradients(
                                        &nodal_coordinates
                                    )[p][a][m][i];
                                    nodal_coordinates[b][n] -= EPSILON;
                                    finite_difference -= $element::<$constitutive_model>::calculate_normal_gradients(
                                        &nodal_coordinates
                                    )[p][a][m][i];
                                    finite_difference/EPSILON
                                }).collect()
                            ).collect()
                        ).collect()
                    ).collect()
                ).collect()
            ).collect()
        }
        fn get_reference_normals(is_transformed: bool) -> ReferenceNormals<P>
        {
            if is_transformed
            {
                $element::<$constitutive_model>::calculate_reference_normals(
                    &(get_rotation_reference_configuration() * get_reference_coordinates())
                )
            }
            else
            {
                $element::<$constitutive_model>::calculate_reference_normals(
                    &get_reference_coordinates()
                )
            }
        }
    }
}
pub(crate) use setup_for_test_composite_surface_element_with_constitutive_model;

macro_rules! test_composite_surface_element_with_constitutive_model
{
    ($element: ident, $constitutive_model: ident) =>
    {
        setup_for_test_composite_surface_element_with_constitutive_model!($element, $constitutive_model);
        mod bases
        {
            use super::*;
            mod deformed
            {
                use super::*;
                #[test]
                fn objectivity()
                {
                    get_bases(true, false).iter()
                    .zip(get_bases(true, true).iter())
                    .for_each(|(basis, res_basis)|
                        basis.iter()
                        .zip(res_basis.iter())
                        .for_each(|(basis_m, res_basis_m)|
                            basis_m.iter()
                            .zip((get_rotation_current_configuration().transpose() * res_basis_m
                            ).iter())
                            .for_each(|(basis_m_i, res_basis_m_i)|
                                assert_eq_within_tols(basis_m_i, res_basis_m_i)
                            )
                        )
                    )
                }
            }
            mod undeformed
            {
                use super::*;
                #[test]
                fn objectivity()
                {
                    get_bases(false, false).iter()
                    .zip(get_bases(false, true).iter())
                    .for_each(|(basis, res_basis)|
                        basis.iter()
                        .zip(res_basis.iter())
                        .for_each(|(basis_m, res_basis_m)|
                            basis_m.iter()
                            .zip((get_rotation_reference_configuration().transpose() * res_basis_m.convert()
                            ).iter())
                            .for_each(|(basis_m_i, res_basis_m_i)|
                                assert_eq_within_tols(basis_m_i, res_basis_m_i)
                            )
                        )
                    )
                }
            }
        }
        mod dual_bases
        {
            use super::*;
            mod deformed
            {
                use super::*;
                #[test]
                fn basis()
                {
                    get_bases(true, false).iter()
                    .zip(get_dual_bases(true, false).iter())
                    .for_each(|(basis, dual_basis)|
                        basis.iter()
                        .enumerate()
                        .for_each(|(m, basis_m)|
                            dual_basis.iter()
                            .enumerate()
                            .for_each(|(n, dual_basis_n)|
                                assert_eq_within_tols(
                                    &(basis_m * dual_basis_n), 
                                    &((m == n) as u8 as Scalar)
                                )
                            )
                        )
                    )
                }
                #[test]
                fn objectivity()
                {
                    get_dual_bases(true, false).iter()
                    .zip(get_dual_bases(true, true).iter())
                    .for_each(|(dual_basis, res_dual_basis)|
                        dual_basis.iter()
                        .zip(res_dual_basis.iter())
                        .for_each(|(dual_basis_m, res_dual_basis_m)|
                            dual_basis_m.iter()
                            .zip((get_rotation_current_configuration().transpose() * res_dual_basis_m
                            ).iter())
                            .for_each(|(dual_basis_m_i, res_dual_basis_m_i)|
                                assert_eq_within_tols(dual_basis_m_i, res_dual_basis_m_i)
                            )
                        )
                    )
                }
            }
            mod undeformed
            {
                use super::*;
                #[test]
                fn basis()
                {
                    get_bases(false, false).iter()
                    .zip(get_dual_bases(false, false).iter())
                    .for_each(|(basis, dual_basis)|
                        basis.iter()
                        .enumerate()
                        .for_each(|(m, basis_m)|
                            dual_basis.iter()
                            .enumerate()
                            .for_each(|(n, dual_basis_n)|
                                assert_eq_within_tols(
                                    &(basis_m * dual_basis_n), 
                                    &((m == n) as u8 as Scalar)
                                )
                            )
                        )
                    )
                }
                #[test]
                fn objectivity()
                {
                    get_dual_bases(false, false).iter()
                    .zip(get_dual_bases(false, true).iter())
                    .for_each(|(dual_basis, res_dual_basis)|
                        dual_basis.iter()
                        .zip(res_dual_basis.iter())
                        .for_each(|(dual_basis_m, res_dual_basis_m)|
                            dual_basis_m.iter()
                            .zip((get_rotation_reference_configuration().transpose() * res_dual_basis_m.convert()
                            ).iter())
                            .for_each(|(dual_basis_m_i, res_dual_basis_m_i)|
                                assert_eq_within_tols(dual_basis_m_i, res_dual_basis_m_i)
                            )
                        )
                    )
                }
            }
        }
        mod normals
        {
            use super::*;
            mod deformed
            {
                use super::*;
                #[test]
                fn finite_difference()
                {
                    get_normal_gradients(true, false).iter()
                    .zip(get_normal_gradients_from_finite_difference(true).iter())
                    .for_each(|(normal_gradient, fd_normal_gradient)|
                        normal_gradient.iter()
                        .zip(fd_normal_gradient.iter())
                        .for_each(|(normal_gradient_a, fd_normal_gradient_a)|
                            normal_gradient_a.iter()
                            .zip(fd_normal_gradient_a.iter())
                            .for_each(|(normal_gradient_a_m, fd_normal_gradient_a_m)|
                                normal_gradient_a_m.iter()
                                .zip(fd_normal_gradient_a_m.iter())
                                .for_each(|(normal_gradient_a_m_i, fd_normal_gradient_a_m_i)|
                                    assert!(
                                        (normal_gradient_a_m_i/fd_normal_gradient_a_m_i - 1.0).abs() < EPSILON ||
                                        (normal_gradient_a_m_i.abs() < EPSILON && fd_normal_gradient_a_m_i.abs() < EPSILON)
                                    )
                                )
                            )
                        )
                    )
                }
                #[test]
                fn normal()
                {
                    get_bases(true, false).iter()
                    .zip(get_normals(true, false).iter())
                    .for_each(|(basis, normal)|{
                        assert_eq_within_tols(
                            &(&basis[0] * normal), &0.0
                        );
                        assert_eq_within_tols(
                            &(&basis[1] * normal), &0.0
                        );
                    })
                }
                #[test]
                fn normalized()
                {
                    get_normals(true, false).iter()
                    .for_each(|normal|
                        assert_eq_within_tols(
                            &normal.norm(), &1.0
                        )
                    )
                }
                #[test]
                fn objectivity()
                {
                    get_normals(true, false).iter()
                    .zip(get_normals(true, true).iter())
                    .for_each(|(normal, res_normal)|
                        normal.iter()
                        .zip((get_rotation_current_configuration().transpose() * res_normal).iter())
                        .for_each(|(normal_i, res_normal_i)|
                            assert_eq_within_tols(normal_i, res_normal_i)
                        )
                    )
                }
            }
            mod undeformed
            {
                use super::*;
                #[test]
                fn finite_difference()
                {
                    get_normal_gradients(false, false).iter()
                    .zip(get_normal_gradients_from_finite_difference(false).iter())
                    .for_each(|(normal_gradient, fd_normal_gradient)|
                        normal_gradient.iter()
                        .zip(fd_normal_gradient.iter())
                        .for_each(|(normal_gradient_a, fd_normal_gradient_a)|
                            normal_gradient_a.iter()
                            .zip(fd_normal_gradient_a.iter())
                            .for_each(|(normal_gradient_a_m, fd_normal_gradient_a_m)|
                                normal_gradient_a_m.iter()
                                .zip(fd_normal_gradient_a_m.iter())
                                .for_each(|(normal_gradient_a_m_i, fd_normal_gradient_a_m_i)|
                                    assert!(
                                        (normal_gradient_a_m_i/fd_normal_gradient_a_m_i - 1.0).abs() < EPSILON ||
                                        (normal_gradient_a_m_i.abs() < EPSILON && fd_normal_gradient_a_m_i.abs() < EPSILON)
                                    )
                                )
                            )
                        )
                    )
                }
                #[test]
                fn normal()
                {
                    get_bases(false, false).iter()
                    .zip(get_normals(false, false).iter())
                    .for_each(|(basis, normal)|{
                        assert_eq_within_tols(
                            &(&basis[0] * normal), &0.0
                        );
                        assert_eq_within_tols(
                            &(&basis[1] * normal), &0.0
                        );
                    })
                }
                #[test]
                fn normalized()
                {
                    get_normals(false, false).iter()
                    .for_each(|normal|
                        assert_eq_within_tols(
                            &normal.norm(), &1.0
                        )
                    )
                }
                #[test]
                fn objectivity()
                {
                    get_normals(false, false).iter()
                    .zip(get_normals(false, true).iter())
                    .for_each(|(normal, res_normal)|
                        normal.iter()
                        .zip((get_rotation_reference_configuration().transpose() * res_normal.convert()).iter())
                        .for_each(|(normal_i, res_normal_i)|
                            assert_eq_within_tols(normal_i, res_normal_i)
                        )
                    )
                }
            }
        }
        mod normal_gradients
        {
            use super::*;
            mod deformed
            {
                use super::*;
                #[test]
                fn finite_difference()
                {
                    get_normal_tangents(true, false).iter()
                    .zip(get_normal_tangents_from_finite_difference(true).iter())
                    .for_each(|(normal_tangent, fd_normal_tangent)|
                        normal_tangent.iter()
                        .zip(fd_normal_tangent.iter())
                        .for_each(|(normal_tangent_a, fd_normal_tangent_a)|
                            normal_tangent_a.iter()
                            .zip(fd_normal_tangent_a.iter())
                            .for_each(|(normal_tangent_ab, fd_normal_tangent_ab)|
                                normal_tangent_ab.iter()
                                .zip(fd_normal_tangent_ab.iter())
                                .for_each(|(normal_tangent_ab_m, fd_normal_tangent_ab_m)|
                                    normal_tangent_ab_m.iter()
                                    .zip(fd_normal_tangent_ab_m.iter())
                                    .for_each(|(normal_tangent_ab_m_i, fd_normal_tangent_ab_m_i)|
                                        normal_tangent_ab_m_i.iter()
                                        .zip(fd_normal_tangent_ab_m_i.iter())
                                        .for_each(|(normal_tangent_ab_mn_i, fd_normal_tangent_ab_mn_i)|
                                            assert!(
                                                (normal_tangent_ab_mn_i/fd_normal_tangent_ab_mn_i - 1.0).abs() < EPSILON ||
                                                (normal_tangent_ab_mn_i.abs() < EPSILON && fd_normal_tangent_ab_mn_i.abs() < EPSILON)
                                            )
                                        )
                                    )
                                )
                            )
                        )
                    )
                }
                #[test]
                fn objectivity()
                {
                    get_normal_gradients(true, false).iter()
                    .zip(get_normal_gradients(true, true).iter())
                    .for_each(|(normal_gradient, res_normal_gradient)|
                        normal_gradient.iter()
                        .zip(res_normal_gradient.iter())
                        .for_each(|(normal_gradient_a, res_normal_gradient_a)|
                            normal_gradient_a.iter()
                            .zip((
                                get_rotation_current_configuration().transpose() *
                                res_normal_gradient_a *
                                get_rotation_current_configuration()
                            ).iter())
                            .for_each(|(normal_gradient_a_i, res_normal_gradient_a_i)|
                                normal_gradient_a_i.iter()
                                .zip(res_normal_gradient_a_i.iter())
                                .for_each(|(normal_gradient_a_i_j, res_normal_gradient_a_i_j)|
                                    assert_eq_within_tols(normal_gradient_a_i_j, res_normal_gradient_a_i_j)
                                )
                            )
                        )
                    )
                }
            }
            mod undeformed
            {
                use super::*;
                #[test]
                fn finite_difference()
                {
                    get_normal_tangents(false, false).iter()
                    .zip(get_normal_tangents_from_finite_difference(false).iter())
                    .for_each(|(normal_tangent, fd_normal_tangent)|
                        normal_tangent.iter()
                        .zip(fd_normal_tangent.iter())
                        .for_each(|(normal_tangent_a, fd_normal_tangent_a)|
                            normal_tangent_a.iter()
                            .zip(fd_normal_tangent_a.iter())
                            .for_each(|(normal_tangent_ab, fd_normal_tangent_ab)|
                                normal_tangent_ab.iter()
                                .zip(fd_normal_tangent_ab.iter())
                                .for_each(|(normal_tangent_ab_i, fd_normal_tangent_ab_i)|
                                    normal_tangent_ab_i.iter()
                                    .zip(fd_normal_tangent_ab_i.iter())
                                    .for_each(|(normal_tangent_ab_i_m, fd_normal_tangent_ab_i_m)|
                                        normal_tangent_ab_i_m.iter()
                                        .zip(fd_normal_tangent_ab_i_m.iter())
                                        .for_each(|(normal_tangent_ab_i_mn, fd_normal_tangent_ab_i_mn)|
                                            assert!(
                                                (normal_tangent_ab_i_mn/fd_normal_tangent_ab_i_mn - 1.0).abs() < EPSILON ||
                                                (normal_tangent_ab_i_mn.abs() < EPSILON && fd_normal_tangent_ab_i_mn.abs() < EPSILON)
                                            )
                                        )
                                    )
                                )
                            )
                        )
                    )
                }
                #[test]
                fn objectivity()
                {
                    get_normal_gradients(false, false).iter()
                    .zip(get_normal_gradients(false, true).iter())
                    .for_each(|(normal_gradient, res_normal_gradient)|
                        normal_gradient.iter()
                        .zip(res_normal_gradient.iter())
                        .for_each(|(normal_gradient_a, res_normal_gradient_a)|
                            normal_gradient_a.iter()
                            .zip((
                                get_rotation_reference_configuration().transpose() *
                                res_normal_gradient_a.convert() *
                                get_rotation_reference_configuration()
                            ).iter())
                            .for_each(|(normal_gradient_a_i, res_normal_gradient_a_i)|
                                normal_gradient_a_i.iter()
                                .zip(res_normal_gradient_a_i.iter())
                                .for_each(|(normal_gradient_a_i_j, res_normal_gradient_a_i_j)|
                                    assert_eq_within_tols(normal_gradient_a_i_j, res_normal_gradient_a_i_j)
                                )
                            )
                        )
                    )
                }
            }
        }
        mod normal_rates
        {
            use super::*;
            mod deformed
            {
                use super::*;
                #[test]
                fn finite_difference()
                {
                    get_normal_rates(true, false).iter()
                    .zip(get_normal_rates_from_finite_difference(true).iter())
                    .for_each(|(normal_rate, fd_normal_rate)|
                        normal_rate.iter()
                        .zip(fd_normal_rate.iter())
                        .for_each(|(normal_rate_i, fd_normal_rate_i)|
                            assert!(
                                (normal_rate_i/fd_normal_rate_i - 1.0).abs() < EPSILON ||
                                normal_rate_i.abs() < EPSILON
                            )
                        )
                    )
                }
                #[test]
                fn objectivity()
                {
                    get_normal_rates(true, false).iter()
                    .zip(get_normal_rates(true, true).iter()
                    .zip(get_normals(true, true).iter()))
                    .for_each(|(normal_rate, (res_normal_rate, res_normal))|
                        normal_rate.iter()
                        .zip((
                            get_rotation_current_configuration().transpose() *
                            res_normal_rate + 
                            get_rotation_rate_current_configuration().transpose() *
                            res_normal
                        ).iter())
                        .for_each(|(normal_rate_i, res_normal_rate_i)|
                            assert_eq_within_tols(normal_rate_i, res_normal_rate_i)
                        )
                    )
                }
            }
            mod undeformed
            {
                use super::*;
                #[test]
                fn finite_difference()
                {
                    get_normal_rates(false, false).iter()
                    .zip(get_normal_rates_from_finite_difference(false).iter())
                    .for_each(|(normal_rate, fd_normal_rate)|
                        normal_rate.iter()
                        .zip(fd_normal_rate.iter())
                        .for_each(|(normal_rate_i, fd_normal_rate_i)|
                            assert!(
                                (normal_rate_i/fd_normal_rate_i - 1.0).abs() < EPSILON ||
                                normal_rate_i.abs() < EPSILON
                            )
                        )
                    )
                }
                #[test]
                fn objectivity()
                {
                    get_normal_rates(false, false).iter()
                    .zip(get_normal_rates(false, true).iter())
                    .for_each(|(normal_rate, res_normal_rate)|
                        normal_rate.iter()
                        .zip((
                            get_rotation_reference_configuration().transpose() *
                            res_normal_rate.convert()
                        ).iter())
                        .for_each(|(normal_rate_i, res_normal_rate_i)|
                            assert_eq_within_tols(normal_rate_i, res_normal_rate_i)
                        )
                    )
                }
            }
        }
        mod normal_tangents
        {
            use super::*;
            mod deformed
            {
                use super::*;
                #[test]
                fn objectivity()
                {
                    let rotation_transpose = get_rotation_current_configuration().transpose();
                    get_normal_tangents(true, false).iter()
                    .zip(get_normal_tangents(true, true).iter())
                    .for_each(|(normal_tangents, res_normal_tangents)|
                        normal_tangents.iter()
                        .zip(res_normal_tangents.iter())
                        .for_each(|(normal_tangent_a, res_normal_tangent_a)|
                            normal_tangent_a.iter()
                            .zip(res_normal_tangent_a.iter())
                            .for_each(|(normal_tangent_ab, res_normal_tangent_ab)|
                                rotation_transpose.iter()
                                .map(|rotation_transpose_m|
                                    rotation_transpose.iter()
                                    .map(|rotation_transpose_n|
                                        rotation_transpose.iter()
                                        .map(|rotation_transpose_k|
                                            rotation_transpose_m.iter()
                                            .zip(res_normal_tangent_ab.iter())
                                            .map(|(rotation_transpose_mo, res_normal_tangent_ab_o)|
                                                rotation_transpose_n.iter()
                                                .zip(res_normal_tangent_ab_o.iter())
                                                .map(|(rotation_transpose_np, res_normal_tangent_ab_op)|
                                                    rotation_transpose_k.iter()
                                                    .zip(res_normal_tangent_ab_op.iter())
                                                    .map(|(rotation_transpose_kq, res_normal_tangent_ab_opq)|
                                                        res_normal_tangent_ab_opq * rotation_transpose_mo *
                                                        rotation_transpose_np * rotation_transpose_kq
                                                    ).sum::<Scalar>()
                                                ).sum::<Scalar>()
                                            ).sum()
                                        ).collect()
                                    ).collect()
                                ).collect::<crate::math::TensorRank3<3, 1, 1, 1>>()
                                .iter()
                                .zip(normal_tangent_ab.iter())
                                .for_each(|(rez_normal_tangent_ab_m, normal_tangent_ab_m)|
                                    normal_tangent_ab_m.iter()
                                    .zip(rez_normal_tangent_ab_m.iter())
                                    .for_each(|(normal_tangent_ab_mn, rez_normal_tangent_ab_mn)|
                                        normal_tangent_ab_mn.iter()
                                        .zip(rez_normal_tangent_ab_mn.iter())
                                        .for_each(|(normal_tangent_ab_mn_k, rez_normal_tangent_ab_mn_k)|
                                            assert_eq_within_tols(
                                                normal_tangent_ab_mn_k, rez_normal_tangent_ab_mn_k
                                            )
                                        )
                                    )
                                )
                            )
                        )
                    )
                }
            }
            mod undeformed
            {
                use super::*;
                #[test]
                fn objectivity()
                {
                    let rotation_transpose = get_rotation_reference_configuration().transpose();
                    get_normal_tangents(false, false).iter()
                    .zip(get_normal_tangents(false, true).iter())
                    .for_each(|(normal_tangents, res_normal_tangents)|
                        normal_tangents.iter()
                        .zip(res_normal_tangents.iter())
                        .for_each(|(normal_tangent_a, res_normal_tangent_a)|
                            normal_tangent_a.iter()
                            .zip(res_normal_tangent_a.iter())
                            .for_each(|(normal_tangent_ab, res_normal_tangent_ab)|
                                rotation_transpose.iter()
                                .map(|rotation_transpose_m|
                                    rotation_transpose.iter()
                                    .map(|rotation_transpose_n|
                                        rotation_transpose.iter()
                                        .map(|rotation_transpose_k|
                                            rotation_transpose_m.iter()
                                            .zip(res_normal_tangent_ab.iter())
                                            .map(|(rotation_transpose_mo, res_normal_tangent_ab_o)|
                                                rotation_transpose_n.iter()
                                                .zip(res_normal_tangent_ab_o.iter())
                                                .map(|(rotation_transpose_np, res_normal_tangent_ab_op)|
                                                    rotation_transpose_k.iter()
                                                    .zip(res_normal_tangent_ab_op.iter())
                                                    .map(|(rotation_transpose_kq, res_normal_tangent_ab_opq)|
                                                        res_normal_tangent_ab_opq * rotation_transpose_mo *
                                                        rotation_transpose_np * rotation_transpose_kq
                                                    ).sum::<Scalar>()
                                                ).sum::<Scalar>()
                                            ).sum()
                                        ).collect()
                                    ).collect()
                                ).collect::<crate::math::TensorRank3<3, 1, 1, 1>>()
                                .iter()
                                .zip(normal_tangent_ab.iter())
                                .for_each(|(rez_normal_tangent_ab_m, normal_tangent_ab_m)|
                                    normal_tangent_ab_m.iter()
                                    .zip(rez_normal_tangent_ab_m.iter())
                                    .for_each(|(normal_tangent_ab_mn, rez_normal_tangent_ab_mn)|
                                        normal_tangent_ab_mn.iter()
                                        .zip(rez_normal_tangent_ab_mn.iter())
                                        .for_each(|(normal_tangent_ab_mn_k, rez_normal_tangent_ab_mn_k)|
                                            assert_eq_within_tols(
                                                normal_tangent_ab_mn_k, rez_normal_tangent_ab_mn_k
                                            )
                                        )
                                    )
                                )
                            )
                        )
                    )
                }
            }
        }
        mod reference_normals
        {
            use super::*;
            #[test]
            fn normals()
            {
                get_dual_bases(false, false).iter()
                .zip(get_reference_normals(false).iter())
                .for_each(|(dual_basis, reference_normal)|{
                    assert_eq_within_tols(
                        &(&dual_basis[0].convert() * reference_normal), &0.0
                    );
                    assert_eq_within_tols(
                        &(&dual_basis[1].convert() * reference_normal), &0.0
                    );
                })
            }
            #[test]
            fn objectivity()
            {
                get_reference_normals(false).iter()
                .zip(get_reference_normals(true).iter())
                .for_each(|(normal, res_normal)|
                    normal.iter()
                    .zip((get_rotation_reference_configuration().transpose() * res_normal).iter())
                    .for_each(|(normal_i, res_normal_i)|
                        assert_eq_within_tols(normal_i, res_normal_i)
                    )
                )
            }
        }
    }
}
pub(crate) use test_composite_surface_element_with_constitutive_model;