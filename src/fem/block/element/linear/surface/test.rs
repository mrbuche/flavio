pub const THICKNESS: crate::mechanics::Scalar = 1.23;

macro_rules! test_linear_surface_element {
    ($element: ident) => {
        crate::fem::block::element::test::setup_for_surface_elements!($element);
        crate::fem::block::element::linear::test::test_linear_element_inner!($element);
        crate::fem::block::element::linear::surface::test::test_linear_surface_element_inner!(
            $element
        );
    };
}
pub(crate) use test_linear_surface_element;

macro_rules! test_linear_surface_element_inner
{
    ($element: ident) =>
    {
        mod linear_surface_element
        {
            use crate::
            {
                fem::block::element::linear::surface::test::test_linear_surface_element_with_constitutive_model,
                math::Convert,
                test::assert_eq_within_tols
            };
            use super::*;
            mod elastic
            {
                use crate::
                {
                    constitutive::solid::elastic::
                    {
                        AlmansiHamel,
                        test::ALMANSIHAMELPARAMETERS
                    }
                };
                use super::*;
                mod almansi_hamel
                {
                    use super::*;
                    test_linear_surface_element_with_constitutive_model!($element, AlmansiHamel, ALMANSIHAMELPARAMETERS);
                }
            }
            mod hyperelastic
            {
                use crate::
                {
                    constitutive::solid::hyperelastic::
                    {
                        ArrudaBoyce,
                        Fung,
                        Gent,
                        MooneyRivlin,
                        NeoHookean,
                        SaintVenantKirchoff,
                        Yeoh,
                        test::
                        {
                            ARRUDABOYCEPARAMETERS,
                            FUNGPARAMETERS,
                            GENTPARAMETERS,
                            MOONEYRIVLINPARAMETERS,
                            NEOHOOKEANPARAMETERS,
                            SAINTVENANTKIRCHOFFPARAMETERS,
                            YEOHPARAMETERS
                        }
                    }
                };
                use super::*;
                mod arruda_boyce
                {
                    use super::*;
                    test_linear_surface_element_with_constitutive_model!($element, ArrudaBoyce, ARRUDABOYCEPARAMETERS);
                }
                mod fung
                {
                    use super::*;
                    test_linear_surface_element_with_constitutive_model!($element, Fung, FUNGPARAMETERS);
                }
                mod gent
                {
                    use super::*;
                    test_linear_surface_element_with_constitutive_model!($element, Gent, GENTPARAMETERS);
                }
                mod mooney_rivlin
                {
                    use super::*;
                    test_linear_surface_element_with_constitutive_model!($element, MooneyRivlin, MOONEYRIVLINPARAMETERS);
                }
                mod neo_hookean
                {
                    use super::*;
                    test_linear_surface_element_with_constitutive_model!($element, NeoHookean, NEOHOOKEANPARAMETERS);
                }
                mod saint_venant_kirchoff
                {
                    use super::*;
                    test_linear_surface_element_with_constitutive_model!($element, SaintVenantKirchoff, SAINTVENANTKIRCHOFFPARAMETERS);
                }
                mod yeoh
                {
                    use super::*;
                    test_linear_surface_element_with_constitutive_model!($element, Yeoh, YEOHPARAMETERS);
                }
            }
            mod elastic_hyperviscous
            {
                use crate::
                {
                    constitutive::solid::elastic_hyperviscous::
                    {
                        AlmansiHamel,
                        test::ALMANSIHAMELPARAMETERS
                    }
                };
                use super::*;
                mod almansi_hamel
                {
                    use super::*;
                    test_linear_surface_element_with_constitutive_model!($element, AlmansiHamel, ALMANSIHAMELPARAMETERS);
                }
            }
            mod hyperviscoelastic
            {
                use crate::
                {
                    constitutive::solid::hyperviscoelastic::
                    {
                        SaintVenantKirchoff,
                        test::SAINTVENANTKIRCHOFFPARAMETERS
                    }
                };
                use super::*;
                mod saint_venant_kirchoff
                {
                    use super::*;
                    test_linear_surface_element_with_constitutive_model!($element, SaintVenantKirchoff, SAINTVENANTKIRCHOFFPARAMETERS);
                }
            }
        }
    }
}
pub(crate) use test_linear_surface_element_inner;

macro_rules! setup_for_test_finite_element_with_elastic_constitutive_model {
    ($element: ident, $constitutive_model: ident, $constitutive_model_parameters: ident) => {
        fn get_element<'a>() -> $element<$constitutive_model<'a>> {
            $element::new(
                $constitutive_model_parameters,
                get_reference_coordinates(),
                &crate::fem::block::element::linear::surface::test::THICKNESS,
            )
        }
        fn get_element_transformed<'a>() -> $element<$constitutive_model<'a>> {
            $element::new(
                $constitutive_model_parameters,
                get_reference_coordinates_transformed(),
                &crate::fem::block::element::linear::surface::test::THICKNESS,
            )
        }
    };
}
pub(crate) use setup_for_test_finite_element_with_elastic_constitutive_model;

macro_rules! setup_for_test_finite_element_block_with_elastic_constitutive_model {
    ($block: ident, $element: ident, $constitutive_model: ident, $constitutive_model_parameters: ident) => {
        fn get_block<'a>() -> $block<D, E, $element<$constitutive_model<'a>>, G, N> {
            $block::<D, E, $element<$constitutive_model<'a>>, G, N>::new(
                $constitutive_model_parameters,
                get_connectivity(),
                get_reference_coordinates_block(),
                crate::fem::block::element::linear::surface::test::THICKNESS,
            )
        }
        fn get_block_transformed<'a>() -> $block<D, E, $element<$constitutive_model<'a>>, G, N> {
            $block::<D, E, $element<$constitutive_model<'a>>, G, N>::new(
                $constitutive_model_parameters,
                get_connectivity(),
                get_reference_coordinates_transformed_block(),
                crate::fem::block::element::linear::surface::test::THICKNESS,
            )
        }
    };
}
pub(crate) use setup_for_test_finite_element_block_with_elastic_constitutive_model;

macro_rules! setup_for_test_linear_element_with_constitutive_model {
    ($element: ident, $constitutive_model: ident, $constitutive_model_parameters: ident) => {
        fn get_element<'a>() -> $element<$constitutive_model<'a>> {
            $element::new(
                $constitutive_model_parameters,
                get_reference_coordinates(),
                &crate::fem::block::element::linear::surface::test::THICKNESS,
            )
        }
        fn get_element_transformed<'a>() -> $element<$constitutive_model<'a>> {
            $element::new(
                $constitutive_model_parameters,
                get_reference_coordinates_transformed(),
                &crate::fem::block::element::linear::surface::test::THICKNESS,
            )
        }
    };
}
pub(crate) use setup_for_test_linear_element_with_constitutive_model;

macro_rules! setup_for_test_linear_surface_element_with_constitutive_model {
    ($element: ident, $constitutive_model: ident, $constitutive_model_parameters: ident) => {
        fn get_basis(is_deformed: bool, is_transformed: bool) -> Basis<1> {
            if is_transformed {
                if is_deformed {
                    $element::<$constitutive_model>::calculate_basis(
                        &(get_rotation_current_configuration() * get_coordinates()),
                    )
                } else {
                    $element::<$constitutive_model>::calculate_basis(
                        &(get_rotation_reference_configuration() * get_reference_coordinates())
                            .into(),
                    )
                }
            } else {
                if is_deformed {
                    $element::<$constitutive_model>::calculate_basis(&get_coordinates())
                } else {
                    $element::<$constitutive_model>::calculate_basis(
                        &get_reference_coordinates().into(),
                    )
                }
            }
        }
        fn get_dual_basis(is_deformed: bool, is_transformed: bool) -> Basis<1> {
            if is_transformed {
                if is_deformed {
                    $element::<$constitutive_model>::calculate_dual_basis(
                        &(get_rotation_current_configuration() * get_coordinates()),
                    )
                } else {
                    $element::<$constitutive_model>::calculate_dual_basis(
                        &(get_rotation_reference_configuration() * get_reference_coordinates())
                            .into(),
                    )
                }
            } else {
                if is_deformed {
                    $element::<$constitutive_model>::calculate_dual_basis(&get_coordinates())
                } else {
                    $element::<$constitutive_model>::calculate_dual_basis(
                        &get_reference_coordinates().into(),
                    )
                }
            }
        }
        fn get_normal(is_deformed: bool, is_transformed: bool) -> Normal {
            if is_transformed {
                if is_deformed {
                    $element::<$constitutive_model>::calculate_normal(
                        &(get_rotation_current_configuration() * get_coordinates()),
                    )
                } else {
                    $element::<$constitutive_model>::calculate_normal(
                        &(get_rotation_reference_configuration() * get_reference_coordinates())
                            .into(),
                    )
                }
            } else {
                if is_deformed {
                    $element::<$constitutive_model>::calculate_normal(&get_coordinates())
                } else {
                    $element::<$constitutive_model>::calculate_normal(
                        &get_reference_coordinates().into(),
                    )
                }
            }
        }
        fn get_normal_gradients(is_deformed: bool, is_transformed: bool) -> NormalGradients<O> {
            if is_transformed {
                if is_deformed {
                    $element::<$constitutive_model>::calculate_normal_gradients(
                        &(get_rotation_current_configuration() * get_coordinates()),
                    )
                } else {
                    $element::<$constitutive_model>::calculate_normal_gradients(
                        &(get_rotation_reference_configuration() * get_reference_coordinates())
                            .into(),
                    )
                }
            } else {
                if is_deformed {
                    $element::<$constitutive_model>::calculate_normal_gradients(&get_coordinates())
                } else {
                    $element::<$constitutive_model>::calculate_normal_gradients(
                        &get_reference_coordinates().into(),
                    )
                }
            }
        }
        fn get_normal_gradients_from_finite_difference(is_deformed: bool) -> NormalGradients<O> {
            let mut finite_difference = 0.0;
            (0..O)
                .map(|a| {
                    (0..3)
                        .map(|m| {
                            (0..3)
                                .map(|i| {
                                    let mut nodal_coordinates = if is_deformed {
                                        get_coordinates()
                                    } else {
                                        get_reference_coordinates().into()
                                    };
                                    nodal_coordinates[a][m] += 0.5 * EPSILON;
                                    finite_difference =
                                        $element::<$constitutive_model>::calculate_normal(
                                            &nodal_coordinates,
                                        )[i];
                                    nodal_coordinates[a][m] -= EPSILON;
                                    finite_difference -=
                                        $element::<$constitutive_model>::calculate_normal(
                                            &nodal_coordinates,
                                        )[i];
                                    finite_difference / EPSILON
                                })
                                .collect()
                        })
                        .collect()
                })
                .collect()
        }
        fn get_normal_rate(is_deformed: bool, is_transformed: bool) -> NormalRate {
            if is_transformed {
                if is_deformed {
                    $element::<$constitutive_model>::calculate_normal_rate(
                        &(get_rotation_current_configuration() * get_coordinates()),
                        &(get_rotation_current_configuration() * get_velocities()
                            + get_rotation_rate_current_configuration() * get_coordinates()),
                    )
                } else {
                    $element::<$constitutive_model>::calculate_normal_rate(
                        &(get_rotation_reference_configuration() * get_reference_coordinates())
                            .into(),
                        &NodalVelocities::zero(),
                    )
                }
            } else {
                if is_deformed {
                    $element::<$constitutive_model>::calculate_normal_rate(
                        &get_coordinates(),
                        &get_velocities(),
                    )
                } else {
                    $element::<$constitutive_model>::calculate_normal_rate(
                        &get_reference_coordinates().into(),
                        &NodalVelocities::zero(),
                    )
                }
            }
        }
        fn get_normal_rate_from_finite_difference(is_deformed: bool) -> NormalRate {
            let mut finite_difference = 0.0;
            (0..3)
                .map(|i| {
                    get_velocities()
                        .iter()
                        .enumerate()
                        .map(|(a, velocity_a)| {
                            velocity_a
                                .iter()
                                .enumerate()
                                .map(|(k, velocity_a_k)| {
                                    let mut nodal_coordinates = if is_deformed {
                                        get_coordinates()
                                    } else {
                                        get_reference_coordinates().into()
                                    };
                                    nodal_coordinates[a][k] += 0.5 * EPSILON;
                                    finite_difference =
                                        $element::<$constitutive_model>::calculate_normal(
                                            &nodal_coordinates,
                                        )[i];
                                    nodal_coordinates[a][k] -= EPSILON;
                                    finite_difference -=
                                        $element::<$constitutive_model>::calculate_normal(
                                            &nodal_coordinates,
                                        )[i];
                                    finite_difference / EPSILON * velocity_a_k
                                })
                                .sum::<Scalar>()
                        })
                        .sum()
                })
                .collect()
        }
        fn get_normal_tangents(is_deformed: bool, is_transformed: bool) -> NormalTangents<O> {
            if is_transformed {
                if is_deformed {
                    $element::<$constitutive_model>::calculate_normal_tangents(
                        &(get_rotation_current_configuration() * get_coordinates()),
                    )
                } else {
                    $element::<$constitutive_model>::calculate_normal_tangents(
                        &(get_rotation_reference_configuration() * get_reference_coordinates())
                            .into(),
                    )
                }
            } else {
                if is_deformed {
                    $element::<$constitutive_model>::calculate_normal_tangents(&get_coordinates())
                } else {
                    $element::<$constitutive_model>::calculate_normal_tangents(
                        &get_reference_coordinates().into(),
                    )
                }
            }
        }
        fn get_normal_tangents_from_finite_difference(is_deformed: bool) -> NormalTangents<O> {
            let mut finite_difference = 0.0;
            (0..O)
                .map(|a| {
                    (0..O)
                        .map(|b| {
                            (0..3)
                                .map(|m| {
                                    (0..3)
                                        .map(|n| {
                                            (0..3)
                                                .map(|i| {
                                                    let mut nodal_coordinates = if is_deformed {
                                                        get_coordinates()
                                                    } else {
                                                        get_reference_coordinates().into()
                                                    };
                                                    nodal_coordinates[b][n] += 0.5 * EPSILON;
                                                    finite_difference = $element::<
                                                        $constitutive_model,
                                                    >::calculate_normal_gradients(
                                                        &nodal_coordinates,
                                                    )[a][m][i];
                                                    nodal_coordinates[b][n] -= EPSILON;
                                                    finite_difference -= $element::<
                                                        $constitutive_model,
                                                    >::calculate_normal_gradients(
                                                        &nodal_coordinates,
                                                    )[a][m][i];
                                                    finite_difference / EPSILON
                                                })
                                                .collect()
                                        })
                                        .collect()
                                })
                                .collect()
                        })
                        .collect()
                })
                .collect()
        }
        #[test]
        fn size() {
            assert_eq!(
                std::mem::size_of::<$element::<$constitutive_model>>(),
                std::mem::size_of::<$constitutive_model>()
                    + std::mem::size_of::<GradientVectors<N>>()
                    + std::mem::size_of::<ReferenceNormal>()
                    + std::mem::size_of::<Scalar>()
            )
        }
    };
}
pub(crate) use setup_for_test_linear_surface_element_with_constitutive_model;

macro_rules! test_linear_surface_element_with_constitutive_model
{
    ($element: ident, $constitutive_model: ident, $constitutive_model_parameters: ident) =>
    {
        setup_for_test_linear_element_with_constitutive_model!($element, $constitutive_model, $constitutive_model_parameters);
        setup_for_test_linear_surface_element_with_constitutive_model!($element, $constitutive_model, $constitutive_model_parameters);
        mod basis
        {
            use super::*;
            mod deformed
            {
                use super::*;
                #[test]
                fn objectivity()
                {
                    get_basis(true, false).iter()
                    .zip(get_basis(true, true).iter())
                    .for_each(|(basis_m, res_basis_m)|
                        basis_m.iter()
                        .zip((get_rotation_current_configuration().transpose() * res_basis_m
                        ).iter())
                        .for_each(|(basis_m_i, res_basis_m_i)|
                            assert_eq_within_tols(basis_m_i, res_basis_m_i)
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
                    get_basis(false, false).iter()
                    .zip(get_basis(false, true).iter())
                    .for_each(|(basis_m, res_basis_m)|
                        basis_m.iter()
                        .zip((get_rotation_reference_configuration().transpose() * res_basis_m.convert()
                        ).iter())
                        .for_each(|(basis_m_i, res_basis_m_i)|
                            assert_eq_within_tols(basis_m_i, res_basis_m_i)
                        )
                    )
                }
            }
        }
        mod dual_basis
        {
            use super::*;
            mod deformed
            {
                use super::*;
                #[test]
                fn basis()
                {
                    get_basis(true, false).iter()
                    .enumerate()
                    .for_each(|(m, basis_m)|
                        get_dual_basis(true, false).iter()
                        .enumerate()
                        .for_each(|(n, dual_basis_n)|
                            assert_eq_within_tols(
                                &(basis_m * dual_basis_n),
                                &((m == n) as u8 as Scalar)
                            )
                        )
                    )
                }
                #[test]
                fn objectivity()
                {
                    get_dual_basis(true, false).iter()
                    .zip(get_dual_basis(true, true).iter())
                    .for_each(|(dual_basis_m, res_dual_basis_m)|
                        dual_basis_m.iter()
                        .zip((get_rotation_current_configuration().transpose() * res_dual_basis_m
                        ).iter())
                        .for_each(|(dual_basis_m_i, res_dual_basis_m_i)|
                            assert_eq_within_tols(dual_basis_m_i, res_dual_basis_m_i)
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
                    get_basis(false, false).iter()
                    .enumerate()
                    .for_each(|(m, basis_m)|
                        get_dual_basis(false, false).iter()
                        .enumerate()
                        .for_each(|(n, dual_basis_n)|
                            assert_eq_within_tols(
                                &(basis_m * dual_basis_n),
                                &((m == n) as u8 as Scalar)
                            )
                        )
                    )
                }
                #[test]
                fn objectivity()
                {
                    get_dual_basis(false, false).iter()
                    .zip(get_dual_basis(false, true).iter())
                    .for_each(|(dual_basis_m, res_dual_basis_m)|
                        dual_basis_m.iter()
                        .zip((get_rotation_reference_configuration().transpose() * res_dual_basis_m.convert()
                        ).iter())
                        .for_each(|(dual_basis_m_i, res_dual_basis_m_i)|
                            assert_eq_within_tols(dual_basis_m_i, res_dual_basis_m_i)
                        )
                    )
                }
            }
        }
        mod normal
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
                    .for_each(|(normal_gradient_a, fd_normal_gradient_a)|
                        normal_gradient_a.iter()
                        .zip(fd_normal_gradient_a.iter())
                        .for_each(|(normal_gradient_a_m, fd_normal_gradient_a_m)|
                            normal_gradient_a_m.iter()
                            .zip(fd_normal_gradient_a_m.iter())
                            .for_each(|(normal_gradient_a_m_i, fd_normal_gradient_a_m_i)|
                                assert!(
                                    (normal_gradient_a_m_i/fd_normal_gradient_a_m_i - 1.0).abs() < EPSILON
                                )
                            )
                        )
                    )
                }
                #[test]
                fn normal()
                {
                    let basis = get_basis(true, false);
                    let normal = get_normal(true, false);
                    assert_eq_within_tols(
                        &(&basis[0] * &normal), &0.0
                    );
                    assert_eq_within_tols(
                        &(&basis[1] * &normal), &0.0
                    );
                }
                #[test]
                fn normalized()
                {
                    assert_eq_within_tols(
                        &get_normal(true, false).norm(), &1.0
                    )
                }
                #[test]
                fn objectivity()
                {
                    get_normal(true, false).iter()
                    .zip((
                        get_rotation_current_configuration().transpose() *
                        get_normal(true, true)
                    ).iter())
                    .for_each(|(normal_i, res_normal_i)|
                        assert_eq_within_tols(normal_i, res_normal_i)
                    )
                }
            }
            mod undeformed
            {
                use super::*;
                #[test]
                fn finite_difference()
                {
                    get_normal_gradients(true, false).iter()
                    .zip(get_normal_gradients_from_finite_difference(true).iter())
                    .for_each(|(normal_gradient_a, fd_normal_gradient_a)|
                        normal_gradient_a.iter()
                        .zip(fd_normal_gradient_a.iter())
                        .for_each(|(normal_gradient_a_i, fd_normal_gradient_a_i)|
                            normal_gradient_a_i.iter()
                            .zip(fd_normal_gradient_a_i.iter())
                            .for_each(|(normal_gradient_a_i_j, fd_normal_gradient_a_i_j)|
                                assert!(
                                    (normal_gradient_a_i_j/fd_normal_gradient_a_i_j - 1.0).abs() < EPSILON ||
                                    (normal_gradient_a_i_j.abs() < EPSILON && fd_normal_gradient_a_i_j.abs() < EPSILON)
                                )
                            )
                        )
                    )
                }
                #[test]
                fn normal()
                {
                    let basis = get_basis(false, false);
                    let normal = get_normal(false, false);
                    assert_eq_within_tols(
                        &(&basis[0] * &normal), &0.0
                    );
                    assert_eq_within_tols(
                        &(&basis[1] * &normal), &0.0
                    );
                }
                #[test]
                fn normalized()
                {
                    assert_eq_within_tols(
                        &get_normal(false, false).norm(), &1.0
                    )
                }
                #[test]
                fn objectivity()
                {
                    get_normal(false, false).iter()
                    .zip((
                        get_rotation_reference_configuration().transpose() *
                        get_normal(false, true).convert()
                    ).iter())
                    .for_each(|(normal_i, res_normal_i)|
                        assert_eq_within_tols(normal_i, res_normal_i)
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
                                            (normal_tangent_ab_mn_i/fd_normal_tangent_ab_mn_i - 1.0).abs() < EPSILON
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
                }
                #[test]
                fn objectivity()
                {
                    get_normal_gradients(false, false).iter()
                    .zip(get_normal_gradients(false, true).iter())
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
                }
            }
        }
        mod normal_rate
        {
            use super::*;
            mod deformed
            {
                use super::*;
                #[test]
                fn finite_difference()
                {
                    get_normal_rate(true, false).iter()
                    .zip(get_normal_rate_from_finite_difference(true).iter())
                    .for_each(|(normal_rate_i, fd_normal_rate_i)|
                        assert!(
                            (normal_rate_i/fd_normal_rate_i - 1.0).abs() < EPSILON ||
                            normal_rate_i.abs() < EPSILON
                        )
                    )
                }
                #[test]
                fn objectivity()
                {
                    get_normal_rate(true, false).iter()
                    .zip((
                        get_rotation_current_configuration().transpose() *
                        get_normal_rate(true, true) +
                        get_rotation_rate_current_configuration().transpose() *
                        get_normal(true, true)
                    ).iter())
                    .for_each(|(normal_rate_i, res_normal_rate_i)|
                        assert_eq_within_tols(normal_rate_i, res_normal_rate_i)
                    )
                }
            }
            mod undeformed
            {
                use super::*;
                #[test]
                fn finite_difference()
                {
                    get_normal_rate(false, false).iter()
                    .zip(get_normal_rate_from_finite_difference(false).iter())
                    .for_each(|(normal_rate_i, fd_normal_rate_i)|
                        assert!(
                            (normal_rate_i/fd_normal_rate_i - 1.0).abs() < EPSILON ||
                            normal_rate_i.abs() < EPSILON
                        )
                    )
                }
                #[test]
                fn objectivity()
                {
                    get_normal_rate(false, false).iter()
                    .zip((
                        get_rotation_reference_configuration().transpose() *
                        get_normal_rate(false, true).convert()
                    ).iter())
                    .for_each(|(normal_rate_i, res_normal_rate_i)|
                        assert_eq_within_tols(normal_rate_i, res_normal_rate_i)
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
                }
            }
        }
        mod reference_normal
        {
            use super::*;
            #[test]
            fn normal()
            {
                let basis = get_dual_basis(false, false);
                let element = get_element();
                assert_eq_within_tols(
                    &(&basis[0].convert() * element.get_reference_normal()), &0.0
                );
                assert_eq_within_tols(
                    &(&basis[1].convert() * element.get_reference_normal()), &0.0
                );
            }
            #[test]
            fn objectivity()
            {
                get_element().get_reference_normal().iter()
                .zip((
                    get_rotation_reference_configuration().transpose() *
                    get_element_transformed().get_reference_normal()
                ).iter())
                .for_each(|(normal_i, res_normal_i)|
                    assert_eq_within_tols(normal_i, res_normal_i)
                )
            }
        }
    }
}
pub(crate) use test_linear_surface_element_with_constitutive_model;
