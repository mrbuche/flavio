macro_rules! test_finite_element
{
    ($element: ident) =>
    {
        mod element
        {
            use crate::
            {
                EPSILON,
                fem::block::element::test::
                {
                    test_finite_element_with_elastic_constitutive_model,
                    test_finite_element_with_hyperelastic_constitutive_model,
                    test_finite_element_with_elastic_hyperviscous_constitutive_model,
                    test_finite_element_with_hyperviscoelastic_constitutive_model
                },
                math::
                {
                    Convert,
                    TensorRank2
                },
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
                    test_finite_element_with_elastic_constitutive_model!($element, AlmansiHamel, ALMANSIHAMELPARAMETERS);
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
                    test_finite_element_with_hyperelastic_constitutive_model!($element, ArrudaBoyce, ARRUDABOYCEPARAMETERS);
                }
                mod fung
                {
                    use super::*;
                    test_finite_element_with_hyperelastic_constitutive_model!($element, Fung, FUNGPARAMETERS);
                }
                mod gent
                {
                    use super::*;
                    test_finite_element_with_hyperelastic_constitutive_model!($element, Gent, GENTPARAMETERS);
                }
                mod mooney_rivlin
                {
                    use super::*;
                    test_finite_element_with_hyperelastic_constitutive_model!($element, MooneyRivlin, MOONEYRIVLINPARAMETERS);
                }
                mod neo_hookean
                {
                    use super::*;
                    test_finite_element_with_hyperelastic_constitutive_model!($element, NeoHookean, NEOHOOKEANPARAMETERS);
                }
                mod saint_venant_kirchoff
                {
                    use super::*;
                    test_finite_element_with_hyperelastic_constitutive_model!($element, SaintVenantKirchoff, SAINTVENANTKIRCHOFFPARAMETERS);
                }
                mod yeoh
                {
                    use super::*;
                    test_finite_element_with_hyperelastic_constitutive_model!($element, Yeoh, YEOHPARAMETERS);
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
                    test_finite_element_with_elastic_hyperviscous_constitutive_model!($element, AlmansiHamel, ALMANSIHAMELPARAMETERS);
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
                    test_finite_element_with_hyperviscoelastic_constitutive_model!($element, SaintVenantKirchoff, SAINTVENANTKIRCHOFFPARAMETERS);
                }
            }
        }
    }
}
pub(crate) use test_finite_element;

macro_rules! setup_for_element_tests_any_element
{
    ($element: ident) =>
    {
        use crate::mechanics::test::
        {
            get_rotation_current_configuration,
            get_rotation_reference_configuration,
            get_rotation_rate_current_configuration,
            get_translation_current_configuration,
            get_translation_reference_configuration,
            get_translation_rate_current_configuration
        };
        fn get_coordinates_transformed() -> NodalCoordinates<N>
        {
            get_coordinates().iter()
            .map(|coordinate|
                get_rotation_current_configuration() * coordinate
                + get_translation_current_configuration()
            ).collect()
        }
        fn get_reference_coordinates_transformed() -> ReferenceNodalCoordinates<N>
        {
            get_reference_coordinates().iter()
            .map(|reference_coordinate|
                get_rotation_reference_configuration() * reference_coordinate
                + get_translation_reference_configuration()
            ).collect()
        }
        fn get_velocities_transformed() -> NodalVelocities<N>
        {
            get_coordinates().iter()
            .zip(get_velocities().iter())
            .map(|(coordinate, velocity)|
                get_rotation_current_configuration() * velocity
                + get_rotation_rate_current_configuration() * coordinate
                + get_translation_rate_current_configuration()
            ).collect()
        }
    }
}
pub(crate) use setup_for_element_tests_any_element;

macro_rules! setup_for_elements
{
    ($element: ident) =>
    {
        use crate::
        {
            constitutive::solid::elastic::AlmansiHamel,
            mechanics::test::
            {
                get_deformation_gradient,
                get_deformation_gradient_rate
            }
        };
        fn get_coordinates() -> NodalCoordinates<N>
        {
            get_deformation_gradient() * get_reference_coordinates()
        }
        fn get_velocities() -> NodalVelocities<N>
        {
            get_deformation_gradient_rate() * get_reference_coordinates()
        }
        #[test]
        fn size()
        {
            assert_eq!(
                std::mem::size_of::<Tetrahedron::<AlmansiHamel>>(),
                std::mem::size_of::<AlmansiHamel>()
                + std::mem::size_of::<GradientVectors<N>>()
            )
        }
        crate::fem::block::element::test::setup_for_element_tests_any_element!($element);
    }
}
pub(crate) use setup_for_elements;

macro_rules! setup_for_surface_or_localization_elements
{
    ($element: ident) =>
    {
        use crate::
        {
            constitutive::solid::elastic::AlmansiHamel,
            mechanics::
            {
                RotationCurrentConfiguration,
                RotationRateCurrentConfiguration
            },
            EPSILON
        };
        fn get_deformation_gradient() -> DeformationGradient
        {
            get_deformation_gradient_rotation() * get_deformation_gradient_special()
        }
        fn get_deformation_gradient_rate() -> DeformationGradientRate
        {
            get_deformation_gradient_rotation_rate() * get_deformation_gradient_special()
        }
        fn get_deformation_gradient_rotation() -> RotationCurrentConfiguration
        {
            crate::mechanics::test::get_rotation_reference_configuration().convert().transpose()
        }
        fn get_deformation_gradient_rotation_rate() -> RotationRateCurrentConfiguration
        {
            crate::mechanics::FrameSpin::new([
                [ 0.0, -0.3,  0.1],
                [ 0.3,  0.0,  0.5],
                [-0.1, -0.5,  0.0]
            ]) * get_deformation_gradient_rotation()
        }
        #[test]
        fn size()
        {
            assert_eq!(
                std::mem::size_of::<$element::<AlmansiHamel>>(),
                std::mem::size_of::<AlmansiHamel>()
                + std::mem::size_of::<GradientVectors<N>>()
                + std::mem::size_of::<ReferenceNormal>()
            )
        }
        crate::fem::block::element::test::setup_for_element_tests_any_element!($element);
    }
}
pub(crate) use setup_for_surface_or_localization_elements;

macro_rules! setup_for_surface_elements
{
    ($element: ident) =>
    {
        fn get_coordinates() -> NodalCoordinates<N>
        {
            get_deformation_gradient() * get_reference_coordinates()
        }
        fn get_deformation_gradient_special() -> DeformationGradient
        {
            DeformationGradient::new([
                [0.62, 0.20, 0.00],
                [0.32, 0.98, 0.00],
                [0.00, 0.00, 1.00]
            ])
        }
        fn get_velocities() -> NodalVelocities<N>
        {
            get_deformation_gradient_rate() * get_reference_coordinates()
        }
        crate::fem::block::element::test::setup_for_surface_or_localization_elements!($element);
    }
}
pub(crate) use setup_for_surface_elements;

macro_rules! setup_for_localization_elements
{
    ($element: ident) =>
    {
        fn get_coordinates() -> NodalCoordinates<N>
        {
            get_deformation_gradient_rotation() * get_coordinates_unrotated()
        }
        fn get_coordinates_unrotated() -> NodalCoordinates<N>
        {
            let jump = get_jump();
            let mut coordinates = get_deformation_gradient_surface() * get_reference_coordinates();
            coordinates.iter_mut().skip(O)
            .for_each(|coordinate_top_a|
                *coordinate_top_a += &jump
            );
            coordinates
        }
        fn get_deformation_gradient_special() -> DeformationGradient
        {
            let jump = get_jump();
            let mut deformation_gradient = get_deformation_gradient_surface();
            deformation_gradient[0][2] = jump[0];
            deformation_gradient[1][2] = jump[1];
            deformation_gradient[2][2] = jump[2] + 1.0;
            deformation_gradient
        }
        fn get_deformation_gradient_surface() -> DeformationGradient
        {
            DeformationGradient::new([
                [0.62, 0.20, 0.00],
                [0.32, 0.98, 0.00],
                [0.00, 0.00, 1.00]
            ])
        }
        fn get_jump() -> Vector<1>
        {
            Vector::new([1.11, 1.22, 1.33])
        }
        fn get_velocities() -> NodalVelocities<N>
        {
            get_deformation_gradient_rotation_rate() * get_coordinates_unrotated()
        }
        crate::fem::block::element::test::setup_for_surface_or_localization_elements!($element);
    }
}
pub(crate) use setup_for_localization_elements;

macro_rules! test_nodal_forces_and_nodal_stiffnesses
{
    ($element: ident, $constitutive_model: ident, $constitutive_model_parameters: ident) =>
    {
        fn get_element<'a>() -> $element<$constitutive_model<'a>>
        {
            $element::new(
                $constitutive_model_parameters,
                get_reference_coordinates()
            )
        }
        fn get_element_transformed<'a>() -> $element<$constitutive_model<'a>>
        {
            $element::<$constitutive_model>::new
            (
                $constitutive_model_parameters,
                get_reference_coordinates_transformed()
            )
        }
        #[test]
        fn integration_weights_sum_to_one()
        {
            assert_eq!(get_element().get_integration_weights().iter().sum::<Scalar>(), 1.0)
        }
        mod nodal_forces
        {
            use super::*;
            mod deformed
            {
                use super::*;
                #[test]
                fn finite_difference()
                {
                    get_nodal_stiffnesses(true, false).iter()
                    .zip(get_finite_difference_of_nodal_forces(true).iter())
                    .for_each(|(nodal_stiffness_a, fd_nodal_stiffness_a)|
                        nodal_stiffness_a.iter()
                        .zip(fd_nodal_stiffness_a.iter())
                        .for_each(|(nodal_stiffness_ab, fd_nodal_stiffness_ab)|
                            nodal_stiffness_ab.iter()
                            .zip(fd_nodal_stiffness_ab.iter())
                            .for_each(|(nodal_stiffness_ab_i, fd_nodal_stiffness_ab_i)|
                                nodal_stiffness_ab_i.iter()
                                .zip(fd_nodal_stiffness_ab_i.iter())
                                .for_each(|(nodal_stiffness_ab_ij, fd_nodal_stiffness_ab_ij)|
                                    assert!(
                                        (nodal_stiffness_ab_ij/fd_nodal_stiffness_ab_ij - 1.0).abs() < EPSILON
                                    )
                                )
                            )
                        )
                    )
                }
                #[test]
                fn objectivity()
                {
                    get_nodal_forces(true, false, true).iter()
                    .zip(get_nodal_forces(true, true, true).iter())
                    .for_each(|(nodal_force, res_nodal_force)|
                        nodal_force.iter()
                        .zip(res_nodal_force.iter())
                        .for_each(|(nodal_force_i, res_nodal_force_i)|
                            assert_eq_within_tols(
                                nodal_force_i, res_nodal_force_i
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
                    get_nodal_stiffnesses(false, false).iter()
                    .zip(get_finite_difference_of_nodal_forces(false).iter())
                    .for_each(|(nodal_stiffness_a, fd_nodal_stiffness_a)|
                        nodal_stiffness_a.iter()
                        .zip(fd_nodal_stiffness_a.iter())
                        .for_each(|(nodal_stiffness_ab, fd_nodal_stiffness_ab)|
                            nodal_stiffness_ab.iter()
                            .zip(fd_nodal_stiffness_ab.iter())
                            .for_each(|(nodal_stiffness_ab_i, fd_nodal_stiffness_ab_i)|
                                nodal_stiffness_ab_i.iter()
                                .zip(fd_nodal_stiffness_ab_i.iter())
                                .for_each(|(nodal_stiffness_ab_ij, fd_nodal_stiffness_ab_ij)|
                                    assert!(
                                        (nodal_stiffness_ab_ij/fd_nodal_stiffness_ab_ij - 1.0).abs() < EPSILON ||
                                        (nodal_stiffness_ab_ij.abs() < EPSILON && fd_nodal_stiffness_ab_ij.abs() < EPSILON)
                                    )
                                )
                            )
                        )
                    )
                }
                #[test]
                fn objectivity()
                {
                    get_nodal_forces(false, false, false).iter()
                    .for_each(|nodal_force|
                        nodal_force.iter()
                        .for_each(|nodal_force_i|
                            assert_eq_within_tols(
                                nodal_force_i, &0.0
                            )
                        )
                    )
                }
                #[test]
                fn zero()
                {
                    get_nodal_forces(false, false, false).iter()
                    .for_each(|nodal_force|
                        nodal_force.iter()
                        .for_each(|nodal_force_i|
                            assert_eq_within_tols(
                                nodal_force_i, &0.0
                            )
                        )
                    )
                }
            }
        }
        mod nodal_stiffnesses
        {
            use super::*;
            mod deformed
            {
                use super::*;
                #[test]
                fn objectivity()
                {
                    get_nodal_stiffnesses(true, false).iter()
                    .zip(get_nodal_stiffnesses(true, true).iter())
                    .for_each(|(nodal_stiffness_a, res_nodal_stiffness_a)|
                        nodal_stiffness_a.iter()
                        .zip(res_nodal_stiffness_a.iter())
                        .for_each(|(nodal_stiffness_ab, res_nodal_stiffness_ab)|
                            nodal_stiffness_ab.iter()
                            .zip(res_nodal_stiffness_ab.iter())
                            .for_each(|(nodal_stiffness_ab_i, res_nodal_stiffness_ab_i)|
                                nodal_stiffness_ab_i.iter()
                                .zip(res_nodal_stiffness_ab_i.iter())
                                .for_each(|(nodal_stiffness_ab_ij, res_nodal_stiffness_ab_ij)|
                                    assert_eq_within_tols(
                                        nodal_stiffness_ab_ij, res_nodal_stiffness_ab_ij
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
                    get_nodal_stiffnesses(false, false).iter()
                    .zip(get_nodal_stiffnesses(false, true).iter())
                    .for_each(|(nodal_stiffness_a, res_nodal_stiffness_a)|
                        nodal_stiffness_a.iter()
                        .zip(res_nodal_stiffness_a.iter())
                        .for_each(|(nodal_stiffness_ab, res_nodal_stiffness_ab)|
                            nodal_stiffness_ab.iter()
                            .zip(res_nodal_stiffness_ab.iter())
                            .for_each(|(nodal_stiffness_ab_i, res_nodal_stiffness_ab_i)|
                                nodal_stiffness_ab_i.iter()
                                .zip(res_nodal_stiffness_ab_i.iter())
                                .for_each(|(nodal_stiffness_ab_ij, res_nodal_stiffness_ab_ij)|
                                    assert_eq_within_tols(
                                        nodal_stiffness_ab_ij, res_nodal_stiffness_ab_ij
                                    )
                                )
                            )
                        )
                    )
                }
                #[test]
                fn symmetry()
                {
                    let nodal_stiffness = get_nodal_stiffnesses(false, false);
                    nodal_stiffness.iter()
                    .enumerate()
                    .for_each(|(a, nodal_stiffness_a)|
                        nodal_stiffness_a.iter()
                        .enumerate()
                        .for_each(|(b, nodal_stiffness_ab)|
                            nodal_stiffness_ab.iter()
                            .enumerate()
                            .zip(nodal_stiffness_ab.transpose().iter())
                            .for_each(|((i, nodal_stiffness_ab_i), nodal_stiffness_ab_j)|
                                nodal_stiffness_ab_i.iter()
                                .enumerate()
                                .zip(nodal_stiffness_ab_j.iter())
                                .for_each(|((j, nodal_stiffness_ab_ij), nodal_stiffness_ab_ji)|
                                    if a == b
                                    {
                                        assert_eq_within_tols(
                                            nodal_stiffness_ab_ij,
                                            &nodal_stiffness_ab_ji
                                        )
                                    }
                                    else if i == j
                                    {
                                        assert_eq_within_tols(
                                            nodal_stiffness_ab_ij, &nodal_stiffness[b][a][i][j]
                                        )
                                    }
                                )
                            )
                        )
                    )
                }
            }
        }
    }
}
pub(crate) use test_nodal_forces_and_nodal_stiffnesses;

macro_rules! test_helmholtz_free_energy
{
    ($element: ident, $constitutive_model: ident, $constitutive_model_parameters: ident) =>
    {
        fn get_helmholtz_free_energy(is_deformed: bool, is_rotated: bool) -> Scalar
        {
            if is_rotated
            {
                if is_deformed
                {
                    get_element_transformed()
                    .calculate_helmholtz_free_energy(
                        &get_coordinates_transformed()
                    )
                }
                else
                {
                    get_element_transformed()
                    .calculate_helmholtz_free_energy(
                        &get_reference_coordinates_transformed().convert()
                    )
                }
            }
            else
            {
                if is_deformed
                {
                    get_element()
                    .calculate_helmholtz_free_energy(
                        &get_coordinates()
                    )
                }
                else
                {
                    get_element()
                    .calculate_helmholtz_free_energy(
                        &get_reference_coordinates().convert()
                    )
                }
            }
        }
        fn get_finite_difference_of_helmholtz_free_energy(is_deformed: bool) -> NodalForces<N>
        {
            let element = get_element();
            let mut finite_difference = 0.0;
            (0..N).map(|node|
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
                    nodal_coordinates[node][i] += 0.5 * EPSILON;
                    finite_difference = element.calculate_helmholtz_free_energy(
                        &nodal_coordinates
                    );
                    nodal_coordinates[node][i] -= EPSILON;
                    finite_difference -= element.calculate_helmholtz_free_energy(
                        &nodal_coordinates
                    );
                    finite_difference/EPSILON
                }).collect()
            ).collect()
        }
        #[test]
        fn nodal_stiffnesses_deformed_symmetry()
        {
            let nodal_stiffness = get_nodal_stiffnesses(true, false);
            nodal_stiffness.iter()
            .enumerate()
            .for_each(|(a, nodal_stiffness_a)|
                nodal_stiffness_a.iter()
                .enumerate()
                .for_each(|(b, nodal_stiffness_ab)|
                    nodal_stiffness_ab.iter()
                    .enumerate()
                    .zip(nodal_stiffness_ab.transpose().iter())
                    .for_each(|((i, nodal_stiffness_ab_i), nodal_stiffness_ab_j)|
                        nodal_stiffness_ab_i.iter()
                        .enumerate()
                        .zip(nodal_stiffness_ab_j.iter())
                        .for_each(|((j, nodal_stiffness_ab_ij), nodal_stiffness_ab_ji)|
                            if a == b
                            {
                                assert_eq_within_tols(
                                    nodal_stiffness_ab_ij,
                                    &nodal_stiffness_ab_ji
                                )
                            }
                            else if i == j
                            {
                                assert_eq_within_tols(
                                    nodal_stiffness_ab_ij, &nodal_stiffness[b][a][i][j]
                                )
                            }
                        )
                    )
                )
            )
        }
        mod helmholtz_free_energy
        {
            use super::*;
            mod deformed
            {
                use super::*;
                #[test]
                fn finite_difference()
                {
                    get_nodal_forces(true, false, false).iter()
                    .zip(get_finite_difference_of_helmholtz_free_energy(true).iter())
                    .for_each(|(nodal_force, fd_nodal_force)|
                        nodal_force.iter()
                        .zip(fd_nodal_force.iter())
                        .for_each(|(nodal_force_i, fd_nodal_force_i)|
                            assert!(
                                (nodal_force_i/fd_nodal_force_i - 1.0).abs() < EPSILON ||
                                (nodal_force_i.abs() < EPSILON && fd_nodal_force_i.abs() < EPSILON)
                            )
                        )
                    )
                }
                #[test]
                fn minimized()
                {
                    let element = get_element();
                    let nodal_forces = get_nodal_forces(true, false, false);
                    let minimum = get_helmholtz_free_energy(true, false) - nodal_forces.dot(&get_coordinates());
                    let mut perturbed_coordinates = get_coordinates();
                    (0..N).for_each(|node|
                        (0..3).for_each(|i|{
                            perturbed_coordinates = get_coordinates();
                            perturbed_coordinates[node][i] += 0.5 * EPSILON;
                            assert!(
                                element.calculate_helmholtz_free_energy(
                                    &perturbed_coordinates
                                ) - nodal_forces.dot(
                                    &perturbed_coordinates
                                ) > minimum
                            );
                            perturbed_coordinates[node][i] -= EPSILON;
                            assert!(
                                element.calculate_helmholtz_free_energy(
                                    &perturbed_coordinates
                                ) - nodal_forces.dot(
                                    &perturbed_coordinates
                                ) > minimum
                            );
                        })
                    )
                }
                #[test]
                fn objectivity()
                {
                    assert_eq_within_tols(
                        &get_helmholtz_free_energy(true, false),
                        &get_helmholtz_free_energy(true, true)
                    )
                }
                #[test]
                fn positive()
                {
                    assert!(
                        get_helmholtz_free_energy(true, false) > 0.0
                    )
                }
            }
            mod undeformed
            {
                use super::*;
                #[test]
                fn finite_difference()
                {
                    get_finite_difference_of_helmholtz_free_energy(false).iter()
                    .for_each(|fd_nodal_force|
                        fd_nodal_force.iter()
                        .for_each(|fd_nodal_force_i|
                            assert!(
                                fd_nodal_force_i.abs() < EPSILON
                            )
                        )
                    )
                }
                #[test]
                fn minimized()
                {
                    let element = get_element();
                    let minimum = get_helmholtz_free_energy(false, false);
                    let mut perturbed_coordinates = get_reference_coordinates();
                    (0..N).for_each(|node|
                        (0..3).for_each(|i|{
                            perturbed_coordinates = get_reference_coordinates();
                            perturbed_coordinates[node][i] += 0.5 * EPSILON;
                            assert!(
                                element.calculate_helmholtz_free_energy(
                                    &perturbed_coordinates.convert()
                                ) > minimum
                            );
                            perturbed_coordinates[node][i] -= EPSILON;
                            assert!(
                                element.calculate_helmholtz_free_energy(
                                    &perturbed_coordinates.convert()
                                ) > minimum
                            );
                        })
                    )
                }
                #[test]
                fn objectivity()
                {
                    assert_eq_within_tols(
                        &get_helmholtz_free_energy(false, true), &0.0
                    )
                }
                #[test]
                fn zero()
                {
                    assert_eq!(
                        get_helmholtz_free_energy(false, false), 0.0
                    )
                }
            }
        }
    }
}
pub(crate) use test_helmholtz_free_energy;

macro_rules! test_finite_element_with_elastic_constitutive_model
{
    ($element: ident, $constitutive_model: ident, $constitutive_model_parameters: ident) =>
    {
        fn get_nodal_forces(is_deformed: bool, is_rotated: bool, _: bool) -> NodalForces<N>
        {
            if is_rotated
            {
                if is_deformed
                {
                    get_rotation_current_configuration().transpose() *
                    get_element_transformed().calculate_nodal_forces(
                        &get_coordinates_transformed()
                    )
                }
                else
                {
                    get_element().calculate_nodal_forces(
                        &get_reference_coordinates_transformed().convert()
                    )
                }
            }
            else
            {
                if is_deformed
                {
                    get_element().calculate_nodal_forces(
                        &get_coordinates()
                    )
                }
                else
                {
                    get_element().calculate_nodal_forces(
                        &get_reference_coordinates().convert()
                    )
                }
            }
        }
        fn get_nodal_stiffnesses(is_deformed: bool, is_rotated: bool) -> NodalStiffnesses<N>
        {
            if is_rotated
            {
                if is_deformed
                {
                    get_rotation_current_configuration().transpose() *
                    get_element_transformed().calculate_nodal_stiffnesses(
                        &get_coordinates_transformed()
                    ) * get_rotation_current_configuration()
                }
                else
                {
                    let converted: TensorRank2<3, 1, 1> = get_rotation_reference_configuration().convert();
                    converted.transpose() *
                    get_element_transformed().calculate_nodal_stiffnesses(
                        &get_reference_coordinates_transformed().convert()
                    ) * converted
                }
            }
            else
            {
                if is_deformed
                {
                    get_element().calculate_nodal_stiffnesses(
                        &get_coordinates()
                    )
                }
                else
                {
                    get_element().calculate_nodal_stiffnesses(
                        &get_reference_coordinates().convert()
                    )
                }
            }
        }
        fn get_finite_difference_of_nodal_forces(is_deformed: bool) -> NodalStiffnesses<N>
        {
            let element = get_element();
            let mut finite_difference = 0.0;
            (0..N).map(|node_a|
                (0..N).map(|node_b|
                    (0..3).map(|i|
                        (0..3).map(|j|{
                            let mut nodal_coordinates = 
                            if is_deformed
                            {
                                get_coordinates()
                            }
                            else
                            {
                                get_reference_coordinates().convert()
                            };
                            nodal_coordinates[node_a][i] += 0.5 * EPSILON;
                            finite_difference = element.calculate_nodal_forces(
                                &nodal_coordinates
                            )[node_b][j];
                            nodal_coordinates[node_a][i] -= EPSILON;
                            finite_difference -= element.calculate_nodal_forces(
                                &nodal_coordinates
                            )[node_b][j];
                            finite_difference/EPSILON
                        }).collect()
                    ).collect()
                ).collect()
            ).collect()
        }
        crate::fem::block::element::test::test_nodal_forces_and_nodal_stiffnesses!(
            $element, $constitutive_model, $constitutive_model_parameters
        );
    }
}
pub(crate) use test_finite_element_with_elastic_constitutive_model;

macro_rules! test_finite_element_with_hyperelastic_constitutive_model
{
    ($element: ident, $constitutive_model: ident, $constitutive_model_parameters: ident) =>
    {
        crate::fem::block::element::test::test_finite_element_with_elastic_constitutive_model!(
            $element, $constitutive_model, $constitutive_model_parameters
        );
        crate::fem::block::element::test::test_helmholtz_free_energy!(
            $element, $constitutive_model, $constitutive_model_parameters
        );
    }
}
pub(crate) use test_finite_element_with_hyperelastic_constitutive_model;

macro_rules! test_finite_element_with_viscoelastic_constitutive_model
{
    ($element: ident, $constitutive_model: ident, $constitutive_model_parameters: ident) =>
    {
        fn get_nodal_forces(is_deformed: bool, is_rotated: bool, is_xtra: bool) -> NodalForces<N>
        {
            if is_xtra
            {
                if is_rotated
                {
                    if is_deformed
                    {
                        get_rotation_current_configuration().transpose() *
                        get_element_transformed().calculate_nodal_forces(
                            &get_coordinates_transformed(), &get_velocities_transformed()
                        )
                    }
                    else
                    {
                        get_element().calculate_nodal_forces(
                            &get_reference_coordinates_transformed().convert(), &NodalVelocities::zero()
                        )
                    }
                }
                else
                {
                    if is_deformed
                    {
                        get_element().calculate_nodal_forces(
                            &get_coordinates(), &get_velocities()
                        )
                    }
                    else
                    {
                        get_element().calculate_nodal_forces(
                            &get_reference_coordinates().convert(), &NodalVelocities::zero()
                        )
                    }
                }
            }
            else
            {
                if is_rotated
                {
                    if is_deformed
                    {
                        get_rotation_current_configuration().transpose() *
                        get_element_transformed().calculate_nodal_forces(
                            &get_coordinates_transformed(), &NodalVelocities::zero()
                        )
                    }
                    else
                    {
                        get_element().calculate_nodal_forces(
                            &get_reference_coordinates_transformed().convert(), &NodalVelocities::zero()
                        )
                    }
                }
                else
                {
                    if is_deformed
                    {
                        get_element().calculate_nodal_forces(
                            &get_coordinates(), &NodalVelocities::zero()
                        )
                    }
                    else
                    {
                        get_element().calculate_nodal_forces(
                            &get_reference_coordinates().convert(), &NodalVelocities::zero()
                        )
                    }
                }
            }
        }
        fn get_nodal_stiffnesses(is_deformed: bool, is_rotated: bool) -> NodalStiffnesses<N>
        {
            if is_rotated
            {
                if is_deformed
                {
                    get_rotation_current_configuration().transpose() *
                    get_element_transformed().calculate_nodal_stiffnesses(
                        &get_coordinates_transformed(), &get_velocities_transformed()
                    ) * get_rotation_current_configuration()
                }
                else
                {
                    let converted: TensorRank2<3, 1, 1> = get_rotation_reference_configuration().convert();
                    converted.transpose() *
                    get_element_transformed().calculate_nodal_stiffnesses(
                        &get_reference_coordinates_transformed().convert(), &NodalVelocities::zero()
                    ) * converted
                }
            }
            else
            {
                if is_deformed
                {
                    get_element().calculate_nodal_stiffnesses(
                        &get_coordinates(), &get_velocities()
                    )
                }
                else
                {
                    get_element().calculate_nodal_stiffnesses(
                        &get_reference_coordinates().convert(), &NodalVelocities::zero()
                    )
                }
            }
        }
        fn get_finite_difference_of_nodal_forces(is_deformed: bool) -> NodalStiffnesses<N>
        {
            let element = get_element();
            let mut finite_difference = 0.0;
            (0..N).map(|node_a|
                (0..N).map(|node_b|
                    (0..3).map(|i|
                        (0..3).map(|j|{
                            let nodal_coordinates = 
                            if is_deformed
                            {
                                get_coordinates()
                            }
                            else
                            {
                                get_reference_coordinates().convert()
                            };
                            let mut nodal_velocities = 
                            if is_deformed
                            {
                                get_velocities()
                            }
                            else
                            {
                                NodalVelocities::zero()
                            };
                            nodal_velocities[node_a][i] += 0.5 * EPSILON;
                            finite_difference = element.calculate_nodal_forces(
                                &nodal_coordinates, &nodal_velocities
                            )[node_b][j];
                            nodal_velocities[node_a][i] -= EPSILON;
                            finite_difference -= element.calculate_nodal_forces(
                                &nodal_coordinates, &nodal_velocities
                            )[node_b][j];
                            finite_difference/EPSILON
                        }).collect()
                    ).collect()
                ).collect()
            ).collect()
        }
        crate::fem::block::element::test::test_nodal_forces_and_nodal_stiffnesses!(
            $element, $constitutive_model, $constitutive_model_parameters
        );
    }
}
pub(crate) use test_finite_element_with_viscoelastic_constitutive_model;

macro_rules! test_finite_element_with_elastic_hyperviscous_constitutive_model
{
    ($element: ident, $constitutive_model: ident, $constitutive_model_parameters: ident) =>
    {
        crate::fem::block::element::test::test_finite_element_with_viscoelastic_constitutive_model!(
            $element, $constitutive_model, $constitutive_model_parameters
        );
        fn get_viscous_dissipation(is_deformed: bool, is_rotated: bool) -> Scalar
        {
            if is_rotated
            {
                if is_deformed
                {
                    get_element_transformed()
                    .calculate_viscous_dissipation(
                        &get_coordinates_transformed(), &get_velocities_transformed()
                    )
                }
                else
                {
                    get_element_transformed()
                    .calculate_viscous_dissipation(
                        &get_reference_coordinates_transformed().convert(), &NodalVelocities::zero()
                    )
                }
            }
            else
            {
                if is_deformed
                {
                    get_element()
                    .calculate_viscous_dissipation(
                        &get_coordinates(), &get_velocities()
                    )
                }
                else
                {
                    get_element()
                    .calculate_viscous_dissipation(
                        &get_reference_coordinates().convert(), &NodalVelocities::zero()
                    )
                }
            }
        }
        fn get_dissipation_potential(is_deformed: bool, is_rotated: bool) -> Scalar
        {
            if is_rotated
            {
                if is_deformed
                {
                    get_element_transformed()
                    .calculate_dissipation_potential(
                        &get_coordinates_transformed(), &get_velocities_transformed()
                    )
                }
                else
                {
                    get_element_transformed()
                    .calculate_dissipation_potential(
                        &get_reference_coordinates_transformed().convert(), &NodalVelocities::zero()
                    )
                }
            }
            else
            {
                if is_deformed
                {
                    get_element()
                    .calculate_dissipation_potential(
                        &get_coordinates(), &get_velocities()
                    )
                }
                else
                {
                    get_element()
                    .calculate_dissipation_potential(
                        &get_reference_coordinates().convert(), &NodalVelocities::zero()
                    )
                }
            }
        }
        fn get_finite_difference_of_viscous_dissipation(is_deformed: bool) -> NodalForces<N>
        {
            let element = get_element();
            let mut finite_difference = 0.0;
            (0..N).map(|node|
                (0..3).map(|i|{
                    let nodal_coordinates = 
                    if is_deformed
                    {
                        get_coordinates()
                    }
                    else
                    {
                        get_reference_coordinates().convert()
                    };
                    let mut nodal_velocities = 
                    if is_deformed
                    {
                        get_velocities()
                    }
                    else
                    {
                        NodalVelocities::zero()
                    };
                    nodal_velocities[node][i] += 0.5 * EPSILON;
                    finite_difference = element.calculate_viscous_dissipation(
                        &nodal_coordinates, &nodal_velocities
                    );
                    nodal_velocities[node][i] -= EPSILON;
                    finite_difference -= element.calculate_viscous_dissipation(
                        &nodal_coordinates, &nodal_velocities
                    );
                    finite_difference/EPSILON
                }).collect()
            ).collect()
        }
        fn get_finite_difference_of_dissipation_potential(is_deformed: bool) -> NodalForces<N>
        {
            let element = get_element();
            let mut finite_difference = 0.0;
            (0..N).map(|node|
                (0..3).map(|i|{
                    let nodal_coordinates = 
                    if is_deformed
                    {
                        get_coordinates()
                    }
                    else
                    {
                        get_reference_coordinates().convert()
                    };
                    let mut nodal_velocities = 
                    if is_deformed
                    {
                        get_velocities()
                    }
                    else
                    {
                        NodalVelocities::zero()
                    };
                    nodal_velocities[node][i] += 0.5 * EPSILON;
                    finite_difference = element.calculate_dissipation_potential(
                        &nodal_coordinates, &nodal_velocities
                    );
                    nodal_velocities[node][i] -= EPSILON;
                    finite_difference -= element.calculate_dissipation_potential(
                        &nodal_coordinates, &nodal_velocities
                    );
                    finite_difference/EPSILON
                }).collect()
            ).collect()
        }
        mod viscous_dissipation
        {
            use super::*;
            mod deformed
            {
                use super::*;
                #[test]
                fn finite_difference()
                {
                    (get_nodal_forces(true, false, true) - get_nodal_forces(true, false, false)).iter()
                    .zip(get_finite_difference_of_viscous_dissipation(true).iter())
                    .for_each(|(nodal_force, fd_nodal_force)|
                        nodal_force.iter()
                        .zip(fd_nodal_force.iter())
                        .for_each(|(nodal_force_i, fd_nodal_force_i)|
                            assert!(
                                (nodal_force_i/fd_nodal_force_i - 1.0).abs() < EPSILON ||
                                (nodal_force_i.abs() < EPSILON && fd_nodal_force_i.abs() < EPSILON)
                            )
                        )
                    )
                }
                #[test]
                fn minimized()
                {
                    let element = get_element();
                    let nodal_forces = get_nodal_forces(true, false, true) - get_nodal_forces(true, false, false);
                    let minimum = get_viscous_dissipation(true, false) - nodal_forces.dot(&get_velocities());
                    let mut perturbed_velocities = get_velocities();
                    (0..N).for_each(|node|
                        (0..3).for_each(|i|{
                            perturbed_velocities = get_velocities();
                            perturbed_velocities[node][i] += 0.5 * EPSILON;
                            assert!(
                                element.calculate_viscous_dissipation(
                                    &get_coordinates(), &perturbed_velocities
                                ) - nodal_forces.dot(
                                    &perturbed_velocities
                                ) > minimum
                            );
                            perturbed_velocities[node][i] -= EPSILON;
                            assert!(
                                element.calculate_viscous_dissipation(
                                    &get_coordinates(), &perturbed_velocities
                                ) - nodal_forces.dot(
                                    &perturbed_velocities
                                ) > minimum
                            );
                        })
                    )
                }
                #[test]
                fn objectivity()
                {
                    assert_eq_within_tols(
                        &get_viscous_dissipation(true, false),
                        &get_viscous_dissipation(true, true)
                    )
                }
                #[test]
                fn positive()
                {
                    assert!(
                        get_viscous_dissipation(true, false) > 0.0
                    )
                }
            }
            mod undeformed
            {
                use super::*;
                #[test]
                fn finite_difference()
                {
                    get_finite_difference_of_viscous_dissipation(false).iter()
                    .for_each(|fd_nodal_force|
                        fd_nodal_force.iter()
                        .for_each(|fd_nodal_force_i|
                            assert!(
                                fd_nodal_force_i.abs() < EPSILON
                            )
                        )
                    )
                }
                #[test]
                fn minimized()
                {
                    let element = get_element();
                    let minimum = get_viscous_dissipation(false, false);
                    let mut perturbed_velocities = NodalVelocities::zero();
                    (0..N).for_each(|node|
                        (0..3).for_each(|i|{
                            perturbed_velocities = NodalVelocities::zero();
                            perturbed_velocities[node][i] += 0.5 * EPSILON;
                            assert!(
                                element.calculate_viscous_dissipation(
                                    &get_reference_coordinates().convert(), &perturbed_velocities
                                ) > minimum
                            );
                            perturbed_velocities[node][i] -= EPSILON;
                            assert!(
                                element.calculate_viscous_dissipation(
                                    &get_reference_coordinates().convert(), &perturbed_velocities
                                ) > minimum
                            );
                        })
                    )
                }
                #[test]
                fn objectivity()
                {
                    assert_eq_within_tols(
                        &get_viscous_dissipation(false, true), &0.0
                    )
                }
                #[test]
                fn zero()
                {
                    assert_eq!(
                        get_viscous_dissipation(false, false), 0.0
                    )
                }
            }
        }
        mod dissipation_potential
        {
            use super::*;
            mod deformed
            {
                use super::*;
                #[test]
                fn finite_difference()
                {
                    get_nodal_forces(true, false, true).iter()
                    .zip(get_finite_difference_of_dissipation_potential(true).iter())
                    .for_each(|(nodal_force, fd_nodal_force)|
                        nodal_force.iter()
                        .zip(fd_nodal_force.iter())
                        .for_each(|(nodal_force_i, fd_nodal_force_i)|
                            assert!(
                                (nodal_force_i/fd_nodal_force_i - 1.0).abs() < EPSILON ||
                                (nodal_force_i.abs() < EPSILON && fd_nodal_force_i.abs() < EPSILON)
                            )
                        )
                    )
                }
                #[test]
                fn minimized()
                {
                    let element = get_element();
                    let nodal_forces = get_nodal_forces(true, false, true);
                    let minimum = get_dissipation_potential(true, false) - nodal_forces.dot(&get_velocities());
                    let mut perturbed_velocities = get_velocities();
                    (0..N).for_each(|node|
                        (0..3).for_each(|i|{
                            perturbed_velocities = get_velocities();
                            perturbed_velocities[node][i] += 0.5 * EPSILON;
                            assert!(
                                element.calculate_dissipation_potential(
                                    &get_coordinates(), &perturbed_velocities
                                ) - nodal_forces.dot(
                                    &perturbed_velocities
                                ) > minimum
                            );
                            perturbed_velocities[node][i] -= EPSILON;
                            assert!(
                                element.calculate_dissipation_potential(
                                    &get_coordinates(), &perturbed_velocities
                                ) - nodal_forces.dot(
                                    &perturbed_velocities
                                ) > minimum
                            );
                        })
                    )
                }
                #[test]
                fn objectivity()
                {
                    assert_eq_within_tols(
                        &get_dissipation_potential(true, false),
                        &get_dissipation_potential(true, true)
                    )
                }
            }
            mod undeformed
            {
                use super::*;
                #[test]
                fn finite_difference()
                {
                    get_finite_difference_of_dissipation_potential(false).iter()
                    .for_each(|fd_nodal_force|
                        fd_nodal_force.iter()
                        .for_each(|fd_nodal_force_i|
                            assert!(
                                fd_nodal_force_i.abs() < EPSILON
                            )
                        )
                    )
                }
                #[test]
                fn minimized()
                {
                    let element = get_element();
                    let minimum = get_dissipation_potential(false, false);
                    let mut perturbed_velocities = NodalVelocities::zero();
                    (0..N).for_each(|node|
                        (0..3).for_each(|i|{
                            perturbed_velocities = NodalVelocities::zero();
                            perturbed_velocities[node][i] += 0.5 * EPSILON;
                            assert!(
                                element.calculate_dissipation_potential(
                                    &get_reference_coordinates().convert(), &perturbed_velocities
                                ) > minimum
                            );
                            perturbed_velocities[node][i] -= EPSILON;
                            assert!(
                                element.calculate_dissipation_potential(
                                    &get_reference_coordinates().convert(), &perturbed_velocities
                                ) > minimum
                            );
                        })
                    )
                }
                #[test]
                fn objectivity()
                {
                    assert_eq_within_tols(
                        &get_dissipation_potential(false, true), &0.0
                    )
                }
                #[test]
                fn zero()
                {
                    assert_eq!(
                        get_dissipation_potential(false, false), 0.0
                    )
                }
            }
        }
    }
}
pub(crate) use test_finite_element_with_elastic_hyperviscous_constitutive_model;

macro_rules! test_finite_element_with_hyperviscoelastic_constitutive_model
{
    ($element: ident, $constitutive_model: ident, $constitutive_model_parameters: ident) =>
    {
        crate::fem::block::element::test::test_finite_element_with_elastic_hyperviscous_constitutive_model!(
            $element, $constitutive_model, $constitutive_model_parameters
        );
        crate::fem::block::element::test::test_helmholtz_free_energy!(
            $element, $constitutive_model, $constitutive_model_parameters
        );
        #[test]
        fn dissipation_potential_deformed_positive()
        {
            assert!(
                get_dissipation_potential(true, false) > 0.0
            )
        }
    }
}
pub(crate) use test_finite_element_with_hyperviscoelastic_constitutive_model;