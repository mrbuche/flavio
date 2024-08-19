use crate::{
    constitutive::solid::elastic::test::ALMANSIHAMELPARAMETERS as ALMANSIHAMELPARAMETERSELASTIC,
    mechanics::Scalar,
};
pub const ALMANSIHAMELPARAMETERS: &[Scalar; 4] = &[
    ALMANSIHAMELPARAMETERSELASTIC[0],
    ALMANSIHAMELPARAMETERSELASTIC[1],
    11.0,
    1.0,
];

macro_rules! calculate_viscous_dissipation_from_deformation_gradient_rate_simple {
    ($constitutive_model_constructed: expr, $deformation_gradient_rate: expr) => {
        $constitutive_model_constructed
            .calculate_viscous_dissipation(
                &DeformationGradient::identity(),
                $deformation_gradient_rate,
            )
            .expect("the unexpected")
    };
}
pub(crate) use calculate_viscous_dissipation_from_deformation_gradient_rate_simple;

macro_rules! calculate_viscous_dissipation_from_deformation_gradient_and_deformation_gradient_rate {
    ($constitutive_model_constructed: expr, $deformation_gradient: expr, $deformation_gradient_rate: expr) => {
        $constitutive_model_constructed
            .calculate_viscous_dissipation($deformation_gradient, $deformation_gradient_rate)
            .expect("the unexpected")
    };
}
pub(crate) use calculate_viscous_dissipation_from_deformation_gradient_and_deformation_gradient_rate;

macro_rules! calculate_dissipation_potential_from_deformation_gradient_and_deformation_gradient_rate {
    ($constitutive_model_constructed: expr, $deformation_gradient: expr, $deformation_gradient_rate: expr) => {
        $constitutive_model_constructed
            .calculate_dissipation_potential($deformation_gradient, $deformation_gradient_rate)
            .expect("the unexpected")
    };
}
pub(crate) use calculate_dissipation_potential_from_deformation_gradient_and_deformation_gradient_rate;

macro_rules! use_viscoelastic_macros
{
    () =>
    {
        use crate::constitutive::solid::viscoelastic::test::
        {
            calculate_cauchy_stress_from_deformation_gradient,
            calculate_cauchy_stress_from_deformation_gradient_simple,
            calculate_cauchy_stress_from_deformation_gradient_rotated,
            calculate_cauchy_stress_from_deformation_gradient_and_deformation_gradient_rate,
            calculate_cauchy_rate_tangent_stiffness_from_deformation_gradient_and_deformation_gradient_rate,
            calculate_first_piola_kirchoff_stress_from_deformation_gradient,
            calculate_first_piola_kirchoff_stress_from_deformation_gradient_simple,
            calculate_first_piola_kirchoff_stress_from_deformation_gradient_rotated,
            calculate_first_piola_kirchoff_stress_from_deformation_gradient_rate_simple,
            calculate_first_piola_kirchoff_stress_from_deformation_gradient_and_deformation_gradient_rate,
            calculate_first_piola_kirchoff_rate_tangent_stiffness_from_deformation_gradient_and_deformation_gradient_rate,
            calculate_second_piola_kirchoff_stress_from_deformation_gradient,
            calculate_second_piola_kirchoff_stress_from_deformation_gradient_simple,
            calculate_second_piola_kirchoff_stress_from_deformation_gradient_rotated,
            calculate_second_piola_kirchoff_stress_from_deformation_gradient_and_deformation_gradient_rate,
            calculate_second_piola_kirchoff_rate_tangent_stiffness_from_deformation_gradient_and_deformation_gradient_rate
        };
    }
}
pub(crate) use use_viscoelastic_macros;

macro_rules! test_solid_elastic_hyperviscous_constitutive_model
{
    ($constitutive_model: ident, $constitutive_model_parameters: expr, $constitutive_model_constructed: expr) =>
    {
        crate::constitutive::solid::elastic::test::test_solid_constitutive_construction!(
            $constitutive_model, $constitutive_model_parameters, $constitutive_model_constructed
        );
        crate::constitutive::solid::elastic::test::test_solid_constitutive_model_no_tangents!(
            $constitutive_model_constructed
        );
        crate::constitutive::solid::viscoelastic::test::test_solid_viscous_constitutive_model!(
            $constitutive_model, $constitutive_model_parameters, $constitutive_model_constructed
        );
        crate::constitutive::solid::elastic_hyperviscous::test::test_solid_elastic_hyperviscous_specifics!(
            $constitutive_model, $constitutive_model_parameters, $constitutive_model_constructed
        );
    }
}
pub(crate) use test_solid_elastic_hyperviscous_constitutive_model;

macro_rules! test_solid_elastic_hyperviscous_specifics
{
    ($constitutive_model: ident, $constitutive_model_parameters: expr, $constitutive_model_constructed: expr) =>
    {
        mod elastic_hyperviscous
        {
            use super::*;
            mod viscous_dissipation // eventually should go in fluid/hyperviscous/test.rs
            {
                use super::*;
                fn calculate_first_piola_kirchoff_stress_from_finite_difference_of_viscous_dissipation(is_deformed: bool) -> FirstPiolaKirchoffStress
                {
                    let mut first_piola_kirchoff_stress = FirstPiolaKirchoffStress::zero();
                    for i in 0..3
                    {
                        for j in 0..3
                        {
                            let mut deformation_gradient_rate_plus =
                                if is_deformed
                                {
                                    get_deformation_gradient_rate()
                                }
                                else
                                {
                                    DeformationGradientRate::zero()
                                };
                            deformation_gradient_rate_plus[i][j] += 0.5*EPSILON;
                            let helmholtz_free_energy_density_plus =
                            calculate_viscous_dissipation_from_deformation_gradient_rate_simple!(
                                $constitutive_model_constructed, &deformation_gradient_rate_plus
                            );
                            let mut deformation_gradient_rate_minus =
                                if is_deformed
                                {
                                    get_deformation_gradient_rate()
                                }
                                else
                                {
                                    DeformationGradientRate::zero()
                                };
                            deformation_gradient_rate_minus[i][j] -= 0.5*EPSILON;
                            let helmholtz_free_energy_density_minus =
                            calculate_viscous_dissipation_from_deformation_gradient_rate_simple!(
                                $constitutive_model_constructed, &deformation_gradient_rate_minus
                            );
                            first_piola_kirchoff_stress[i][j] = (
                                helmholtz_free_energy_density_plus - helmholtz_free_energy_density_minus
                            )/EPSILON;
                        }
                    }
                    first_piola_kirchoff_stress
                }
                mod deformed
                {
                    use super::*;
                    #[test]
                    fn finite_difference()
                    {
                        calculate_first_piola_kirchoff_stress_from_deformation_gradient_rate_simple!(
                            $constitutive_model_constructed, &get_deformation_gradient_rate()
                        ).iter().zip(
                            calculate_first_piola_kirchoff_stress_from_finite_difference_of_viscous_dissipation(true).iter()
                        ).for_each(|(first_piola_kirchoff_stress_i, fd_first_piola_kirchoff_stress_i)|
                            first_piola_kirchoff_stress_i.iter()
                            .zip(fd_first_piola_kirchoff_stress_i.iter())
                            .for_each(|(first_piola_kirchoff_stress_ij, fd_first_piola_kirchoff_stress_ij)|
                                assert!((first_piola_kirchoff_stress_ij/fd_first_piola_kirchoff_stress_ij - 1.0).abs() < EPSILON)
                            )
                        )
                    }
                    #[test]
                    fn minimized()
                    {
                        let first_piola_kirchoff_stress =
                        calculate_first_piola_kirchoff_stress_from_deformation_gradient_rate_simple!(
                            $constitutive_model_constructed, &get_deformation_gradient_rate()
                        );
                        let minimum =
                        calculate_viscous_dissipation_from_deformation_gradient_rate_simple!(
                            $constitutive_model_constructed, &get_deformation_gradient_rate()
                        ) - first_piola_kirchoff_stress.full_contraction(
                            &get_deformation_gradient_rate()
                        );
                        let mut perturbed_deformation_gradient_rate = get_deformation_gradient_rate();
                        (0..3).for_each(|i|
                            (0..3).for_each(|j|{
                                perturbed_deformation_gradient_rate = get_deformation_gradient_rate();
                                perturbed_deformation_gradient_rate[i][j] += 0.5 * EPSILON;
                                assert!(
                                    calculate_viscous_dissipation_from_deformation_gradient_rate_simple!(
                                        $constitutive_model_constructed, &perturbed_deformation_gradient_rate
                                    ) - first_piola_kirchoff_stress.full_contraction(
                                        &perturbed_deformation_gradient_rate
                                    ) > minimum
                                );
                                perturbed_deformation_gradient_rate[i][j] -= EPSILON;
                                assert!(
                                    calculate_viscous_dissipation_from_deformation_gradient_rate_simple!(
                                        $constitutive_model_constructed, &perturbed_deformation_gradient_rate
                                    ) - first_piola_kirchoff_stress.full_contraction(
                                        &perturbed_deformation_gradient_rate
                                    ) > minimum
                                );
                            })
                        )
                    }
                    #[test]
                    fn objectivity()
                    {
                        assert_eq_within_tols(
                            &calculate_viscous_dissipation_from_deformation_gradient_and_deformation_gradient_rate!(
                                $constitutive_model_constructed, &get_deformation_gradient(), &get_deformation_gradient_rate()
                            ),
                            &calculate_viscous_dissipation_from_deformation_gradient_and_deformation_gradient_rate!(
                                $constitutive_model_constructed, &get_deformation_gradient_rotated(), &get_deformation_gradient_rate_rotated()
                            )
                        )
                    }
                    #[test]
                    fn positive()
                    {
                        assert!(
                            calculate_viscous_dissipation_from_deformation_gradient_rate_simple!(
                                $constitutive_model_constructed,  &get_deformation_gradient_rate()
                            ) > 0.0
                        )
                    }
                }
                mod undeformed
                {
                    use super::*;
                    #[test]
                    fn finite_difference()
                    {
                        calculate_first_piola_kirchoff_stress_from_finite_difference_of_viscous_dissipation(false).iter()
                        .for_each(|fd_first_piola_kirchoff_stress_i|
                            fd_first_piola_kirchoff_stress_i.iter()
                            .for_each(|fd_first_piola_kirchoff_stress_ij|
                                assert!(fd_first_piola_kirchoff_stress_ij.abs() < EPSILON)
                            )
                        )
                    }
                    #[test]
                    fn minimized()
                    {
                        let minimum =
                        calculate_viscous_dissipation_from_deformation_gradient_rate_simple!(
                            $constitutive_model_constructed, &DeformationGradientRate::zero()
                        );
                        let mut perturbed_deformation_gradient_rate = DeformationGradientRate::zero();
                        (0..3).for_each(|i|
                            (0..3).for_each(|j|{
                                perturbed_deformation_gradient_rate = DeformationGradientRate::zero();
                                perturbed_deformation_gradient_rate[i][j] += 0.5 * EPSILON;
                                assert!(
                                    calculate_viscous_dissipation_from_deformation_gradient_rate_simple!(
                                        $constitutive_model_constructed, &perturbed_deformation_gradient_rate
                                    ) > minimum
                                );
                                perturbed_deformation_gradient_rate[i][j] -= EPSILON;
                                assert!(
                                    calculate_viscous_dissipation_from_deformation_gradient_rate_simple!(
                                        $constitutive_model_constructed, &perturbed_deformation_gradient_rate
                                    ) > minimum
                                );
                            })
                        )
                    }
                    #[test]
                    fn zero()
                    {
                        assert_eq!(
                            calculate_viscous_dissipation_from_deformation_gradient_rate_simple!(
                                $constitutive_model_constructed,  &DeformationGradientRate::zero()
                            ), 0.0
                        )
                    }
                }
            }
            mod dissipation_potential
            {
                use super::*;
                fn calculate_first_piola_kirchoff_stress_from_finite_difference_of_dissipation_potential(is_deformed: bool) -> FirstPiolaKirchoffStress
                {
                    let deformation_gradient =
                        if is_deformed
                        {
                            get_deformation_gradient()
                        }
                        else
                        {
                            DeformationGradient::identity()
                        };
                    let mut first_piola_kirchoff_stress = FirstPiolaKirchoffStress::zero();
                    for i in 0..3
                    {
                        for j in 0..3
                        {
                            let mut deformation_gradient_rate_plus =
                                if is_deformed
                                {
                                    get_deformation_gradient_rate()
                                }
                                else
                                {
                                    DeformationGradientRate::zero()
                                };
                            deformation_gradient_rate_plus[i][j] += 0.5*EPSILON;
                            let helmholtz_free_energy_density_plus =
                            calculate_dissipation_potential_from_deformation_gradient_and_deformation_gradient_rate!(
                                $constitutive_model_constructed, &deformation_gradient, &deformation_gradient_rate_plus
                            );
                            let mut deformation_gradient_rate_minus =
                                if is_deformed
                                {
                                    get_deformation_gradient_rate()
                                }
                                else
                                {
                                    DeformationGradientRate::zero()
                                };
                            deformation_gradient_rate_minus[i][j] -= 0.5*EPSILON;
                            let helmholtz_free_energy_density_minus =
                            calculate_dissipation_potential_from_deformation_gradient_and_deformation_gradient_rate!(
                                $constitutive_model_constructed, &deformation_gradient, &deformation_gradient_rate_minus
                            );
                            first_piola_kirchoff_stress[i][j] = (
                                helmholtz_free_energy_density_plus - helmholtz_free_energy_density_minus
                            )/EPSILON;
                        }
                    }
                    first_piola_kirchoff_stress
                }
                mod deformed
                {
                    use super::*;
                    #[test]
                    fn finite_difference()
                    {
                        calculate_first_piola_kirchoff_stress_from_deformation_gradient_and_deformation_gradient_rate!(
                            $constitutive_model_constructed, &get_deformation_gradient(), &get_deformation_gradient_rate()
                        ).iter().zip(
                            calculate_first_piola_kirchoff_stress_from_finite_difference_of_dissipation_potential(true).iter()
                        ).for_each(|(first_piola_kirchoff_stress_i, fd_first_piola_kirchoff_stress_i)|
                            first_piola_kirchoff_stress_i.iter()
                            .zip(fd_first_piola_kirchoff_stress_i.iter())
                            .for_each(|(first_piola_kirchoff_stress_ij, fd_first_piola_kirchoff_stress_ij)|
                                assert!((first_piola_kirchoff_stress_ij/fd_first_piola_kirchoff_stress_ij - 1.0).abs() < EPSILON)
                            )
                        )
                    }
                    #[test]
                    fn minimized()
                    {
                        let first_piola_kirchoff_stress =
                        calculate_first_piola_kirchoff_stress_from_deformation_gradient_and_deformation_gradient_rate!(
                            $constitutive_model_constructed, &get_deformation_gradient(), &get_deformation_gradient_rate()
                        );
                        let minimum =
                        calculate_dissipation_potential_from_deformation_gradient_and_deformation_gradient_rate!(
                            $constitutive_model_constructed, &get_deformation_gradient(), &get_deformation_gradient_rate()
                        ) - first_piola_kirchoff_stress.full_contraction(
                            &get_deformation_gradient_rate()
                        );
                        let mut perturbed_deformation_gradient_rate = get_deformation_gradient_rate();
                        (0..3).for_each(|i|
                            (0..3).for_each(|j|{
                                perturbed_deformation_gradient_rate = get_deformation_gradient_rate();
                                perturbed_deformation_gradient_rate[i][j] += 0.5 * EPSILON;
                                assert!(
                                    calculate_dissipation_potential_from_deformation_gradient_and_deformation_gradient_rate!(
                                        $constitutive_model_constructed, &get_deformation_gradient(), &perturbed_deformation_gradient_rate
                                    ) - first_piola_kirchoff_stress.full_contraction(
                                        &perturbed_deformation_gradient_rate
                                    ) > minimum
                                );
                                perturbed_deformation_gradient_rate[i][j] -= EPSILON;
                                assert!(
                                    calculate_dissipation_potential_from_deformation_gradient_and_deformation_gradient_rate!(
                                        $constitutive_model_constructed, &get_deformation_gradient(), &perturbed_deformation_gradient_rate
                                    ) - first_piola_kirchoff_stress.full_contraction(
                                        &perturbed_deformation_gradient_rate
                                    ) > minimum
                                );
                            })
                        )
                    }
                    #[test]
                    fn objectivity()
                    {
                        assert_eq_within_tols(
                            &calculate_dissipation_potential_from_deformation_gradient_and_deformation_gradient_rate!(
                                $constitutive_model_constructed, &get_deformation_gradient(), &get_deformation_gradient_rate()
                            ),
                            &calculate_dissipation_potential_from_deformation_gradient_and_deformation_gradient_rate!(
                                $constitutive_model_constructed, &get_deformation_gradient_rotated(), &get_deformation_gradient_rate_rotated()
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
                        calculate_first_piola_kirchoff_stress_from_finite_difference_of_dissipation_potential(false).iter()
                        .for_each(|fd_first_piola_kirchoff_stress_i|
                            fd_first_piola_kirchoff_stress_i.iter()
                            .for_each(|fd_first_piola_kirchoff_stress_ij|
                                assert!(fd_first_piola_kirchoff_stress_ij.abs() < EPSILON)
                            )
                        )
                    }
                    #[test]
                    fn minimized()
                    {
                        let minimum =
                        calculate_dissipation_potential_from_deformation_gradient_and_deformation_gradient_rate!(
                            $constitutive_model_constructed, &DeformationGradient::identity(), &DeformationGradientRate::zero()
                        );
                        let mut perturbed_deformation_gradient_rate = DeformationGradientRate::zero();
                        (0..3).for_each(|i|
                            (0..3).for_each(|j|{
                                perturbed_deformation_gradient_rate = DeformationGradientRate::zero();
                                perturbed_deformation_gradient_rate[i][j] += 0.5 * EPSILON;
                                assert!(
                                    calculate_dissipation_potential_from_deformation_gradient_and_deformation_gradient_rate!(
                                        $constitutive_model_constructed, &DeformationGradient::identity(), &perturbed_deformation_gradient_rate
                                    ) > minimum
                                );
                                perturbed_deformation_gradient_rate[i][j] -= EPSILON;
                                assert!(
                                    calculate_dissipation_potential_from_deformation_gradient_and_deformation_gradient_rate!(
                                        $constitutive_model_constructed, &DeformationGradient::identity(), &perturbed_deformation_gradient_rate
                                    ) > minimum
                                );
                            })
                        )
                    }
                    #[test]
                    fn zero()
                    {
                        assert_eq!(
                            calculate_dissipation_potential_from_deformation_gradient_and_deformation_gradient_rate!(
                                $constitutive_model_constructed, &DeformationGradient::identity(), &DeformationGradientRate::zero()
                            ), 0.0
                        )
                    }
                }
            }
            mod first_piola_kirchoff_rate_tangent_stiffness
            {
                use super::*;
                mod deformed
                {
                    use super::*;
                    #[test]
                    fn symmetry()
                    {
                        let first_piola_kirchoff_rate_tangent_stiffness =
                        calculate_first_piola_kirchoff_rate_tangent_stiffness_from_deformation_gradient_and_deformation_gradient_rate!(
                            $constitutive_model_constructed, &get_deformation_gradient(), &get_deformation_gradient_rate()
                        );
                        first_piola_kirchoff_rate_tangent_stiffness.iter().enumerate()
                        .for_each(|(i, first_piola_kirchoff_rate_tangent_stiffness_i)|
                            first_piola_kirchoff_rate_tangent_stiffness_i.iter().enumerate()
                            .for_each(|(j, first_piola_kirchoff_rate_tangent_stiffness_ij)|
                                first_piola_kirchoff_rate_tangent_stiffness_ij.iter()
                                .zip(first_piola_kirchoff_rate_tangent_stiffness.iter())
                                .for_each(|(first_piola_kirchoff_rate_tangent_stiffness_ijk, first_piola_kirchoff_rate_tangent_stiffness_k)|
                                    first_piola_kirchoff_rate_tangent_stiffness_ijk.iter()
                                    .zip(first_piola_kirchoff_rate_tangent_stiffness_k.iter())
                                    .for_each(|(first_piola_kirchoff_rate_tangent_stiffness_ijkl, first_piola_kirchoff_rate_tangent_stiffness_kl)|
                                        assert_eq_within_tols(
                                            first_piola_kirchoff_rate_tangent_stiffness_ijkl, &first_piola_kirchoff_rate_tangent_stiffness_kl[i][j]
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
                    fn symmetry()
                    {
                        let first_piola_kirchoff_rate_tangent_stiffness =
                        calculate_first_piola_kirchoff_rate_tangent_stiffness_from_deformation_gradient_and_deformation_gradient_rate!(
                            $constitutive_model_constructed, &DeformationGradient::identity(), &DeformationGradientRate::zero()
                        );
                        first_piola_kirchoff_rate_tangent_stiffness.iter().enumerate()
                        .for_each(|(i, first_piola_kirchoff_rate_tangent_stiffness_i)|
                            first_piola_kirchoff_rate_tangent_stiffness_i.iter().enumerate()
                            .for_each(|(j, first_piola_kirchoff_rate_tangent_stiffness_ij)|
                                first_piola_kirchoff_rate_tangent_stiffness_ij.iter()
                                .zip(first_piola_kirchoff_rate_tangent_stiffness.iter())
                                .for_each(|(first_piola_kirchoff_rate_tangent_stiffness_ijk, first_piola_kirchoff_rate_tangent_stiffness_k)|
                                    first_piola_kirchoff_rate_tangent_stiffness_ijk.iter()
                                    .zip(first_piola_kirchoff_rate_tangent_stiffness_k.iter())
                                    .for_each(|(first_piola_kirchoff_rate_tangent_stiffness_ijkl, first_piola_kirchoff_rate_tangent_stiffness_kl)|
                                        assert_eq_within_tols(
                                            first_piola_kirchoff_rate_tangent_stiffness_ijkl, &first_piola_kirchoff_rate_tangent_stiffness_kl[i][j]
                                        )
                                    )
                                )
                            )
                        )
                    }
                }
            }
        }
    }
}
pub(crate) use test_solid_elastic_hyperviscous_specifics;
