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

macro_rules! test_solve {
    ($constitutive_model_constructed: expr) => {
        #[test]
        fn temporary() -> Result<(), crate::math::test::TestError> {
            let rate = 0.1;
            let function = |_| rate;
            let (deformation_gradients, cauchy_stresses) =
                $constitutive_model_constructed.solve_uniaxial(function, [0.0, 0.2, 0.5])?;
            deformation_gradients
                .iter()
                .for_each(|deformation_gradient| {
                    println!();
                    deformation_gradient
                        .iter()
                        .for_each(|row| println!("{:?}", (row[0], row[1], row[2])))
                });
            println!();
            cauchy_stresses.iter().for_each(|cauchy_stress| {
                println!();
                cauchy_stress
                    .iter()
                    .for_each(|row| println!("{:?}", (row[0], row[1], row[2])))
            });
            todo!()
        }
        #[test]
        fn solve_uniaxial_compression() -> Result<(), crate::math::test::TestError> {
            let (deformation_gradient_rate, cauchy_stress) = $constitutive_model_constructed
                .solve_uniaxial_inner_inner(&DeformationGradient::identity(), &-4.4)?;
            assert!(cauchy_stress[0][0] < 0.0);
            crate::math::test::assert_eq_within_tols(
                &(cauchy_stress[1][1] / cauchy_stress[0][0]),
                &0.0,
            )?;
            crate::math::test::assert_eq_within_tols(
                &(cauchy_stress[2][2] / cauchy_stress[0][0]),
                &0.0,
            )?;
            assert!(cauchy_stress.is_diagonal());
            crate::math::test::assert_eq(
                &deformation_gradient_rate[1][1],
                &deformation_gradient_rate[2][2],
            )?;
            assert!(deformation_gradient_rate.is_diagonal());
            Ok(())
        }
        #[test]
        fn solve_uniaxial_tension() -> Result<(), crate::math::test::TestError> {
            let (deformation_gradient_rate, cauchy_stress) = $constitutive_model_constructed
                .solve_uniaxial_inner_inner(&DeformationGradient::identity(), &4.4)?;
            assert!(cauchy_stress[0][0] > 0.0);
            crate::math::test::assert_eq_within_tols(
                &(cauchy_stress[1][1] / cauchy_stress[0][0]),
                &0.0,
            )?;
            crate::math::test::assert_eq_within_tols(
                &(cauchy_stress[2][2] / cauchy_stress[0][0]),
                &0.0,
            )?;
            assert!(cauchy_stress.is_diagonal());
            crate::math::test::assert_eq(
                &deformation_gradient_rate[1][1],
                &deformation_gradient_rate[2][2],
            )?;
            assert!(deformation_gradient_rate.is_diagonal());
            Ok(())
        }
        #[test]
        fn solve_uniaxial_undeformed() -> Result<(), crate::math::test::TestError> {
            let (deformation_gradient_rate, cauchy_stress) = $constitutive_model_constructed
                .solve_uniaxial_inner_inner(&DeformationGradient::identity(), &0.0)?;
            assert!(cauchy_stress.is_zero());
            assert!(deformation_gradient_rate.is_zero());
            Ok(())
        }
    };
}
pub(crate) use test_solve;

macro_rules! calculate_viscous_dissipation_from_deformation_gradient_rate_simple {
    ($constitutive_model_constructed: expr, $deformation_gradient_rate: expr) => {
        $constitutive_model_constructed.calculate_viscous_dissipation(
            &DeformationGradient::identity(),
            $deformation_gradient_rate,
        )
    };
}
pub(crate) use calculate_viscous_dissipation_from_deformation_gradient_rate_simple;

macro_rules! calculate_viscous_dissipation_from_deformation_gradient_and_deformation_gradient_rate {
    ($constitutive_model_constructed: expr, $deformation_gradient: expr, $deformation_gradient_rate: expr) => {
        $constitutive_model_constructed
            .calculate_viscous_dissipation($deformation_gradient, $deformation_gradient_rate)
    };
}
pub(crate) use calculate_viscous_dissipation_from_deformation_gradient_and_deformation_gradient_rate;

macro_rules! calculate_dissipation_potential_from_deformation_gradient_and_deformation_gradient_rate {
    ($constitutive_model_constructed: expr, $deformation_gradient: expr, $deformation_gradient_rate: expr) => {
        $constitutive_model_constructed
            .calculate_dissipation_potential($deformation_gradient, $deformation_gradient_rate)
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
            calculate_second_piola_kirchoff_rate_tangent_stiffness_from_deformation_gradient_and_deformation_gradient_rate,
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
                fn calculate_first_piola_kirchoff_stress_from_finite_difference_of_viscous_dissipation(is_deformed: bool) ->  Result<FirstPiolaKirchoffStress, TestError>
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
                            )?;
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
                            )?;
                            first_piola_kirchoff_stress[i][j] = (
                                helmholtz_free_energy_density_plus - helmholtz_free_energy_density_minus
                            )/EPSILON;
                        }
                    }
                    Ok(first_piola_kirchoff_stress)
                }
                mod deformed
                {
                    use super::*;
                    #[test]
                    fn finite_difference() -> Result<(), TestError>
                    {
                        calculate_first_piola_kirchoff_stress_from_deformation_gradient_rate_simple!(
                            $constitutive_model_constructed, &get_deformation_gradient_rate()
                        )?.iter().zip(
                            calculate_first_piola_kirchoff_stress_from_finite_difference_of_viscous_dissipation(true)?.iter()
                        ).for_each(|(first_piola_kirchoff_stress_i, fd_first_piola_kirchoff_stress_i)|
                            first_piola_kirchoff_stress_i.iter()
                            .zip(fd_first_piola_kirchoff_stress_i.iter())
                            .for_each(|(first_piola_kirchoff_stress_ij, fd_first_piola_kirchoff_stress_ij)|
                                assert!((first_piola_kirchoff_stress_ij/fd_first_piola_kirchoff_stress_ij - 1.0).abs() < EPSILON)
                            )
                        );
                        Ok(())
                    }
                    #[test]
                    fn minimized() -> Result<(), TestError>
                    {
                        let first_piola_kirchoff_stress =
                        calculate_first_piola_kirchoff_stress_from_deformation_gradient_rate_simple!(
                            $constitutive_model_constructed, &get_deformation_gradient_rate()
                        )?;
                        let minimum =
                        calculate_viscous_dissipation_from_deformation_gradient_rate_simple!(
                            $constitutive_model_constructed, &get_deformation_gradient_rate()
                        )? - first_piola_kirchoff_stress.full_contraction(
                            &get_deformation_gradient_rate()
                        );
                        let mut perturbed_deformation_gradient_rate = get_deformation_gradient_rate();
                        (0..3).try_for_each(|i|
                            (0..3).try_for_each(|j|{
                                perturbed_deformation_gradient_rate = get_deformation_gradient_rate();
                                perturbed_deformation_gradient_rate[i][j] += 0.5 * EPSILON;
                                assert!(
                                    calculate_viscous_dissipation_from_deformation_gradient_rate_simple!(
                                        $constitutive_model_constructed, &perturbed_deformation_gradient_rate
                                    )? - first_piola_kirchoff_stress.full_contraction(
                                        &perturbed_deformation_gradient_rate
                                    ) > minimum
                                );
                                perturbed_deformation_gradient_rate[i][j] -= EPSILON;
                                assert!(
                                    calculate_viscous_dissipation_from_deformation_gradient_rate_simple!(
                                        $constitutive_model_constructed, &perturbed_deformation_gradient_rate
                                    )? - first_piola_kirchoff_stress.full_contraction(
                                        &perturbed_deformation_gradient_rate
                                    ) > minimum
                                );
                                Ok(())
                            })
                        )
                    }
                    #[test]
                    fn objectivity() -> Result<(), TestError>
                    {
                        assert_eq_within_tols_new(
                            &calculate_viscous_dissipation_from_deformation_gradient_and_deformation_gradient_rate!(
                                $constitutive_model_constructed, &get_deformation_gradient(), &get_deformation_gradient_rate()
                            )?,
                            &calculate_viscous_dissipation_from_deformation_gradient_and_deformation_gradient_rate!(
                                $constitutive_model_constructed, &get_deformation_gradient_rotated(), &get_deformation_gradient_rate_rotated()
                            )?
                        )
                    }
                    #[test]
                    fn positive() -> Result<(), TestError>
                    {
                        assert!(
                            calculate_viscous_dissipation_from_deformation_gradient_rate_simple!(
                                $constitutive_model_constructed,  &get_deformation_gradient_rate()
                            )? > 0.0
                        );
                        Ok(())
                    }
                }
                mod undeformed
                {
                    use super::*;
                    #[test]
                    fn finite_difference() -> Result<(), TestError>
                    {
                        calculate_first_piola_kirchoff_stress_from_finite_difference_of_viscous_dissipation(false)?.iter()
                        .for_each(|fd_first_piola_kirchoff_stress_i|
                            fd_first_piola_kirchoff_stress_i.iter()
                            .for_each(|fd_first_piola_kirchoff_stress_ij|
                                assert!(fd_first_piola_kirchoff_stress_ij.abs() < EPSILON)
                            )
                        );
                        Ok(())
                    }
                    #[test]
                    fn minimized() -> Result<(), TestError>
                    {
                        let minimum =
                        calculate_viscous_dissipation_from_deformation_gradient_rate_simple!(
                            $constitutive_model_constructed, &DeformationGradientRate::zero()
                        )?;
                        let mut perturbed_deformation_gradient_rate = DeformationGradientRate::zero();
                        (0..3).try_for_each(|i|
                            (0..3).try_for_each(|j|{
                                perturbed_deformation_gradient_rate = DeformationGradientRate::zero();
                                perturbed_deformation_gradient_rate[i][j] += 0.5 * EPSILON;
                                assert!(
                                    calculate_viscous_dissipation_from_deformation_gradient_rate_simple!(
                                        $constitutive_model_constructed, &perturbed_deformation_gradient_rate
                                    )? > minimum
                                );
                                perturbed_deformation_gradient_rate[i][j] -= EPSILON;
                                assert!(
                                    calculate_viscous_dissipation_from_deformation_gradient_rate_simple!(
                                        $constitutive_model_constructed, &perturbed_deformation_gradient_rate
                                    )? > minimum
                                );
                                Ok(())
                            })
                        )
                    }
                    #[test]
                    fn zero() -> Result<(), TestError>
                    {
                        assert_eq(
                            &calculate_viscous_dissipation_from_deformation_gradient_rate_simple!(
                                $constitutive_model_constructed,  &DeformationGradientRate::zero()
                            )?, &0.0
                        )
                    }
                }
            }
            mod dissipation_potential
            {
                use super::*;
                fn calculate_first_piola_kirchoff_stress_from_finite_difference_of_dissipation_potential(is_deformed: bool) -> Result<FirstPiolaKirchoffStress, TestError>
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
                            )?;
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
                            )?;
                            first_piola_kirchoff_stress[i][j] = (
                                helmholtz_free_energy_density_plus - helmholtz_free_energy_density_minus
                            )/EPSILON;
                        }
                    }
                    Ok(first_piola_kirchoff_stress)
                }
                mod deformed
                {
                    use super::*;
                    #[test]
                    fn finite_difference() -> Result<(), TestError>
                    {
                        calculate_first_piola_kirchoff_stress_from_deformation_gradient_and_deformation_gradient_rate!(
                            $constitutive_model_constructed, &get_deformation_gradient(), &get_deformation_gradient_rate()
                        )?.iter().zip(
                            calculate_first_piola_kirchoff_stress_from_finite_difference_of_dissipation_potential(true)?.iter()
                        ).for_each(|(first_piola_kirchoff_stress_i, fd_first_piola_kirchoff_stress_i)|
                            first_piola_kirchoff_stress_i.iter()
                            .zip(fd_first_piola_kirchoff_stress_i.iter())
                            .for_each(|(first_piola_kirchoff_stress_ij, fd_first_piola_kirchoff_stress_ij)|
                                assert!((first_piola_kirchoff_stress_ij/fd_first_piola_kirchoff_stress_ij - 1.0).abs() < EPSILON)
                            )
                        );
                        Ok(())
                    }
                    #[test]
                    fn minimized() -> Result<(), TestError>
                    {
                        let first_piola_kirchoff_stress =
                        calculate_first_piola_kirchoff_stress_from_deformation_gradient_and_deformation_gradient_rate!(
                            $constitutive_model_constructed, &get_deformation_gradient(), &get_deformation_gradient_rate()
                        )?;
                        let minimum =
                        calculate_dissipation_potential_from_deformation_gradient_and_deformation_gradient_rate!(
                            $constitutive_model_constructed, &get_deformation_gradient(), &get_deformation_gradient_rate()
                        )? - first_piola_kirchoff_stress.full_contraction(
                            &get_deformation_gradient_rate()
                        );
                        let mut perturbed_deformation_gradient_rate = get_deformation_gradient_rate();
                        (0..3).try_for_each(|i|
                            (0..3).try_for_each(|j|{
                                perturbed_deformation_gradient_rate = get_deformation_gradient_rate();
                                perturbed_deformation_gradient_rate[i][j] += 0.5 * EPSILON;
                                assert!(
                                    calculate_dissipation_potential_from_deformation_gradient_and_deformation_gradient_rate!(
                                        $constitutive_model_constructed, &get_deformation_gradient(), &perturbed_deformation_gradient_rate
                                    )? - first_piola_kirchoff_stress.full_contraction(
                                        &perturbed_deformation_gradient_rate
                                    ) > minimum
                                );
                                perturbed_deformation_gradient_rate[i][j] -= EPSILON;
                                assert!(
                                    calculate_dissipation_potential_from_deformation_gradient_and_deformation_gradient_rate!(
                                        $constitutive_model_constructed, &get_deformation_gradient(), &perturbed_deformation_gradient_rate
                                    )? - first_piola_kirchoff_stress.full_contraction(
                                        &perturbed_deformation_gradient_rate
                                    ) > minimum
                                );
                                Ok(())
                            })
                        )
                    }
                    #[test]
                    fn objectivity() -> Result<(), TestError>
                    {
                        assert_eq_within_tols_new(
                            &calculate_dissipation_potential_from_deformation_gradient_and_deformation_gradient_rate!(
                                $constitutive_model_constructed, &get_deformation_gradient(), &get_deformation_gradient_rate()
                            )?,
                            &calculate_dissipation_potential_from_deformation_gradient_and_deformation_gradient_rate!(
                                $constitutive_model_constructed, &get_deformation_gradient_rotated(), &get_deformation_gradient_rate_rotated()
                            )?
                        )
                    }
                }
                mod undeformed
                {
                    use super::*;
                    #[test]
                    fn finite_difference() -> Result<(), TestError>
                    {
                        calculate_first_piola_kirchoff_stress_from_finite_difference_of_dissipation_potential(false)?.iter()
                        .for_each(|fd_first_piola_kirchoff_stress_i|
                            fd_first_piola_kirchoff_stress_i.iter()
                            .for_each(|fd_first_piola_kirchoff_stress_ij|
                                assert!(fd_first_piola_kirchoff_stress_ij.abs() < EPSILON)
                            )
                        );
                        Ok(())
                    }
                    #[test]
                    fn minimized() -> Result<(), TestError>
                    {
                        let minimum =
                        calculate_dissipation_potential_from_deformation_gradient_and_deformation_gradient_rate!(
                            $constitutive_model_constructed, &DeformationGradient::identity(), &DeformationGradientRate::zero()
                        )?;
                        let mut perturbed_deformation_gradient_rate = DeformationGradientRate::zero();
                        (0..3).try_for_each(|i|
                            (0..3).try_for_each(|j|{
                                perturbed_deformation_gradient_rate = DeformationGradientRate::zero();
                                perturbed_deformation_gradient_rate[i][j] += 0.5 * EPSILON;
                                assert!(
                                    calculate_dissipation_potential_from_deformation_gradient_and_deformation_gradient_rate!(
                                        $constitutive_model_constructed, &DeformationGradient::identity(), &perturbed_deformation_gradient_rate
                                    )? > minimum
                                );
                                perturbed_deformation_gradient_rate[i][j] -= EPSILON;
                                assert!(
                                    calculate_dissipation_potential_from_deformation_gradient_and_deformation_gradient_rate!(
                                        $constitutive_model_constructed, &DeformationGradient::identity(), &perturbed_deformation_gradient_rate
                                    )? > minimum
                                );
                                Ok(())
                            })
                        )
                    }
                    #[test]
                    fn zero() -> Result<(), TestError>
                    {
                        assert_eq(
                            &calculate_dissipation_potential_from_deformation_gradient_and_deformation_gradient_rate!(
                                $constitutive_model_constructed, &DeformationGradient::identity(), &DeformationGradientRate::zero()
                            )?, &0.0
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
                    fn symmetry() -> Result<(), TestError>
                    {
                        let first_piola_kirchoff_rate_tangent_stiffness =
                        calculate_first_piola_kirchoff_rate_tangent_stiffness_from_deformation_gradient_and_deformation_gradient_rate!(
                            $constitutive_model_constructed, &get_deformation_gradient(), &get_deformation_gradient_rate()
                        )?;
                        let result = first_piola_kirchoff_rate_tangent_stiffness.iter().enumerate()
                        .try_for_each(|(i, first_piola_kirchoff_rate_tangent_stiffness_i)|
                            first_piola_kirchoff_rate_tangent_stiffness_i.iter().enumerate()
                            .try_for_each(|(j, first_piola_kirchoff_rate_tangent_stiffness_ij)|
                                first_piola_kirchoff_rate_tangent_stiffness_ij.iter()
                                .zip(first_piola_kirchoff_rate_tangent_stiffness.iter())
                                .try_for_each(|(first_piola_kirchoff_rate_tangent_stiffness_ijk, first_piola_kirchoff_rate_tangent_stiffness_k)|
                                    first_piola_kirchoff_rate_tangent_stiffness_ijk.iter()
                                    .zip(first_piola_kirchoff_rate_tangent_stiffness_k.iter())
                                    .try_for_each(|(first_piola_kirchoff_rate_tangent_stiffness_ijkl, first_piola_kirchoff_rate_tangent_stiffness_kl)|
                                        assert_eq_within_tols_new(
                                            first_piola_kirchoff_rate_tangent_stiffness_ijkl, &first_piola_kirchoff_rate_tangent_stiffness_kl[i][j]
                                        )
                                    )
                                )
                            )
                        );
                        result
                    }
                }
                mod undeformed
                {
                    use super::*;
                    #[test]
                    fn symmetry() -> Result<(), TestError>
                    {
                        let first_piola_kirchoff_rate_tangent_stiffness =
                        calculate_first_piola_kirchoff_rate_tangent_stiffness_from_deformation_gradient_and_deformation_gradient_rate!(
                            $constitutive_model_constructed, &DeformationGradient::identity(), &DeformationGradientRate::zero()
                        )?;
                        let result = first_piola_kirchoff_rate_tangent_stiffness.iter().enumerate()
                        .try_for_each(|(i, first_piola_kirchoff_rate_tangent_stiffness_i)|
                            first_piola_kirchoff_rate_tangent_stiffness_i.iter().enumerate()
                            .try_for_each(|(j, first_piola_kirchoff_rate_tangent_stiffness_ij)|
                                first_piola_kirchoff_rate_tangent_stiffness_ij.iter()
                                .zip(first_piola_kirchoff_rate_tangent_stiffness.iter())
                                .try_for_each(|(first_piola_kirchoff_rate_tangent_stiffness_ijk, first_piola_kirchoff_rate_tangent_stiffness_k)|
                                    first_piola_kirchoff_rate_tangent_stiffness_ijk.iter()
                                    .zip(first_piola_kirchoff_rate_tangent_stiffness_k.iter())
                                    .try_for_each(|(first_piola_kirchoff_rate_tangent_stiffness_ijkl, first_piola_kirchoff_rate_tangent_stiffness_kl)|
                                        assert_eq_within_tols_new(
                                            first_piola_kirchoff_rate_tangent_stiffness_ijkl, &first_piola_kirchoff_rate_tangent_stiffness_kl[i][j]
                                        )
                                    )
                                )
                            )
                        );
                        result
                    }
                }
            }
        }
    }
}
pub(crate) use test_solid_elastic_hyperviscous_specifics;
