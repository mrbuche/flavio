use crate::
{
    constitutive::solid::elastic::test::ALMANSIHAMELPARAMETERS,
    mechanics::Scalar
};

pub const ARRUDABOYCEPARAMETERS: &[Scalar; 3] = &[ALMANSIHAMELPARAMETERS[0], ALMANSIHAMELPARAMETERS[1], 8.0];
pub const FUNGPARAMETERS: &[Scalar; 4] = &[ALMANSIHAMELPARAMETERS[0], ALMANSIHAMELPARAMETERS[1], 1.0, 1.0];
pub const GENTPARAMETERS: &[Scalar; 3] = &[ALMANSIHAMELPARAMETERS[0], ALMANSIHAMELPARAMETERS[1], 23.0];
pub const MOONEYRIVLINPARAMETERS: &[Scalar; 3] = &[ALMANSIHAMELPARAMETERS[0], ALMANSIHAMELPARAMETERS[1], 1.0];
pub const NEOHOOKEANPARAMETERS: &[Scalar; 2] = &[ALMANSIHAMELPARAMETERS[0], ALMANSIHAMELPARAMETERS[1]];
pub const SAINTVENANTKIRCHOFFPARAMETERS: &[Scalar; 2] = &[ALMANSIHAMELPARAMETERS[0], ALMANSIHAMELPARAMETERS[1]];
pub const YEOHPARAMETERS: &[Scalar; 6] = &[ALMANSIHAMELPARAMETERS[0], ALMANSIHAMELPARAMETERS[1], -1.0, 3e-1, -1e-3, 1e-5];

macro_rules! calculate_helmholtz_free_energy_density_from_deformation_gradient_simple
{
    ($constitutive_model_constructed: expr, $deformation_gradient: expr) =>
    {
        $constitutive_model_constructed.calculate_helmholtz_free_energy_density(
            $deformation_gradient
        )
    }
}
pub(crate) use calculate_helmholtz_free_energy_density_from_deformation_gradient_simple;

macro_rules! use_elastic_macros
{
    () =>
    {
        use crate::constitutive::solid::elastic::test::
        {
            calculate_cauchy_stress_from_deformation_gradient,
            calculate_cauchy_stress_from_deformation_gradient_simple,
            calculate_cauchy_stress_from_deformation_gradient_rotated,
            calculate_cauchy_tangent_stiffness_from_deformation_gradient,
            calculate_cauchy_tangent_stiffness_from_deformation_gradient_rotated,
            calculate_first_piola_kirchoff_stress_from_deformation_gradient,
            calculate_first_piola_kirchoff_stress_from_deformation_gradient_simple,
            calculate_first_piola_kirchoff_stress_from_deformation_gradient_rotated,
            calculate_first_piola_kirchoff_tangent_stiffness_from_deformation_gradient,
            calculate_first_piola_kirchoff_tangent_stiffness_from_deformation_gradient_simple,
            calculate_second_piola_kirchoff_stress_from_deformation_gradient,
            calculate_second_piola_kirchoff_stress_from_deformation_gradient_simple,
            calculate_second_piola_kirchoff_stress_from_deformation_gradient_rotated,
            calculate_second_piola_kirchoff_tangent_stiffness_from_deformation_gradient
        };
    }
}
pub(crate) use use_elastic_macros;

macro_rules! test_solid_hyperelastic_constitutive_model
{
    ($constitutive_model: ident, $constitutive_model_parameters: expr, $constitutive_model_constructed: expr) =>
    {
        crate::constitutive::solid::elastic::test::test_solid_constitutive_model!(
            $constitutive_model,
            $constitutive_model_parameters,
            $constitutive_model_constructed
        );
        mod hyperelastic
        {
            use super::*;
            fn calculate_first_piola_kirchoff_stress_from_finite_difference_of_helmholtz_free_energy_density(is_deformed: bool) -> FirstPiolaKirchoffStress
            {
                let mut first_piola_kirchoff_stress = FirstPiolaKirchoffStress::zero();
                for i in 0..3
                {
                    for j in 0..3
                    {
                        let mut deformation_gradient_plus = 
                            if is_deformed
                            {
                                get_deformation_gradient()
                            }
                            else
                            {
                                DeformationGradient::identity()
                            };
                        deformation_gradient_plus[i][j] += 0.5*EPSILON;
                        let helmholtz_free_energy_density_plus =
                        calculate_helmholtz_free_energy_density_from_deformation_gradient_simple!(
                            $constitutive_model_constructed, &deformation_gradient_plus
                        );
                        let mut deformation_gradient_minus = 
                            if is_deformed
                            {
                                get_deformation_gradient()
                            }
                            else
                            {
                                DeformationGradient::identity()
                            };
                        deformation_gradient_minus[i][j] -= 0.5*EPSILON;
                        let helmholtz_free_energy_density_minus =
                        calculate_helmholtz_free_energy_density_from_deformation_gradient_simple!(
                            $constitutive_model_constructed, &deformation_gradient_minus
                        );
                        first_piola_kirchoff_stress[i][j] = (
                            helmholtz_free_energy_density_plus - helmholtz_free_energy_density_minus
                        )/EPSILON;
                    }
                }
                first_piola_kirchoff_stress
            }
            mod helmholtz_free_energy_density
            {
                use super::*;
                mod deformed
                {
                    use super::*;
                    #[test]
                    fn finite_difference()
                    {
                        calculate_first_piola_kirchoff_stress_from_deformation_gradient_simple!(
                            $constitutive_model_constructed, &get_deformation_gradient()
                        ).iter().zip(
                            calculate_first_piola_kirchoff_stress_from_finite_difference_of_helmholtz_free_energy_density(true).iter()
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
                        calculate_first_piola_kirchoff_stress_from_deformation_gradient_simple!(
                            $constitutive_model_constructed, &get_deformation_gradient()
                        );
                        let minimum =
                        calculate_helmholtz_free_energy_density_from_deformation_gradient_simple!(
                            $constitutive_model_constructed, &get_deformation_gradient()
                        ) - first_piola_kirchoff_stress.full_contraction(
                            &get_deformation_gradient()
                        );
                        let mut perturbed_deformation_gradient = get_deformation_gradient();
                        (0..3).for_each(|i|
                            (0..3).for_each(|j|{
                                perturbed_deformation_gradient = get_deformation_gradient();
                                perturbed_deformation_gradient[i][j] += 0.5 * EPSILON;
                                assert!(
                                    calculate_helmholtz_free_energy_density_from_deformation_gradient_simple!(
                                        $constitutive_model_constructed, &perturbed_deformation_gradient
                                    ) - first_piola_kirchoff_stress.full_contraction(
                                        &perturbed_deformation_gradient
                                    ) > minimum
                                );
                                perturbed_deformation_gradient[i][j] -= EPSILON;
                                assert!(
                                    calculate_helmholtz_free_energy_density_from_deformation_gradient_simple!(
                                        $constitutive_model_constructed, &perturbed_deformation_gradient
                                    ) - first_piola_kirchoff_stress.full_contraction(
                                        &perturbed_deformation_gradient
                                    ) > minimum
                                );
                            })
                        )
                    }
                    #[test]
                    fn objectivity()
                    {
                        assert_eq_within_tols(
                            &calculate_helmholtz_free_energy_density_from_deformation_gradient_simple!(
                                $constitutive_model_constructed,  &get_deformation_gradient()
                            ),
                            &calculate_helmholtz_free_energy_density_from_deformation_gradient_simple!(
                                $constitutive_model_constructed,  &get_deformation_gradient_rotated()
                            )
                        )
                    }
                    #[test]
                    fn positive()
                    {
                        assert!(
                            calculate_helmholtz_free_energy_density_from_deformation_gradient_simple!(
                                $constitutive_model_constructed,  &get_deformation_gradient()
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
                        calculate_first_piola_kirchoff_stress_from_finite_difference_of_helmholtz_free_energy_density(false).iter()
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
                        calculate_helmholtz_free_energy_density_from_deformation_gradient_simple!(
                            $constitutive_model_constructed, &DeformationGradient::identity()
                        );
                        let mut perturbed_deformation_gradient = DeformationGradient::identity();
                        (0..3).for_each(|i|
                            (0..3).for_each(|j|{
                                perturbed_deformation_gradient = DeformationGradient::identity();
                                perturbed_deformation_gradient[i][j] += 0.5 * EPSILON;
                                assert!(
                                    calculate_helmholtz_free_energy_density_from_deformation_gradient_simple!(
                                        $constitutive_model_constructed, &perturbed_deformation_gradient
                                    ) > minimum
                                );
                                perturbed_deformation_gradient[i][j] -= EPSILON;
                                assert!(
                                    calculate_helmholtz_free_energy_density_from_deformation_gradient_simple!(
                                        $constitutive_model_constructed, &perturbed_deformation_gradient
                                    ) > minimum
                                );
                            })
                        )
                    }
                    #[test]
                    fn zero()
                    {
                        assert_eq!(
                            calculate_helmholtz_free_energy_density_from_deformation_gradient_simple!(
                                $constitutive_model_constructed,  &DeformationGradient::identity()
                            ), 0.0
                        )
                    }
                }
            }
            mod first_piola_kirchoff_tangent_stiffness
            {
                use super::*;
                mod deformed
                {
                    use super::*;
                    #[test]
                    fn symmetry()
                    {
                        let first_piola_kirchoff_tangent_stiffness =
                        calculate_first_piola_kirchoff_tangent_stiffness_from_deformation_gradient_simple!(
                            $constitutive_model_constructed, &get_deformation_gradient()
                        );
                        first_piola_kirchoff_tangent_stiffness.iter().enumerate()
                        .for_each(|(i, first_piola_kirchoff_tangent_stiffness_i)|
                            first_piola_kirchoff_tangent_stiffness_i.iter().enumerate()
                            .for_each(|(j, first_piola_kirchoff_tangent_stiffness_ij)|
                                first_piola_kirchoff_tangent_stiffness_ij.iter()
                                .zip(first_piola_kirchoff_tangent_stiffness.iter())
                                .for_each(|(first_piola_kirchoff_tangent_stiffness_ijk, first_piola_kirchoff_tangent_stiffness_k)|
                                    first_piola_kirchoff_tangent_stiffness_ijk.iter()
                                    .zip(first_piola_kirchoff_tangent_stiffness_k.iter())
                                    .for_each(|(first_piola_kirchoff_tangent_stiffness_ijkl, first_piola_kirchoff_tangent_stiffness_kl)|
                                        assert_eq_within_tols(
                                            first_piola_kirchoff_tangent_stiffness_ijkl, &first_piola_kirchoff_tangent_stiffness_kl[i][j]
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
                        let first_piola_kirchoff_tangent_stiffness =
                        calculate_first_piola_kirchoff_tangent_stiffness_from_deformation_gradient_simple!(
                            $constitutive_model_constructed, &DeformationGradient::identity()
                        );
                        first_piola_kirchoff_tangent_stiffness.iter().enumerate()
                        .for_each(|(i, first_piola_kirchoff_tangent_stiffness_i)|
                            first_piola_kirchoff_tangent_stiffness_i.iter().enumerate()
                            .for_each(|(j, first_piola_kirchoff_tangent_stiffness_ij)|
                                first_piola_kirchoff_tangent_stiffness_ij.iter()
                                .zip(first_piola_kirchoff_tangent_stiffness.iter())
                                .for_each(|(first_piola_kirchoff_tangent_stiffness_ijk, first_piola_kirchoff_tangent_stiffness_k)|
                                    first_piola_kirchoff_tangent_stiffness_ijk.iter()
                                    .zip(first_piola_kirchoff_tangent_stiffness_k.iter())
                                    .for_each(|(first_piola_kirchoff_tangent_stiffness_ijkl, first_piola_kirchoff_tangent_stiffness_kl)|
                                        assert_eq_within_tols(
                                            first_piola_kirchoff_tangent_stiffness_ijkl, &first_piola_kirchoff_tangent_stiffness_kl[i][j]
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
pub(crate) use test_solid_hyperelastic_constitutive_model;