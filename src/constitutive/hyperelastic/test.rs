macro_rules! test_hyperelastic_constitutive_model
{
    ($hyperelastic_constitutive_model: ident, $hyperelastic_constitutive_model_parameters: expr, $hyperelastic_constitutive_model_constructed: expr) =>
    {
        use crate::constitutive::
        {
            elastic::test::test_elastic_constitutive_model,
            hyperelastic::
            {
                HyperelasticConstitutiveModel,
                test::test_hyperelastic_constitutive_model_constructed
            }
        };
        test_elastic_constitutive_model!($hyperelastic_constitutive_model, $hyperelastic_constitutive_model_parameters, $hyperelastic_constitutive_model_constructed);
        test_hyperelastic_constitutive_model_constructed!($hyperelastic_constitutive_model_constructed);
    }
}
pub(crate) use test_hyperelastic_constitutive_model;
macro_rules! test_hyperelastic_constitutive_model_constructed
{
    ($hyperelastic_constitutive_model_constructed: expr) =>
    {
        mod hyperelastic_only
        {
            use crate::constitutive::FirstPiolaKirchoffStress;
            use super::*;
            fn calculate_first_piola_kirchoff_stress_from_finite_difference_of_helmholtz_free_energy_density(is_deformed: bool) -> FirstPiolaKirchoffStress
            {
                let model = $hyperelastic_constitutive_model_constructed;
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
                        let helmholtz_free_energy_density_plus = model.calculate_helmholtz_free_energy_density(&deformation_gradient_plus);
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
                        let helmholtz_free_energy_density_minus = model.calculate_helmholtz_free_energy_density(&deformation_gradient_minus);
                        first_piola_kirchoff_stress[i][j] = (helmholtz_free_energy_density_plus - helmholtz_free_energy_density_minus)/EPSILON;
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
                        $hyperelastic_constitutive_model_constructed.calculate_first_piola_kirchoff_stress(&get_deformation_gradient()).iter()
                        .zip(calculate_first_piola_kirchoff_stress_from_finite_difference_of_helmholtz_free_energy_density(true).iter())
                        .for_each(|(first_piola_kirchoff_stress_i, fd_first_piola_kirchoff_stress_i)|
                            first_piola_kirchoff_stress_i.iter()
                            .zip(fd_first_piola_kirchoff_stress_i.iter())
                            .for_each(|(first_piola_kirchoff_stress_ij, fd_first_piola_kirchoff_stress_ij)|
                                assert!((first_piola_kirchoff_stress_ij/fd_first_piola_kirchoff_stress_ij - 1.0).abs() < EPSILON)
                            )
                        )
                    }
                    #[test]
                    fn objectivity()
                    {
                        let model = $hyperelastic_constitutive_model_constructed;
                        assert_eq_within_tols(
                            &model.calculate_helmholtz_free_energy_density(&get_deformation_gradient()),
                            &model.calculate_helmholtz_free_energy_density(&get_deformation_gradient_rotated())
                        )
                    }
                    #[test]
                    fn positive()
                    {
                        assert!($hyperelastic_constitutive_model_constructed.calculate_helmholtz_free_energy_density(&get_deformation_gradient()) > 0.0)
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
                    fn zero()
                    {
                        assert_eq!($hyperelastic_constitutive_model_constructed.calculate_helmholtz_free_energy_density(&DeformationGradient::identity()), 0.0)
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
                        let model = $hyperelastic_constitutive_model_constructed;
                        let first_piola_kirchoff_tangent_stiffness = model.calculate_first_piola_kirchoff_tangent_stiffness(&get_deformation_gradient());
                        for i in 0..3
                        {
                            for j in 0..3
                            {
                                for k in 0..3
                                {
                                    for l in 0..3
                                    {
                                        assert_eq_within_tols(
                                            &first_piola_kirchoff_tangent_stiffness[i][j][k][l],
                                            &first_piola_kirchoff_tangent_stiffness[k][l][i][j]
                                        )
                                    }
                                }
                            }
                        }
                    }
                }
                mod undeformed
                {
                    use super::*;
                    #[test]
                    fn symmetry()
                    {
                        let model = $hyperelastic_constitutive_model_constructed;
                        let first_piola_kirchoff_tangent_stiffness = model.calculate_first_piola_kirchoff_tangent_stiffness(&DeformationGradient::identity());
                        for i in 0..3
                        {
                            for j in 0..3
                            {
                                for k in 0..3
                                {
                                    for l in 0..3
                                    {
                                        assert_eq_within_tols(
                                            &first_piola_kirchoff_tangent_stiffness[i][j][k][l],
                                            &first_piola_kirchoff_tangent_stiffness[k][l][i][j]
                                        )
                                    }
                                }
                            }
                        }
                    }
                }

            }
        }
    }
}
pub(crate) use test_hyperelastic_constitutive_model_constructed;