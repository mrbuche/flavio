use crate::
{
    constitutive::solid::viscoelastic::test::ALMANSIHAMELPARAMETERS,
    mechanics::Scalar
};
pub const SAINTVENANTKIRCHOFFPARAMETERS: &[Scalar; 4] = &[ALMANSIHAMELPARAMETERS[0], ALMANSIHAMELPARAMETERS[1], ALMANSIHAMELPARAMETERS[2], ALMANSIHAMELPARAMETERS[3]];

macro_rules! calculate_helmholtz_free_energy_density_from_deformation_gradient
{
    ($constitutive_model_constructed: expr, $deformation_gradient: expr) =>
    {
        $constitutive_model_constructed.calculate_helmholtz_free_energy_density(
            $deformation_gradient
        )
    }
}
pub(crate) use calculate_helmholtz_free_energy_density_from_deformation_gradient;

macro_rules! use_viscoelastic_macros
{
    () =>
    {
        use crate::constitutive::solid::viscoelastic::test::
        {
            calculate_cauchy_stress_from_deformation_gradient,
            calculate_cauchy_tangent_stiffness_from_deformation_gradient,
            calculate_first_piola_kirchoff_stress_from_deformation_gradient,
            calculate_first_piola_kirchoff_tangent_stiffness_from_deformation_gradient,
            calculate_second_piola_kirchoff_stress_from_deformation_gradient,
            calculate_second_piola_kirchoff_tangent_stiffness_from_deformation_gradient
        };
    }
}
pub(crate) use use_viscoelastic_macros;

macro_rules! test_solid_hyperviscoelastic_constitutive_model
{
    ($constitutive_model: ident, $constitutive_model_parameters: expr, $constitutive_model_constructed: expr) =>
    {
        crate::constitutive::solid::hyperelastic::test::test_solid_hyperelastic_constitutive_model!(
            $constitutive_model,
            $constitutive_model_parameters,
            $constitutive_model_constructed
        );
        crate::constitutive::solid::viscoelastic::test::test_solid_viscous_constitutive_model!(
            $constitutive_model,
            $constitutive_model_parameters,
            $constitutive_model_constructed
        );
        mod hyperviscoelastic
        {
            #[test]
            fn fd_for_first_pk()
            {
                todo!()
            }
            #[test]
            fn first_pk_range_tangent_symmetry()
            {
                todo!()
            }
            mod helmholtz_free_energy_density
            {
                mod deformed
                {
                    #[test]
                    fn finite_difference()
                    {
                        todo!()
                    }
                    #[test]
                    fn minimized()
                    {
                        todo!()
                    }
                    #[test]
                    fn objectivity()
                    {
                        todo!()
                    }
                    #[test]
                    fn positive()
                    {
                        todo!()
                    }
                }
                mod undeformed
                {
                    #[test]
                    fn finite_difference()
                    {
                        todo!()
                    }
                    #[test]
                    fn minimized()
                    {
                        todo!()
                    }
                    #[test]
                    fn zero()
                    {
                        todo!()
                    }
                }
            }
            mod viscous_dissipation
            {
                mod deformed
                {
                    #[test]
                    fn finite_difference()
                    {
                        todo!()
                    }
                    #[test]
                    fn minimized()
                    {
                        todo!()
                    }
                    #[test]
                    fn objectivity()
                    {
                        todo!()
                    }
                    #[test]
                    fn positive()
                    {
                        todo!()
                    }
                }
                mod undeformed
                {
                    #[test]
                    fn finite_difference()
                    {
                        todo!()
                    }
                    #[test]
                    fn minimized()
                    {
                        todo!()
                    }
                    #[test]
                    fn zero()
                    {
                        todo!()
                    }
                }
            }
            mod helmholtz_free_energy_density_plus_viscous_dissipation
            {
                mod deformed
                {
                    #[test]
                    fn finite_difference()
                    {
                        todo!()
                    }
                    #[test]
                    fn minimized()
                    {
                        todo!()
                    }
                    #[test]
                    fn objectivity()
                    {
                        todo!()
                    }
                    #[test]
                    fn positive()
                    {
                        todo!()
                    }
                }
                mod undeformed
                {
                    #[test]
                    fn finite_difference()
                    {
                        todo!()
                    }
                    #[test]
                    fn minimized()
                    {
                        todo!()
                    }
                    #[test]
                    fn zero()
                    {
                        todo!()
                    }
                }
            }
            mod first_piola_kirchoff_rate_tangent_stiffness
            {
                mod deformed
                {
                    #[test]
                    fn symmetry()
                    {
                        todo!()
                    }
                }
                mod undeformed
                {
                    #[test]
                    fn symmetry()
                    {
                        todo!()
                    }
                }
            }
        }
    }
}
pub(crate) use test_solid_hyperviscoelastic_constitutive_model;