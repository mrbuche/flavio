use crate::
{
    constitutive::solid::elastic::test::ALMANSIHAMELPARAMETERS as ALMANSIHAMELPARAMETERSELASTIC,
    mechanics::Scalar
};
pub const ALMANSIHAMELPARAMETERS: &[Scalar; 4] = &[ALMANSIHAMELPARAMETERSELASTIC[0], ALMANSIHAMELPARAMETERSELASTIC[1], ALMANSIHAMELPARAMETERSELASTIC[0], ALMANSIHAMELPARAMETERSELASTIC[1]];

macro_rules! calculate_cauchy_stress_from_deformation_gradient
{
    ($constitutive_model_constructed: expr, $deformation_gradient: expr) =>
    {
        $constitutive_model_constructed.calculate_cauchy_stress(
            $deformation_gradient, &DeformationGradientRate::zero()
        )
    }
}
pub(crate) use calculate_cauchy_stress_from_deformation_gradient;

macro_rules! calculate_cauchy_tangent_stiffness_from_deformation_gradient
{
    ($constitutive_model_constructed: expr, $deformation_gradient: expr) =>
    {
        $constitutive_model_constructed.calculate_cauchy_tangent_stiffness(
            $deformation_gradient, &DeformationGradientRate::zero()
        )
    }
}
pub(crate) use calculate_cauchy_tangent_stiffness_from_deformation_gradient;

macro_rules! calculate_first_piola_kirchoff_stress_from_deformation_gradient
{
    ($constitutive_model_constructed: expr, $deformation_gradient: expr) =>
    {
        $constitutive_model_constructed.calculate_first_piola_kirchoff_stress(
            $deformation_gradient, &DeformationGradientRate::zero()
        )
    }
}
pub(crate) use calculate_first_piola_kirchoff_stress_from_deformation_gradient;

macro_rules! calculate_first_piola_kirchoff_tangent_stiffness_from_deformation_gradient
{
    ($constitutive_model_constructed: expr, $deformation_gradient: expr) =>
    {
        $constitutive_model_constructed.calculate_first_piola_kirchoff_tangent_stiffness(
            $deformation_gradient, &DeformationGradientRate::zero()
        )
    }
}
pub(crate) use calculate_first_piola_kirchoff_tangent_stiffness_from_deformation_gradient;

macro_rules! calculate_second_piola_kirchoff_stress_from_deformation_gradient
{
    ($constitutive_model_constructed: expr, $deformation_gradient: expr) =>
    {
        $constitutive_model_constructed.calculate_second_piola_kirchoff_stress(
            $deformation_gradient, &DeformationGradientRate::zero()
        )
    }
}
pub(crate) use calculate_second_piola_kirchoff_stress_from_deformation_gradient;

macro_rules! calculate_second_piola_kirchoff_tangent_stiffness_from_deformation_gradient
{
    ($constitutive_model_constructed: expr, $deformation_gradient: expr) =>
    {
        $constitutive_model_constructed.calculate_second_piola_kirchoff_tangent_stiffness(
            $deformation_gradient, &DeformationGradientRate::zero()
        )
    }
}
pub(crate) use calculate_second_piola_kirchoff_tangent_stiffness_from_deformation_gradient;

macro_rules! test_solid_viscous_constitutive_model
{
    ($constitutive_model: ident, $constitutive_model_parameters: expr, $constitutive_model_constructed: expr) =>
    {
        #[test]
        fn get_bulk_viscosity()
        {
            assert_eq!(get_constitutive_model().get_bulk_viscosity(), &$constitutive_model_parameters[2])
        }
        #[test]
        fn get_shear_viscosity()
        {
            assert_eq!(get_constitutive_model().get_shear_viscosity(), &$constitutive_model_parameters[3])
        }
        #[test]
        fn bulk_viscosity()
        {
            todo!()
        }
        #[test]
        fn shear_viscosity()
        {
            todo!()
        }
        #[test]
        fn repeat_all_elastic_tests_as_relevant_here()
        {
            todo!()
        }
    }
}
pub(crate) use test_solid_viscous_constitutive_model;

macro_rules! test_solid_viscoelastic_constitutive_model
{
    ($constitutive_model: ident, $constitutive_model_parameters: expr, $constitutive_model_constructed: expr) =>
    {
        crate::constitutive::solid::elastic::test::test_solid_elastic_constitutive_model!(
            $constitutive_model,
            $constitutive_model_parameters,
            $constitutive_model_constructed
        );
        crate::constitutive::solid::viscoelastic::test::test_solid_viscous_constitutive_model!(
            $constitutive_model,
            $constitutive_model_parameters,
            $constitutive_model_constructed
        );
    }
}
pub(crate) use test_solid_viscoelastic_constitutive_model;