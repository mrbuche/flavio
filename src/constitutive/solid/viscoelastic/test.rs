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

macro_rules! calculate_first_piola_kirchoff_stress_from_deformation_gradient_rate
{
    ($constitutive_model_constructed: expr, $deformation_gradient_rate: expr) =>
    {
        $constitutive_model_constructed.calculate_first_piola_kirchoff_stress(
            &DeformationGradient::identity(), $deformation_gradient_rate
        )
    }
}
pub(crate) use calculate_first_piola_kirchoff_stress_from_deformation_gradient_rate;

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
            let model = get_constitutive_model();
            let mut deformation_gradient_rate = DeformationGradientRate::zero();
            deformation_gradient_rate += DeformationGradientRate::identity()*(EPSILON/3.0);
            let first_piola_kirchoff_stress = calculate_first_piola_kirchoff_stress_from_deformation_gradient_rate!(&model, &deformation_gradient_rate);
            assert!((3.0*EPSILON*model.get_bulk_modulus()/first_piola_kirchoff_stress.trace() - 1.0).abs() < 3.0*EPSILON);
        }
        #[test]
        fn shear_viscosity()
        {
            let model = get_constitutive_model();
            let mut deformation_gradient_rate = DeformationGradientRate::zero();
            deformation_gradient_rate[0][1] = EPSILON;
            let first_piola_kirchoff_stress = calculate_first_piola_kirchoff_stress_from_deformation_gradient_rate!(&model, &deformation_gradient_rate);
            assert!((EPSILON*model.get_shear_viscosity()/first_piola_kirchoff_stress[0][1] - 1.0).abs() < EPSILON)
        }
        #[test]
        fn repeat_all_elastic_tests_as_relevant_here()
        {
            todo!("need tests for tangent/rate-tangent where the non-important ESV is non-trivial
                   because dF of those strain rates in stress are nonzero")
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