use crate::mechanics::Scalar;

pub const ALMANSIHAMELPARAMETERS: &[Scalar; 2] = &[13.0, 3.0];

#[test]
fn finish_getting_rid_of_loops_and_old_assert_eq_with_tols() {
    todo!()
}

#[test]
fn implement_and_use_something_for_epsilon_fd_comparisons() {
    todo!()
}

macro_rules! calculate_cauchy_stress_from_deformation_gradient {
    ($constitutive_model_constructed: expr, $deformation_gradient: expr) => {
        $constitutive_model_constructed.calculate_cauchy_stress($deformation_gradient)
    };
}
pub(crate) use calculate_cauchy_stress_from_deformation_gradient;

macro_rules! calculate_cauchy_stress_from_deformation_gradient_simple {
    ($constitutive_model_constructed: expr, $deformation_gradient: expr) => {
        $constitutive_model_constructed.calculate_cauchy_stress($deformation_gradient)
    };
}
pub(crate) use calculate_cauchy_stress_from_deformation_gradient_simple;

macro_rules! calculate_cauchy_stress_from_deformation_gradient_rotated {
    ($constitutive_model_constructed: expr, $deformation_gradient: expr) => {
        $constitutive_model_constructed.calculate_cauchy_stress($deformation_gradient)
    };
}
pub(crate) use calculate_cauchy_stress_from_deformation_gradient_rotated;

macro_rules! calculate_cauchy_tangent_stiffness_from_deformation_gradient {
    ($constitutive_model_constructed: expr, $deformation_gradient: expr) => {
        $constitutive_model_constructed.calculate_cauchy_tangent_stiffness($deformation_gradient)
    };
}
pub(crate) use calculate_cauchy_tangent_stiffness_from_deformation_gradient;

macro_rules! calculate_first_piola_kirchoff_stress_from_deformation_gradient {
    ($constitutive_model_constructed: expr, $deformation_gradient: expr) => {
        $constitutive_model_constructed.calculate_first_piola_kirchoff_stress($deformation_gradient)
    };
}
pub(crate) use calculate_first_piola_kirchoff_stress_from_deformation_gradient;

macro_rules! calculate_first_piola_kirchoff_stress_from_deformation_gradient_simple {
    ($constitutive_model_constructed: expr, $deformation_gradient: expr) => {
        $constitutive_model_constructed.calculate_first_piola_kirchoff_stress($deformation_gradient)
    };
}
pub(crate) use calculate_first_piola_kirchoff_stress_from_deformation_gradient_simple;

macro_rules! calculate_first_piola_kirchoff_stress_from_deformation_gradient_rotated {
    ($constitutive_model_constructed: expr, $deformation_gradient: expr) => {
        $constitutive_model_constructed.calculate_first_piola_kirchoff_stress($deformation_gradient)
    };
}
pub(crate) use calculate_first_piola_kirchoff_stress_from_deformation_gradient_rotated;

macro_rules! calculate_first_piola_kirchoff_tangent_stiffness_from_deformation_gradient {
    ($constitutive_model_constructed: expr, $deformation_gradient: expr) => {
        $constitutive_model_constructed
            .calculate_first_piola_kirchoff_tangent_stiffness($deformation_gradient)
    };
}
pub(crate) use calculate_first_piola_kirchoff_tangent_stiffness_from_deformation_gradient;

macro_rules! calculate_first_piola_kirchoff_tangent_stiffness_from_deformation_gradient_simple {
    ($constitutive_model_constructed: expr, $deformation_gradient: expr) => {
        $constitutive_model_constructed
            .calculate_first_piola_kirchoff_tangent_stiffness($deformation_gradient)
    };
}
pub(crate) use calculate_first_piola_kirchoff_tangent_stiffness_from_deformation_gradient_simple;

macro_rules! calculate_second_piola_kirchoff_stress_from_deformation_gradient {
    ($constitutive_model_constructed: expr, $deformation_gradient: expr) => {
        $constitutive_model_constructed
            .calculate_second_piola_kirchoff_stress($deformation_gradient)
    };
}
pub(crate) use calculate_second_piola_kirchoff_stress_from_deformation_gradient;

macro_rules! calculate_second_piola_kirchoff_stress_from_deformation_gradient_simple {
    ($constitutive_model_constructed: expr, $deformation_gradient: expr) => {
        $constitutive_model_constructed
            .calculate_second_piola_kirchoff_stress($deformation_gradient)
    };
}
pub(crate) use calculate_second_piola_kirchoff_stress_from_deformation_gradient_simple;

macro_rules! calculate_second_piola_kirchoff_stress_from_deformation_gradient_rotated {
    ($constitutive_model_constructed: expr, $deformation_gradient: expr) => {
        $constitutive_model_constructed
            .calculate_second_piola_kirchoff_stress($deformation_gradient)
    };
}
pub(crate) use calculate_second_piola_kirchoff_stress_from_deformation_gradient_rotated;

macro_rules! calculate_second_piola_kirchoff_tangent_stiffness_from_deformation_gradient {
    ($constitutive_model_constructed: expr, $deformation_gradient: expr) => {
        $constitutive_model_constructed
            .calculate_second_piola_kirchoff_tangent_stiffness($deformation_gradient)
    };
}
pub(crate) use calculate_second_piola_kirchoff_tangent_stiffness_from_deformation_gradient;

macro_rules! test_solid_constitutive_model {
    ($constitutive_model: ident, $constitutive_model_parameters: expr, $constitutive_model_constructed: expr) => {
        crate::constitutive::solid::elastic::test::test_solid_constitutive_construction!(
            $constitutive_model,
            $constitutive_model_parameters,
            $constitutive_model_constructed
        );
        crate::constitutive::solid::elastic::test::test_constructed_solid_constitutive_model!(
            $constitutive_model_constructed
        );
    };
}
pub(crate) use test_solid_constitutive_model;

macro_rules! test_constructed_solid_constitutive_model {
    ($constitutive_model_constructed: expr) => {
        crate::constitutive::solid::elastic::test::test_solid_constitutive_model_no_tangents!(
            $constitutive_model_constructed
        );
        crate::constitutive::solid::elastic::test::test_solid_constitutive_model_tangents!(
            $constitutive_model_constructed
        );
    };
}
pub(crate) use test_constructed_solid_constitutive_model;

macro_rules! test_solid_constitutive_construction
{
    ($constitutive_model: ident, $constitutive_model_parameters: expr, $constitutive_model_constructed: expr) =>
    {
        fn get_constitutive_model<'a>() -> $constitutive_model<'a>
        {
            $constitutive_model::new($constitutive_model_parameters)
        }
        #[test]
        fn get_bulk_modulus() -> Result<(), TestError>
        {
            assert_eq(get_constitutive_model().get_bulk_modulus(), &$constitutive_model_parameters[0])
        }
        #[test]
        fn get_shear_modulus() -> Result<(), TestError>
        {
            assert_eq(get_constitutive_model().get_shear_modulus(), &$constitutive_model_parameters[1])
        }
        #[test]
        fn bulk_modulus() -> Result<(), TestError>
        {
            let model = get_constitutive_model();
            let deformation_gradient = DeformationGradient::identity()*(1.0 + crate::EPSILON/3.0);
            let first_piola_kirchoff_stress = calculate_first_piola_kirchoff_stress_from_deformation_gradient_simple!(&model, &deformation_gradient)?;
            assert!((3.0*crate::EPSILON*model.get_bulk_modulus()/first_piola_kirchoff_stress.trace() - 1.0).abs() < crate::EPSILON);
            Ok(())
        }
        #[test]
        fn shear_modulus() -> Result<(), TestError>
        {
            let model = get_constitutive_model();
            let mut deformation_gradient = DeformationGradient::identity();
            deformation_gradient[0][1] = crate::EPSILON;
            let first_piola_kirchoff_stress = calculate_first_piola_kirchoff_stress_from_deformation_gradient_simple!(&model, &deformation_gradient)?;
            assert!((crate::EPSILON*model.get_shear_modulus()/first_piola_kirchoff_stress[0][1] - 1.0).abs() < crate::EPSILON);
            Ok(())
        }
        #[test]
        fn size()
        {
            assert_eq!(
                std::mem::size_of::<$constitutive_model>(),
                std::mem::size_of::<crate::constitutive::Parameters>()
            )
        }
    }
}
pub(crate) use test_solid_constitutive_construction;

macro_rules! test_solid_constitutive_model_no_tangents
{
    ($constitutive_model_constructed: expr) =>
    {
        use crate::
        {
            EPSILON,
            math::test::{assert_eq, assert_eq_within_tols as assert_eq_within_tols_new, TestError},
            mechanics::
            {
                CauchyStress,
                FirstPiolaKirchoffStress,
                SecondPiolaKirchoffStress,
                test::
                {
                    get_deformation_gradient,
                    get_deformation_gradient_rotated,
                    get_rotation_current_configuration,
                    get_rotation_reference_configuration
                }
            }
        };
        mod cauchy_stress
        {
            use super::*;
            mod deformed
            {
                use super::*;
                #[test]
                #[should_panic(expected = "Invalid Jacobian")]
                fn invalid_jacobian()
                {
                    let mut deformation_gradient = DeformationGradient::identity();
                    deformation_gradient[0][0] *= -1.0;
                    calculate_cauchy_stress_from_deformation_gradient!(
                        $constitutive_model_constructed, &deformation_gradient
                    ).unwrap();
                }
                #[test]
                fn objectivity() -> Result<(), TestError>
                {
                    let _ = EPSILON;
                    assert_eq_within_tols_new(
                        &calculate_cauchy_stress_from_deformation_gradient!(
                            &$constitutive_model_constructed, &get_deformation_gradient()
                        )?, &(get_rotation_current_configuration().transpose() *
                        calculate_cauchy_stress_from_deformation_gradient_rotated!(
                            &$constitutive_model_constructed, &get_deformation_gradient_rotated()
                        )? * get_rotation_current_configuration())
                    )
                }
                #[test]
                fn symmetry() -> Result<(), TestError>
                {
                    let cauchy_stress = calculate_cauchy_stress_from_deformation_gradient!(&$constitutive_model_constructed, &get_deformation_gradient())?;
                    assert_eq_within_tols_new(&cauchy_stress, &cauchy_stress.transpose())
                }
            }
            mod undeformed
            {
                use super::*;
                #[test]
                fn zero() -> Result<(), TestError>
                {
                    assert_eq_within_tols_new(
                        &calculate_cauchy_stress_from_deformation_gradient_simple!(
                            &$constitutive_model_constructed, &DeformationGradient::identity()
                        )?, &CauchyStress::zero()
                    )
                }
            }
        }
        mod first_piola_kirchoff_stress
        {
            use super::*;
            mod deformed
            {
                #[test]
                #[should_panic(expected = "Invalid Jacobian")]
                fn invalid_jacobian()
                {
                    let mut deformation_gradient = DeformationGradient::identity();
                    deformation_gradient[0][0] *= -1.0;
                    calculate_first_piola_kirchoff_stress_from_deformation_gradient!(
                        $constitutive_model_constructed, &deformation_gradient
                    ).unwrap();
                }
                use super::*;
                #[test]
                fn objectivity() -> Result<(), TestError>
                {
                    assert_eq_within_tols_new(
                        &calculate_first_piola_kirchoff_stress_from_deformation_gradient!(
                            &$constitutive_model_constructed, &get_deformation_gradient()
                        )?,
                        &(get_rotation_current_configuration().transpose() *
                        calculate_first_piola_kirchoff_stress_from_deformation_gradient_rotated!(
                            &$constitutive_model_constructed, &get_deformation_gradient_rotated()
                        )? * get_rotation_reference_configuration())
                    )
                }
            }
            mod undeformed
            {
                use super::*;
                #[test]
                fn zero() -> Result<(), TestError>
                {
                    assert_eq(
                        &calculate_first_piola_kirchoff_stress_from_deformation_gradient_simple!(
                            &$constitutive_model_constructed, &DeformationGradient::identity()
                        )?, &FirstPiolaKirchoffStress::zero()
                    )
                }
            }
        }
        mod second_piola_kirchoff_stress
        {
            use super::*;
            mod deformed
            {
                #[test]
                #[should_panic(expected = "Invalid Jacobian")]
                fn invalid_jacobian()
                {
                    let mut deformation_gradient = DeformationGradient::identity();
                    deformation_gradient[0][0] *= -1.0;
                    calculate_second_piola_kirchoff_stress_from_deformation_gradient!(
                        $constitutive_model_constructed, &deformation_gradient
                    ).unwrap();
                }
                use super::*;
                #[test]
                fn objectivity() -> Result<(), TestError>
                {
                    assert_eq_within_tols_new(
                        &calculate_second_piola_kirchoff_stress_from_deformation_gradient!(
                            &$constitutive_model_constructed, &get_deformation_gradient()
                        )?,
                        &(get_rotation_reference_configuration().transpose() *
                        calculate_second_piola_kirchoff_stress_from_deformation_gradient_rotated!(
                            &$constitutive_model_constructed, &get_deformation_gradient_rotated()
                        )? * get_rotation_reference_configuration())
                    )
                }
            }
            mod undeformed
            {
                use super::*;
                #[test]
                fn zero() -> Result<(), TestError>
                {
                    assert_eq(
                        &calculate_second_piola_kirchoff_stress_from_deformation_gradient_simple!(
                            &$constitutive_model_constructed, &DeformationGradient::identity()
                        )?, &SecondPiolaKirchoffStress::zero()
                    )
                }
            }
        }
    }
}
pub(crate) use test_solid_constitutive_model_no_tangents;

macro_rules! test_solid_constitutive_model_tangents
{
    ($constitutive_model_constructed: expr) =>
    {
        mod tangents
        {
            use crate::
            {
                math::ContractAllIndicesWithFirstIndicesOf,
                mechanics::test::get_deformation_gradient_rotated_undeformed
            };
            use super::*;
            fn calculate_cauchy_tangent_stiffness_from_finite_difference_of_cauchy_stress(is_deformed: bool) -> Result<CauchyTangentStiffness, TestError>
            {
                let model = $constitutive_model_constructed;
                let mut cauchy_tangent_stiffness = CauchyTangentStiffness::zero();
                for k in 0..3
                {
                    for l in 0..3
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
                        deformation_gradient_plus[k][l] += 0.5*EPSILON;
                        let calculate_cauchy_stress_plus = calculate_cauchy_stress_from_deformation_gradient!(&model, &deformation_gradient_plus)?;
                        let mut deformation_gradient_minus =
                            if is_deformed
                            {
                                get_deformation_gradient()
                            }
                            else
                            {
                                DeformationGradient::identity()
                            };
                        deformation_gradient_minus[k][l] -= 0.5*EPSILON;
                        let calculate_cauchy_stress_minus = calculate_cauchy_stress_from_deformation_gradient!(&model, &deformation_gradient_minus)?;
                        for i in 0..3
                        {
                            for j in 0..3
                            {
                                cauchy_tangent_stiffness[i][j][k][l] = (calculate_cauchy_stress_plus[i][j] - calculate_cauchy_stress_minus[i][j])/EPSILON;
                            }
                        }
                    }
                }
                Ok(cauchy_tangent_stiffness)
            }
            fn calculate_first_piola_kirchoff_tangent_stiffness_from_finite_difference_of_first_piola_kirchoff_stress(is_deformed: bool) -> Result<FirstPiolaKirchoffTangentStiffness, TestError>
            {
                let model = $constitutive_model_constructed;
                let mut first_piola_kirchoff_tangent_stiffness = FirstPiolaKirchoffTangentStiffness::zero();
                for k in 0..3
                {
                    for l in 0..3
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
                        deformation_gradient_plus[k][l] += 0.5*EPSILON;
                        let first_piola_kirchoff_stress_plus = calculate_first_piola_kirchoff_stress_from_deformation_gradient!(&model, &deformation_gradient_plus)?;
                        let mut deformation_gradient_minus =
                            if is_deformed
                            {
                                get_deformation_gradient()
                            }
                            else
                            {
                                DeformationGradient::identity()
                            };
                        deformation_gradient_minus[k][l] -= 0.5*EPSILON;
                        let first_piola_kirchoff_stress_minus = calculate_first_piola_kirchoff_stress_from_deformation_gradient!(&model, &deformation_gradient_minus)?;
                        for i in 0..3
                        {
                            for j in 0..3
                            {
                                first_piola_kirchoff_tangent_stiffness[i][j][k][l] = (first_piola_kirchoff_stress_plus[i][j] - first_piola_kirchoff_stress_minus[i][j])/EPSILON;
                            }
                        }
                    }
                }
                Ok(first_piola_kirchoff_tangent_stiffness)
            }
            fn calculate_second_piola_kirchoff_tangent_stiffness_from_finite_difference_of_second_piola_kirchoff_stress(is_deformed: bool) -> Result<SecondPiolaKirchoffTangentStiffness, TestError>
            {
                let mut second_piola_kirchoff_tangent_stiffness = SecondPiolaKirchoffTangentStiffness::zero();
                for k in 0..3
                {
                    for l in 0..3
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
                        deformation_gradient_plus[k][l] += 0.5*EPSILON;
                        let second_piola_kirchoff_stress_plus =
                        calculate_second_piola_kirchoff_stress_from_deformation_gradient!(
                            $constitutive_model_constructed, &deformation_gradient_plus
                        )?;
                        let mut deformation_gradient_minus =
                            if is_deformed
                            {
                                get_deformation_gradient()
                            }
                            else
                            {
                                DeformationGradient::identity()
                            };
                        deformation_gradient_minus[k][l] -= 0.5*EPSILON;
                        let second_piola_kirchoff_stress_minus =
                        calculate_second_piola_kirchoff_stress_from_deformation_gradient!(
                            $constitutive_model_constructed, &deformation_gradient_minus
                        )?;
                        for i in 0..3
                        {
                            for j in 0..3
                            {
                                second_piola_kirchoff_tangent_stiffness[i][j][k][l] = (
                                    second_piola_kirchoff_stress_plus[i][j] - second_piola_kirchoff_stress_minus[i][j]
                                )/EPSILON;
                            }
                        }
                    }
                }
                Ok(second_piola_kirchoff_tangent_stiffness)
            }
            mod cauchy_stress
            {
                use super::*;
                mod deformed
                {
                    use super::*;
                    #[test]
                    fn finite_difference() -> Result<(), TestError>
                    {
                        calculate_cauchy_tangent_stiffness_from_deformation_gradient!(
                            &$constitutive_model_constructed, &get_deformation_gradient()
                        )?.iter().zip(calculate_cauchy_tangent_stiffness_from_finite_difference_of_cauchy_stress(true)?.iter())
                        .for_each(|(cauchy_tangent_stiffness_i, fd_cauchy_tangent_stiffness_i)|
                            cauchy_tangent_stiffness_i.iter()
                            .zip(fd_cauchy_tangent_stiffness_i.iter())
                            .for_each(|(cauchy_tangent_stiffness_ij, fd_cauchy_tangent_stiffness_ij)|
                                cauchy_tangent_stiffness_ij.iter()
                                .zip(fd_cauchy_tangent_stiffness_ij.iter())
                                .for_each(|(cauchy_tangent_stiffness_ijk, fd_cauchy_tangent_stiffness_ijk)|
                                    cauchy_tangent_stiffness_ijk.iter()
                                    .zip(fd_cauchy_tangent_stiffness_ijk.iter())
                                    .for_each(|(cauchy_tangent_stiffness_ijkl, fd_cauchy_tangent_stiffness_ijkl)|
                                        assert!((cauchy_tangent_stiffness_ijkl/fd_cauchy_tangent_stiffness_ijkl - 1.0).abs() < EPSILON)
                                    )
                                )
                            )
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
                        calculate_cauchy_tangent_stiffness_from_deformation_gradient!(
                            &$constitutive_model_constructed, &DeformationGradient::identity()
                        )?.iter().zip(calculate_cauchy_tangent_stiffness_from_finite_difference_of_cauchy_stress(false)?.iter())
                        .for_each(|(cauchy_tangent_stiffness_i, fd_cauchy_tangent_stiffness_i)|
                            cauchy_tangent_stiffness_i.iter()
                            .zip(fd_cauchy_tangent_stiffness_i.iter())
                            .for_each(|(cauchy_tangent_stiffness_ij, fd_cauchy_tangent_stiffness_ij)|
                                cauchy_tangent_stiffness_ij.iter()
                                .zip(fd_cauchy_tangent_stiffness_ij.iter())
                                .for_each(|(cauchy_tangent_stiffness_ijk, fd_cauchy_tangent_stiffness_ijk)|
                                    cauchy_tangent_stiffness_ijk.iter()
                                    .zip(fd_cauchy_tangent_stiffness_ijk.iter())
                                    .for_each(|(cauchy_tangent_stiffness_ijkl, fd_cauchy_tangent_stiffness_ijkl)|
                                        assert!(
                                            (cauchy_tangent_stiffness_ijkl/fd_cauchy_tangent_stiffness_ijkl - 1.0).abs() < EPSILON ||
                                            fd_cauchy_tangent_stiffness_ijkl.abs() < EPSILON
                                        )
                                    )
                                )
                            )
                        );
                        Ok(())
                    }
                }
            }
            mod first_piola_kirchoff_stress
            {
                use super::*;
                mod deformed
                {
                    use super::*;
                    #[test]
                    fn finite_difference() -> Result<(), TestError>
                    {
                        calculate_first_piola_kirchoff_tangent_stiffness_from_deformation_gradient!(
                            &$constitutive_model_constructed, &get_deformation_gradient()
                        )?.iter()
                        .zip(calculate_first_piola_kirchoff_tangent_stiffness_from_finite_difference_of_first_piola_kirchoff_stress(true)?.iter())
                        .for_each(|(first_piola_kirchoff_tangent_stiffness_i, fd_first_piola_kirchoff_tangent_stiffness_i)|
                            first_piola_kirchoff_tangent_stiffness_i.iter()
                            .zip(fd_first_piola_kirchoff_tangent_stiffness_i.iter())
                            .for_each(|(first_piola_kirchoff_tangent_stiffness_ij, fd_first_piola_kirchoff_tangent_stiffness_ij)|
                                first_piola_kirchoff_tangent_stiffness_ij.iter()
                                .zip(fd_first_piola_kirchoff_tangent_stiffness_ij.iter())
                                .for_each(|(first_piola_kirchoff_tangent_stiffness_ijk, fd_first_piola_kirchoff_tangent_stiffness_ijk)|
                                    first_piola_kirchoff_tangent_stiffness_ijk.iter()
                                    .zip(fd_first_piola_kirchoff_tangent_stiffness_ijk.iter())
                                    .for_each(|(first_piola_kirchoff_tangent_stiffness_ijkl, fd_first_piola_kirchoff_tangent_stiffness_ijkl)|
                                        assert!((first_piola_kirchoff_tangent_stiffness_ijkl/fd_first_piola_kirchoff_tangent_stiffness_ijkl - 1.0).abs() < EPSILON)
                                    )
                                )
                            )
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
                        calculate_first_piola_kirchoff_tangent_stiffness_from_deformation_gradient!(
                            &$constitutive_model_constructed, &DeformationGradient::identity()
                        )?.iter()
                        .zip(calculate_first_piola_kirchoff_tangent_stiffness_from_finite_difference_of_first_piola_kirchoff_stress(false)?.iter())
                        .for_each(|(first_piola_kirchoff_tangent_stiffness_i, fd_first_piola_kirchoff_tangent_stiffness_i)|
                            first_piola_kirchoff_tangent_stiffness_i.iter()
                            .zip(fd_first_piola_kirchoff_tangent_stiffness_i.iter())
                            .for_each(|(first_piola_kirchoff_tangent_stiffness_ij, fd_first_piola_kirchoff_tangent_stiffness_ij)|
                                first_piola_kirchoff_tangent_stiffness_ij.iter()
                                .zip(fd_first_piola_kirchoff_tangent_stiffness_ij.iter())
                                .for_each(|(first_piola_kirchoff_tangent_stiffness_ijk, fd_first_piola_kirchoff_tangent_stiffness_ijk)|
                                    first_piola_kirchoff_tangent_stiffness_ijk.iter()
                                    .zip(fd_first_piola_kirchoff_tangent_stiffness_ijk.iter())
                                    .for_each(|(first_piola_kirchoff_tangent_stiffness_ijkl, fd_first_piola_kirchoff_tangent_stiffness_ijkl)|
                                        assert!(
                                            (first_piola_kirchoff_tangent_stiffness_ijkl/fd_first_piola_kirchoff_tangent_stiffness_ijkl - 1.0).abs() < EPSILON ||
                                            fd_first_piola_kirchoff_tangent_stiffness_ijkl.abs() < EPSILON
                                        )
                                    )
                                )
                            )
                        );
                        Ok(())
                    }
                }
            }
            mod second_piola_kirchoff_stress
            {
                use super::*;
                mod deformed
                {
                    use super::*;
                    #[test]
                    fn finite_difference() -> Result<(), TestError>
                    {
                        calculate_second_piola_kirchoff_tangent_stiffness_from_deformation_gradient!(
                            &$constitutive_model_constructed, &get_deformation_gradient()
                        )?.iter()
                        .zip(calculate_second_piola_kirchoff_tangent_stiffness_from_finite_difference_of_second_piola_kirchoff_stress(true)?.iter())
                        .for_each(|(second_piola_kirchoff_tangent_stiffness_i, fd_second_piola_kirchoff_tangent_stiffness_i)|
                            second_piola_kirchoff_tangent_stiffness_i.iter()
                            .zip(fd_second_piola_kirchoff_tangent_stiffness_i.iter())
                            .for_each(|(second_piola_kirchoff_tangent_stiffness_ij, fd_second_piola_kirchoff_tangent_stiffness_ij)|
                                second_piola_kirchoff_tangent_stiffness_ij.iter()
                                .zip(fd_second_piola_kirchoff_tangent_stiffness_ij.iter())
                                .for_each(|(second_piola_kirchoff_tangent_stiffness_ijk, fd_second_piola_kirchoff_tangent_stiffness_ijk)|
                                    second_piola_kirchoff_tangent_stiffness_ijk.iter()
                                    .zip(fd_second_piola_kirchoff_tangent_stiffness_ijk.iter())
                                    .for_each(|(second_piola_kirchoff_tangent_stiffness_ijkl, fd_second_piola_kirchoff_tangent_stiffness_ijkl)|
                                        assert!(
                                            (second_piola_kirchoff_tangent_stiffness_ijkl/fd_second_piola_kirchoff_tangent_stiffness_ijkl - 1.0).abs() < EPSILON ||
                                            fd_second_piola_kirchoff_tangent_stiffness_ijkl.abs() < EPSILON
                                        )
                                    )
                                )
                            )
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
                        calculate_second_piola_kirchoff_tangent_stiffness_from_deformation_gradient!(
                            &$constitutive_model_constructed, &DeformationGradient::identity()
                        )?.iter()
                        .zip(calculate_second_piola_kirchoff_tangent_stiffness_from_finite_difference_of_second_piola_kirchoff_stress(false)?.iter())
                        .for_each(|(second_piola_kirchoff_tangent_stiffness_i, fd_second_piola_kirchoff_tangent_stiffness_i)|
                            second_piola_kirchoff_tangent_stiffness_i.iter()
                            .zip(fd_second_piola_kirchoff_tangent_stiffness_i.iter())
                            .for_each(|(second_piola_kirchoff_tangent_stiffness_ij, fd_second_piola_kirchoff_tangent_stiffness_ij)|
                                second_piola_kirchoff_tangent_stiffness_ij.iter()
                                .zip(fd_second_piola_kirchoff_tangent_stiffness_ij.iter())
                                .for_each(|(second_piola_kirchoff_tangent_stiffness_ijk, fd_second_piola_kirchoff_tangent_stiffness_ijk)|
                                    second_piola_kirchoff_tangent_stiffness_ijk.iter()
                                    .zip(fd_second_piola_kirchoff_tangent_stiffness_ijk.iter())
                                    .for_each(|(second_piola_kirchoff_tangent_stiffness_ijkl, fd_second_piola_kirchoff_tangent_stiffness_ijkl)|
                                        assert!(
                                            (second_piola_kirchoff_tangent_stiffness_ijkl/fd_second_piola_kirchoff_tangent_stiffness_ijkl - 1.0).abs() < EPSILON ||
                                            fd_second_piola_kirchoff_tangent_stiffness_ijkl.abs() < EPSILON
                                        )
                                    )
                                )
                            )
                        );
                        Ok(())
                    }
                }
            }
            mod cauchy_tangent_stiffness
            {
                use super::*;
                mod deformed
                {
                    use super::*;
                    #[test]
                    #[should_panic(expected = "Invalid Jacobian")]
                    fn invalid_jacobian()
                    {
                        let mut deformation_gradient = DeformationGradient::identity();
                        deformation_gradient[0][0] *= -1.0;
                        calculate_cauchy_tangent_stiffness_from_deformation_gradient!(
                            $constitutive_model_constructed, &deformation_gradient
                        ).unwrap();
                    }
                    #[test]
                    fn objectivity() -> Result<(), TestError>
                    {
                        assert_eq_within_tols_new(
                            &calculate_cauchy_tangent_stiffness_from_deformation_gradient!(
                                &$constitutive_model_constructed, &get_deformation_gradient()
                            )?,
                            &calculate_cauchy_tangent_stiffness_from_deformation_gradient!(
                                &$constitutive_model_constructed, &get_deformation_gradient_rotated()
                            )?.contract_all_indices_with_first_indices_of(
                                &get_rotation_current_configuration(),
                                &get_rotation_current_configuration(),
                                &get_rotation_current_configuration(),
                                &get_rotation_reference_configuration()
                            )
                        )
                    }
                    #[test]
                    fn symmetry() -> Result<(), TestError>
                    {
                        let cauchy_tangent_stiffness =
                        calculate_cauchy_tangent_stiffness_from_deformation_gradient!(
                            &$constitutive_model_constructed, &get_deformation_gradient()
                        )?;
                        assert_eq_within_tols_new(
                            &cauchy_tangent_stiffness,
                            &(0..3).map(|i|
                                (0..3).map(|j|
                                    cauchy_tangent_stiffness[j][i].copy()
                                ).collect()
                            ).collect()
                        )
                    }
                }
                mod undeformed
                {
                    use super::*;
                    #[test]
                    fn objectivity() -> Result<(), TestError>
                    {
                        assert_eq_within_tols_new(
                            &calculate_cauchy_tangent_stiffness_from_deformation_gradient!(
                                &$constitutive_model_constructed, &DeformationGradient::identity()
                            )?,
                            &calculate_cauchy_tangent_stiffness_from_deformation_gradient!(
                                &$constitutive_model_constructed, &get_deformation_gradient_rotated_undeformed()
                            )?.contract_all_indices_with_first_indices_of(
                                &get_rotation_current_configuration(),
                                &get_rotation_current_configuration(),
                                &get_rotation_current_configuration(),
                                &get_rotation_reference_configuration()
                            )
                        )
                    }
                    #[test]
                    fn symmetry() -> Result<(), TestError>
                    {
                        let cauchy_tangent_stiffness =
                        calculate_cauchy_tangent_stiffness_from_deformation_gradient!(
                            &$constitutive_model_constructed, &DeformationGradient::identity()
                        )?;
                        assert_eq_within_tols_new(
                            &cauchy_tangent_stiffness,
                            &(0..3).map(|i|
                                (0..3).map(|j|
                                    cauchy_tangent_stiffness[j][i].copy()
                                ).collect()
                            ).collect()
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
                    #[should_panic(expected = "Invalid Jacobian")]
                    fn invalid_jacobian()
                    {
                        let mut deformation_gradient = DeformationGradient::identity();
                        deformation_gradient[0][0] *= -1.0;
                        calculate_first_piola_kirchoff_tangent_stiffness_from_deformation_gradient!(
                            $constitutive_model_constructed, &deformation_gradient
                        ).unwrap();
                    }
                    #[test]
                    fn objectivity() -> Result<(), TestError>
                    {
                        assert_eq_within_tols_new(
                            &calculate_first_piola_kirchoff_tangent_stiffness_from_deformation_gradient!(
                                &$constitutive_model_constructed, &get_deformation_gradient()
                            )?,
                            &calculate_first_piola_kirchoff_tangent_stiffness_from_deformation_gradient!(
                                &$constitutive_model_constructed, &get_deformation_gradient_rotated()
                            )?.contract_all_indices_with_first_indices_of(
                                &get_rotation_current_configuration(),
                                &get_rotation_reference_configuration(),
                                &get_rotation_current_configuration(),
                                &get_rotation_reference_configuration()
                            )
                        )
                    }
                }
                mod undeformed
                {
                    use super::*;
                    #[test]
                    fn objectivity() -> Result<(), TestError>
                    {
                        assert_eq_within_tols_new(
                            &calculate_first_piola_kirchoff_tangent_stiffness_from_deformation_gradient!(
                                &$constitutive_model_constructed, &DeformationGradient::identity()
                            )?,
                            &calculate_first_piola_kirchoff_tangent_stiffness_from_deformation_gradient!(
                                &$constitutive_model_constructed, &get_deformation_gradient_rotated_undeformed()
                            )?.contract_all_indices_with_first_indices_of(
                                &get_rotation_current_configuration(),
                                &get_rotation_reference_configuration(),
                                &get_rotation_current_configuration(),
                                &get_rotation_reference_configuration()
                            )
                        )
                    }
                }
            }
            mod second_piola_kirchoff_tangent_stiffness
            {
                use super::*;
                mod deformed
                {
                    use super::*;
                    #[test]
                    #[should_panic(expected = "Invalid Jacobian")]
                    fn invalid_jacobian()
                    {
                        let mut deformation_gradient = DeformationGradient::identity();
                        deformation_gradient[0][0] *= -1.0;
                        calculate_second_piola_kirchoff_tangent_stiffness_from_deformation_gradient!(
                            $constitutive_model_constructed, &deformation_gradient
                        ).unwrap();
                    }
                    #[test]
                    fn objectivity() -> Result<(), TestError>
                    {
                        assert_eq_within_tols_new(
                            &calculate_second_piola_kirchoff_tangent_stiffness_from_deformation_gradient!(
                                &$constitutive_model_constructed, &get_deformation_gradient()
                            )?,
                            &calculate_second_piola_kirchoff_tangent_stiffness_from_deformation_gradient!(
                                &$constitutive_model_constructed, &get_deformation_gradient_rotated()
                            )?.contract_all_indices_with_first_indices_of(
                                &get_rotation_reference_configuration(),
                                &get_rotation_reference_configuration(),
                                &get_rotation_current_configuration(),
                                &get_rotation_reference_configuration()
                            )
                        )
                    }
                }
                mod undeformed
                {
                    use super::*;
                    #[test]
                    fn objectivity() -> Result<(), TestError>
                    {
                        assert_eq_within_tols_new(
                            &calculate_second_piola_kirchoff_tangent_stiffness_from_deformation_gradient!(
                                &$constitutive_model_constructed, &DeformationGradient::identity()
                            )?,
                            &calculate_second_piola_kirchoff_tangent_stiffness_from_deformation_gradient!(
                                &$constitutive_model_constructed, &get_deformation_gradient_rotated_undeformed()
                            )?.contract_all_indices_with_first_indices_of(
                                &get_rotation_reference_configuration(),
                                &get_rotation_reference_configuration(),
                                &get_rotation_current_configuration(),
                                &get_rotation_reference_configuration()
                            )
                        )
                    }
                }
            }
        }
    }
}
pub(crate) use test_solid_constitutive_model_tangents;

macro_rules! test_solid_elastic_constitutive_model
{
    ($constitutive_model: ident, $constitutive_model_parameters: expr, $constitutive_model_constructed: expr) =>
    {
        crate::constitutive::solid::elastic::test::test_solid_constitutive_model!(
            $constitutive_model,
            $constitutive_model_parameters,
            $constitutive_model_constructed
        );
        mod elastic
        {
            use super::*;
            mod first_piola_kirchoff_tangent_stiffness
            {
                use super::*;
                mod deformed
                {
                    use super::*;
                    #[test]
                    fn non_symmetry() -> Result<(), TestError>
                    {
                        let first_piola_kirchoff_tangent_stiffness =
                        calculate_first_piola_kirchoff_tangent_stiffness_from_deformation_gradient!(
                            &$constitutive_model_constructed, &get_deformation_gradient()
                        )?;
                        assert!(
                            assert_eq_within_tols_new(
                                &first_piola_kirchoff_tangent_stiffness,
                                &(0..3).map(|i|
                                    (0..3).map(|j|
                                        (0..3).map(|k|
                                            (0..3).map(|l|
                                                first_piola_kirchoff_tangent_stiffness[k][l][i][j].copy()
                                            ).collect()
                                        ).collect()
                                    ).collect()
                                ).collect()
                            ).is_err()
                        );
                        Ok(())
                    }
                }
            }
        }
    }
}
pub(crate) use test_solid_elastic_constitutive_model;
