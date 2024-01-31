use crate::
{
    EPSILON,
    constitutive::solid::
    {
        hyperelastic::
        {
            GentModel,
            MooneyRivlinModel,
            NeoHookeanModel,
            YeohModel,
            test::test_hyperelastic_constitutive_model_constructed,
            test::
            {
                GENTPARAMETERS,
                MOONEYRIVLINPARAMETERS,
                NEOHOOKEANPARAMETERS,
                YEOHPARAMETERS
            }
        }
    },
    mechanics::test::
    {
        get_deformation_gradient,
        get_deformation_gradient_rotated,
    },
    test::assert_eq_within_tols
};
use super::*;

mod dual
{
    use super::*;
    test_hyperelastic_constitutive_model_constructed!(
        CompositeHyperelasticConstitutiveModel::construct(
            NeoHookeanModel::new(NEOHOOKEANPARAMETERS),
            NeoHookeanModel::new(NEOHOOKEANPARAMETERS)
        )
    );
}

mod mixed
{
    use super::*;
    test_hyperelastic_constitutive_model_constructed!(
        CompositeHyperelasticConstitutiveModel::construct(
            GentModel::new(GENTPARAMETERS),
            MooneyRivlinModel::new(MOONEYRIVLINPARAMETERS)
        )
    );
}

mod nested
{
    use super::*;
    test_hyperelastic_constitutive_model_constructed!(
        CompositeHyperelasticConstitutiveModel::construct(
            CompositeHyperelasticConstitutiveModel::construct(
                GentModel::new(GENTPARAMETERS),
                YeohModel::new(YEOHPARAMETERS)
            ),
            CompositeHyperelasticConstitutiveModel::construct(
                MooneyRivlinModel::new(MOONEYRIVLINPARAMETERS),
                NeoHookeanModel::new(NEOHOOKEANPARAMETERS)
            )
        )
    );
}

#[test]
fn additive_cauchy_stress()
{
    CompositeHyperelasticConstitutiveModel::construct(
        GentModel::new(GENTPARAMETERS),
        MooneyRivlinModel::new(MOONEYRIVLINPARAMETERS)
    ).calculate_cauchy_stress(&get_deformation_gradient()).iter()
    .zip(
        MooneyRivlinModel::new(MOONEYRIVLINPARAMETERS).calculate_cauchy_stress(&get_deformation_gradient()).iter()
        .zip(GentModel::new(GENTPARAMETERS).calculate_cauchy_stress(&get_deformation_gradient()).iter())
    ).for_each(|(cauchy_stress_i, (cauchy_stress_1_i, cauchy_stress_2_i))|
        cauchy_stress_i.iter()
        .zip(
            cauchy_stress_1_i.iter()
            .zip(cauchy_stress_2_i.iter())
        ).for_each(|(cauchy_stress_ij, (cauchy_stress_1_ij, cauchy_stress_2_ij))|
            assert_eq!(cauchy_stress_ij, &(cauchy_stress_1_ij + cauchy_stress_2_ij))
        )
    )
}

#[test]
fn additive_cauchy_tangent_stiffness()
{
    CompositeHyperelasticConstitutiveModel::construct(
        GentModel::new(GENTPARAMETERS),
        MooneyRivlinModel::new(MOONEYRIVLINPARAMETERS)
    ).calculate_cauchy_tangent_stiffness(&get_deformation_gradient()).iter()
    .zip(
        MooneyRivlinModel::new(MOONEYRIVLINPARAMETERS).calculate_cauchy_tangent_stiffness(&get_deformation_gradient()).iter()
        .zip(GentModel::new(GENTPARAMETERS).calculate_cauchy_tangent_stiffness(&get_deformation_gradient()).iter())
    ).for_each(|(cauchy_tangent_stiffness_i, (cauchy_tangent_stiffness_1_i, cauchy_tangent_stiffness_2_i))|
        cauchy_tangent_stiffness_i.iter()
        .zip(
            cauchy_tangent_stiffness_1_i.iter()
            .zip(cauchy_tangent_stiffness_2_i.iter())
        ).for_each(|(cauchy_tangent_stiffness_ij, (cauchy_tangent_stiffness_1_ij, cauchy_tangent_stiffness_2_ij))|
            cauchy_tangent_stiffness_ij.iter()
            .zip(
                cauchy_tangent_stiffness_1_ij.iter()
                .zip(cauchy_tangent_stiffness_2_ij.iter())
            ).for_each(|(cauchy_tangent_stiffness_ijk, (cauchy_tangent_stiffness_1_ijk, cauchy_tangent_stiffness_2_ijk))|
                cauchy_tangent_stiffness_ijk.iter()
                .zip(
                    cauchy_tangent_stiffness_1_ijk.iter()
                    .zip(cauchy_tangent_stiffness_2_ijk.iter())
                ).for_each(|(cauchy_tangent_stiffness_ijkl, (cauchy_tangent_stiffness_1_ijkl, cauchy_tangent_stiffness_2_ijkl))|
                    assert_eq!(cauchy_tangent_stiffness_ijkl, &(cauchy_tangent_stiffness_1_ijkl + cauchy_tangent_stiffness_2_ijkl))
                )
            )
        )
    )
}

#[test]
fn additive_helmholtz_free_energy_density()
{
    assert_eq!(
        CompositeHyperelasticConstitutiveModel::construct(
            GentModel::new(GENTPARAMETERS),
            MooneyRivlinModel::new(MOONEYRIVLINPARAMETERS)
        ).calculate_helmholtz_free_energy_density(&get_deformation_gradient()),
        MooneyRivlinModel::new(MOONEYRIVLINPARAMETERS).calculate_helmholtz_free_energy_density(&get_deformation_gradient())
        + GentModel::new(GENTPARAMETERS).calculate_helmholtz_free_energy_density(&get_deformation_gradient())
    )
}
