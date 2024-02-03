use crate::
{
    EPSILON,
    constitutive::solid::
    {
        hyperelastic::
        {
            Gent,
            MooneyRivlin,
            NeoHookean,
            Yeoh,
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
        CombinedHyperelastic::construct(
            NeoHookean::new(NEOHOOKEANPARAMETERS),
            NeoHookean::new(NEOHOOKEANPARAMETERS)
        )
    );
}

mod mixed
{
    use super::*;
    test_hyperelastic_constitutive_model_constructed!(
        CombinedHyperelastic::construct(
            Gent::new(GENTPARAMETERS),
            MooneyRivlin::new(MOONEYRIVLINPARAMETERS)
        )
    );
}

mod nested
{
    use super::*;
    test_hyperelastic_constitutive_model_constructed!(
        CombinedHyperelastic::construct(
            CombinedHyperelastic::construct(
                Gent::new(GENTPARAMETERS),
                Yeoh::new(YEOHPARAMETERS)
            ),
            CombinedHyperelastic::construct(
                MooneyRivlin::new(MOONEYRIVLINPARAMETERS),
                NeoHookean::new(NEOHOOKEANPARAMETERS)
            )
        )
    );
}

#[test]
fn additive_cauchy_stress()
{
    CombinedHyperelastic::construct(
        Gent::new(GENTPARAMETERS),
        MooneyRivlin::new(MOONEYRIVLINPARAMETERS)
    ).calculate_cauchy_stress(&get_deformation_gradient()).iter()
    .zip(
        MooneyRivlin::new(MOONEYRIVLINPARAMETERS).calculate_cauchy_stress(&get_deformation_gradient()).iter()
        .zip(Gent::new(GENTPARAMETERS).calculate_cauchy_stress(&get_deformation_gradient()).iter())
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
    CombinedHyperelastic::construct(
        Gent::new(GENTPARAMETERS),
        MooneyRivlin::new(MOONEYRIVLINPARAMETERS)
    ).calculate_cauchy_tangent_stiffness(&get_deformation_gradient()).iter()
    .zip(
        MooneyRivlin::new(MOONEYRIVLINPARAMETERS).calculate_cauchy_tangent_stiffness(&get_deformation_gradient()).iter()
        .zip(Gent::new(GENTPARAMETERS).calculate_cauchy_tangent_stiffness(&get_deformation_gradient()).iter())
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
        CombinedHyperelastic::construct(
            Gent::new(GENTPARAMETERS),
            MooneyRivlin::new(MOONEYRIVLINPARAMETERS)
        ).calculate_helmholtz_free_energy_density(&get_deformation_gradient()),
        MooneyRivlin::new(MOONEYRIVLINPARAMETERS).calculate_helmholtz_free_energy_density(&get_deformation_gradient())
        + Gent::new(GENTPARAMETERS).calculate_helmholtz_free_energy_density(&get_deformation_gradient())
    )
}
