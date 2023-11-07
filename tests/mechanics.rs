#[cfg(feature = "mechanics")]
mod public
{
    use flavio::mechanics::
    {
        CauchyStress,
        CauchyTangentStiffness,
        DeformationGradient,
        FirstPiolaKirchoffStress,
        FirstPiolaKirchoffTangentStiffness,
        LeftCauchyGreenDeformation,
        RightCauchyGreenDeformation,
        Scalar
    };
    #[test]
    fn cauchy_stress()
    {
        let _: CauchyStress;
    }
    #[test]
    fn cauchy_tangent_stiffness()
    {
        let _: CauchyTangentStiffness;
    }
    #[test]
    fn deformation_gradient()
    {
        let _: DeformationGradient;
    }
    #[test]
    fn first_piola_kirchoff_stress()
    {
        let _: FirstPiolaKirchoffStress;
    }
    #[test]
    fn first_piola_kirchoff_tangent_stiffness()
    {
        let _: FirstPiolaKirchoffTangentStiffness;
    }
    #[test]
    fn left_cauchy_green_deformation()
    {
        let _: LeftCauchyGreenDeformation;
    }
    #[test]
    fn right_cauchy_green_deformation()
    {
        let _: RightCauchyGreenDeformation;
    }
    #[test]
    fn scalar()
    {
        let _: Scalar;
    }
}