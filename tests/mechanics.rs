#[cfg(feature = "mechanics")]
mod public
{
    use flavio::mechanics::
    {
        CauchyStress,
        CauchyTangentStiffness,
        DeformationGradient,
        DeformationGradient1,
        DeformationGradient2,
        FirstPiolaKirchoffStress,
        FirstPiolaKirchoffTangentStiffness,
        LeftCauchyGreenDeformation,
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
    fn deformation_gradient_1()
    {
        let _: DeformationGradient1;
    }
    #[test]
    fn deformation_gradient_2()
    {
        let _: DeformationGradient2;
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
    fn scalar()
    {
        let _: Scalar;
    }
}