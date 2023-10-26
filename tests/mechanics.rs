#[cfg(feature = "mechanics")]
mod public
{
    use flavio::mechanics::
    {
        CauchyStress,
        DeformationGradient,
        DeformationGradient1,
        DeformationGradient2,
        FirstPiolaKirchoffStress,
        Scalar
    };
    #[test]
    fn cauchy_stress()
    {
        let _: CauchyStress;
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
    fn scalar()
    {
        let _: Scalar;
    }
}