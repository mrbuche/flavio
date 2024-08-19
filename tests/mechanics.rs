#[cfg(feature = "mechanics")]
mod public {
    use flavio::mechanics::{
        CauchyRateTangentStiffness, CauchyStress, CauchyTangentStiffness, DeformationGradient,
        DeformationGradientRate, FirstPiolaKirchoffRateTangentStiffness, FirstPiolaKirchoffStress,
        FirstPiolaKirchoffTangentStiffness, LeftCauchyGreenDeformation,
        RightCauchyGreenDeformation, Scalar,
    };
    #[test]
    fn cauchy_stress() {
        let _: CauchyStress;
    }
    #[test]
    fn cauchy_tangent_stiffness() {
        let _: CauchyTangentStiffness;
    }
    #[test]
    fn cauchy_rate_tangent_stiffness() {
        let _: CauchyRateTangentStiffness;
    }
    #[test]
    fn deformation_gradient() {
        let _: DeformationGradient;
    }
    #[test]
    fn deformation_gradient_rate() {
        let _: DeformationGradientRate;
    }
    #[test]
    fn first_piola_kirchoff_stress() {
        let _: FirstPiolaKirchoffStress;
    }
    #[test]
    fn first_piola_kirchoff_tangent_stiffness() {
        let _: FirstPiolaKirchoffTangentStiffness;
    }
    #[test]
    fn first_piola_kirchoff_rate_tangent_stiffness() {
        let _: FirstPiolaKirchoffRateTangentStiffness;
    }
    #[test]
    fn left_cauchy_green_deformation() {
        let _: LeftCauchyGreenDeformation;
    }
    #[test]
    fn right_cauchy_green_deformation() {
        let _: RightCauchyGreenDeformation;
    }
    #[test]
    fn scalar() {
        let _: Scalar;
    }
}
