use super::
{
    ConstitutiveModel,
    Fourier,
    TemperatureGradient,
    ThermalConduction,
    super::test::FOURIERPARAMETERS
};
use crate::math::TensorRank1Trait;

#[test]
fn zero()
{
    let model = Fourier::new(FOURIERPARAMETERS);
    let zero_vector = TemperatureGradient::new([0.0, 0.0, 0.0]);
    model.calculate_heat_flux(&zero_vector).iter().for_each(|entry|
        assert_eq!(entry, &0.0)
    )
}