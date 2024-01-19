use super::Scalar;

pub const ALMANSIHAMELPARAMETERS: &[Scalar; 2] = &[NEOHOOKEANPARAMETERS[0], NEOHOOKEANPARAMETERS[1]];
pub const ARRUDABOYCEPARAMETERS: &[Scalar; 3] = &[NEOHOOKEANPARAMETERS[0], NEOHOOKEANPARAMETERS[1], 8.0];
pub const FUNGPARAMETERS: &[Scalar; 3] = &[NEOHOOKEANPARAMETERS[0], NEOHOOKEANPARAMETERS[1], 1.0];
pub const GENTPARAMETERS: &[Scalar; 3] = &[NEOHOOKEANPARAMETERS[0], NEOHOOKEANPARAMETERS[1], 23.0];
pub const MOONEYRIVLINPARAMETERS: &[Scalar; 3] = &[NEOHOOKEANPARAMETERS[0], NEOHOOKEANPARAMETERS[1], 1.0];
pub const NEOHOOKEANPARAMETERS: &[Scalar; 2] = &[13.0, 3.0];
pub const SAINTVENANTKIRCHOFFPARAMETERS: &[Scalar; 2] = &[NEOHOOKEANPARAMETERS[0], NEOHOOKEANPARAMETERS[1]];
pub const YEOHPARAMETERS: &[Scalar; 6] = &[NEOHOOKEANPARAMETERS[0], NEOHOOKEANPARAMETERS[1], -1.0, 3e-1, -1e-3, 1e-5];

#[test]
fn size()
{
    assert_eq!(
        std::mem::size_of::<super::ConstitutiveModelParameters>(), 16
    )
}
