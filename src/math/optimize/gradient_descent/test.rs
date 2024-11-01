use super::{FirstOrder, GradientDescent, TensorRank0};

const TOLERANCE: TensorRank0 = 1e-6;

#[test]
fn linear() {
    assert!(
        GradientDescent {
            ..Default::default()
        }
        .minimize(|x: &TensorRank0| *x, 1.0,)
        .unwrap()
        .abs()
            < TOLERANCE
    )
}

#[test]
fn quadratic() {
    assert!(
        GradientDescent {
            ..Default::default()
        }
        .minimize(|x: &TensorRank0| x.powi(2) / 2.0, 1.0,)
        .unwrap()
        .abs()
            < TOLERANCE
    )
}

#[test]
fn sin() {
    assert!(
        GradientDescent {
            ..Default::default()
        }
        .minimize(|x: &TensorRank0| x.sin(), 1.0,)
        .unwrap()
        .abs()
            < TOLERANCE
    )
}
