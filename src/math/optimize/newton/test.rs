use super::{Newton, Optimize, TensorRank0};

const TOLERANCE: TensorRank0 = 1e-6;

#[test]
fn linear() {
    assert!(
        Newton {
            ..Default::default()
        }
        .minimize(|x: &TensorRank0| *x, |_: &TensorRank0| 1.0, 1.0,)
        .unwrap()
        .abs()
            < TOLERANCE
    )
}

#[test]
fn quadratic() {
    assert!(
        Newton {
            ..Default::default()
        }
        .minimize(|x: &TensorRank0| x.powi(2) / 2.0, |x: &TensorRank0| *x, 1.0,)
        .unwrap()
        .abs()
            < TOLERANCE
    )
}

#[test]
fn sin() {
    assert!(
        Newton {
            ..Default::default()
        }
        .minimize(|x: &TensorRank0| x.sin(), |x: &TensorRank0| x.cos(), 1.0,)
        .unwrap()
        .abs()
            < TOLERANCE
    )
}

#[test]
#[should_panic(expected = "The obtained solution is not a minimum.")]
fn sin_max() {
    Newton {
        ..Default::default()
    }
    .minimize(|x: &TensorRank0| x.sin(), |x: &TensorRank0| x.cos(), 3.0)
    .unwrap();
}
