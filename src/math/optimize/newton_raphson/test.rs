use super::{NewtonRaphson, SecondOrder, TensorRank0};

const TOLERANCE: TensorRank0 = 1e-6;

#[test]
fn linear() {
    assert!(
        NewtonRaphson {
            ..Default::default()
        }
        .minimize(|x: &TensorRank0| *x, |_: &TensorRank0| 1.0, 1.0, None)
        .unwrap()
        .abs()
            < TOLERANCE
    )
}

#[test]
fn quadratic() {
    assert!(
        NewtonRaphson {
            ..Default::default()
        }
        .minimize(
            |x: &TensorRank0| x.powi(2) / 2.0,
            |x: &TensorRank0| *x,
            1.0,
            None
        )
        .unwrap()
        .abs()
            < TOLERANCE
    )
}

#[test]
fn sin() {
    assert!(
        NewtonRaphson {
            ..Default::default()
        }
        .minimize(
            |x: &TensorRank0| x.sin(),
            |x: &TensorRank0| x.cos(),
            1.0,
            None
        )
        .unwrap()
        .abs()
            < TOLERANCE
    )
}

#[test]
#[should_panic(expected = "The obtained solution is not a minimum.")]
fn sin_max() {
    NewtonRaphson {
        ..Default::default()
    }
    .minimize(
        |x: &TensorRank0| x.sin(),
        |x: &TensorRank0| x.cos(),
        3.0,
        None,
    )
    .unwrap();
}
