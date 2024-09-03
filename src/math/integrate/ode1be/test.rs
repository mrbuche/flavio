use super::{super::test::zero_to_tau, ode1be};
use crate::math::{Tensor, TensorRank0, TensorRank0List, TensorRank1, TensorRank1List, Tensors};

const LENGTH: usize = 33;
const TOLERANCE: TensorRank0 = 1e-5;

#[test]
fn first_order_tensor_rank_0() {
    let evaluation_times = zero_to_tau::<LENGTH>();
    let solution: TensorRank0List<LENGTH> =
        ode1be(|t: &f64, _: &TensorRank0| t.cos(), &evaluation_times, 0.0).expect("the unexpected");
    evaluation_times
        .iter()
        .zip(solution.iter())
        .for_each(|(t, y)| {
            assert!((t.sin() - y).abs() < TOLERANCE || (t.sin() / y - 1.0).abs() < TOLERANCE)
        });
}

#[test]
fn second_order_tensor_rank_0() {
    let evaluation_times = zero_to_tau::<LENGTH>();
    let solution: TensorRank1List<2, 1, LENGTH> = ode1be(
        |t: &f64, y: &TensorRank1<2, 1>| TensorRank1::new([y[1], -t.sin()]),
        &evaluation_times,
        TensorRank1::new([0.0, 1.0]),
    )
    .expect("the unexpected");
    evaluation_times
        .iter()
        .zip(solution.iter())
        .for_each(|(t, y)| {
            assert!((t.sin() - y[0]).abs() < TOLERANCE || (t.sin() / y[0] - 1.0).abs() < TOLERANCE)
        });
}

#[test]
fn third_order_tensor_rank_0() {
    let evaluation_times = zero_to_tau::<LENGTH>();
    let solution: TensorRank1List<3, 1, LENGTH> = ode1be(
        |t: &f64, y: &TensorRank1<3, 1>| TensorRank1::new([y[1], y[2], -t.cos()]),
        &evaluation_times,
        TensorRank1::new([0.0, 1.0, 0.0]),
    )
    .expect("the unexpected");
    evaluation_times
        .iter()
        .zip(solution.iter())
        .for_each(|(t, y)| {
            assert!((t.sin() - y[0]).abs() < TOLERANCE || (t.sin() / y[0] - 1.0).abs() < TOLERANCE)
        });
}

#[test]
fn fourth_order_tensor_rank_0() {
    let evaluation_times = zero_to_tau::<LENGTH>();
    let solution: TensorRank1List<4, 1, LENGTH> = ode1be(
        |t: &f64, y: &TensorRank1<4, 1>| TensorRank1::new([y[1], y[2], y[3], t.sin()]),
        &evaluation_times,
        TensorRank1::new([0.0, 1.0, 0.0, -1.0]),
    )
    .expect("the unexpected");
    evaluation_times
        .iter()
        .zip(solution.iter())
        .for_each(|(t, y)| {
            assert!((t.sin() - y[0]).abs() < TOLERANCE || (t.sin() / y[0] - 1.0).abs() < TOLERANCE)
        });
}
