use super::*;
use crate::math::{TensorRank1, TensorRank1List, TensorRank1Trait};

fn ode_23_cos_fun(t: &TensorRank0, y: &TensorRank1<2, 1>) -> TensorRank1<2, 1> {
    TensorRank1::new([y[1], -t.sin()])
}

#[test]
fn ode_23_cos() {
    let evaluation_times = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
    let _: TensorRank1List<2, 1, 2> = ode_23(
        |t: &f64, y: &TensorRank1<2, 1>| ode_23_cos_fun(t, y),
        evaluation_times,
        TensorRank1::new([0.0, 1.0]),
    );
}
