use super::{
    super::{
        super::{
            Tensor, TensorRank0, TensorRank0List, TensorRank1, TensorRank1List, TensorRank2,
            Tensors,
        },
        test::zero_to_tau,
    },
    Implicit, Ode1be,
};

const LENGTH: usize = 33;
const TOLERANCE: TensorRank0 = 1e-5;

#[test]
#[should_panic(expected = "Evaluation times must be strictly increasing.")]
fn evaluation_times_not_strictly_increasing() {
    let mut evaluation_times = zero_to_tau::<LENGTH>();
    evaluation_times[3] = evaluation_times[2];
    let _: TensorRank0List<LENGTH> = Ode1be {
        ..Default::default()
    }
    .integrate(
        |_: &TensorRank0, y: &TensorRank0| -y,
        |_: &TensorRank0, _: &TensorRank0| -1.0,
        0.0,
        1.0,
        &evaluation_times,
    )
    .expect("the unexpected");
}

#[test]
#[should_panic(expected = "Evaluation times precede the initial time.")]
fn evaluation_times_precede_initial_time() {
    let mut evaluation_times = zero_to_tau::<LENGTH>();
    evaluation_times[0] = -1.0;
    let _: TensorRank0List<LENGTH> = Ode1be {
        ..Default::default()
    }
    .integrate(
        |_: &TensorRank0, y: &TensorRank0| -y,
        |_: &TensorRank0, _: &TensorRank0| -1.0,
        0.0,
        1.0,
        &evaluation_times,
    )
    .expect("the unexpected");
}

#[test]
#[should_panic(expected = "Evaluation times must include a final time.")]
fn evaluation_times_no_final_time() {
    let _: TensorRank0List<LENGTH> = Ode1be {
        ..Default::default()
    }
    .integrate(
        |_: &TensorRank0, y: &TensorRank0| -y,
        |_: &TensorRank0, _: &TensorRank0| -1.0,
        0.0,
        1.0,
        &TensorRank0List::new([0.0]),
    )
    .expect("the unexpected");
}

#[test]
#[should_panic(expected = "asdfghjk")]
fn maximum_steps_reached() {
    let frequency = 1e5;
    let evaluation_times = zero_to_tau::<LENGTH>();
    let _: TensorRank0List<LENGTH> = Ode1be {
        max_steps: 3,
        ..Default::default()
    }
    .integrate(
        |t: &TensorRank0, _: &TensorRank0| (t * frequency).cos() * frequency,
        |_: &TensorRank0, _: &TensorRank0| 0.0,
        0.0,
        1.0,
        &evaluation_times,
    )
    .expect("the unexpected");
}

#[test]
fn first_order_tensor_rank_0() {
    let evaluation_times = zero_to_tau::<LENGTH>();
    let solution: TensorRank0List<LENGTH> = Ode1be {
        ..Default::default()
    }
    .integrate(
        |_: &TensorRank0, y: &TensorRank0| -y,
        |_: &TensorRank0, _: &TensorRank0| -1.0,
        0.0,
        1.0,
        &evaluation_times,
    )
    .expect("the unexpected");
    evaluation_times
        .iter()
        .zip(solution.iter())
        .for_each(|(t, y)| {
            assert!(((-t).exp() - y).abs() < TOLERANCE || ((-t).exp() / y - 1.0).abs() < TOLERANCE)
        });
}

#[test]
fn first_order_tensor_rank_0_one_evaluation_time_after_initial_time() {
    let evaluation_times = TensorRank0List::new([1.0]);
    let solution: TensorRank0List<LENGTH> = Ode1be {
        ..Default::default()
    }
    .integrate(
        |_: &TensorRank0, y: &TensorRank0| -y,
        |_: &TensorRank0, _: &TensorRank0| -1.0,
        0.0,
        1.0,
        &evaluation_times,
    )
    .expect("the unexpected");
    evaluation_times
        .iter()
        .zip(solution.iter())
        .for_each(|(t, y)| {
            assert!(((-t).exp() - y).abs() < TOLERANCE || ((-t).exp() / y - 1.0).abs() < TOLERANCE)
        });
}

#[test]
fn first_order_tensor_rank_0_first_evaluation_time() {
    let mut evaluation_times = zero_to_tau::<LENGTH>();
    evaluation_times[0] = 1e-8;
    evaluation_times[3] = evaluation_times[2] + 1e-8;
    evaluation_times[4] = evaluation_times[3] + 1e-8;
    evaluation_times[5] = evaluation_times[4] + 1e-8;
    let solution: TensorRank0List<LENGTH> = Ode1be {
        ..Default::default()
    }
    .integrate(
        |_: &TensorRank0, y: &TensorRank0| -y,
        |_: &TensorRank0, _: &TensorRank0| -1.0,
        0.0,
        1.0,
        &evaluation_times,
    )
    .expect("the unexpected");
    evaluation_times
        .iter()
        .zip(solution.iter())
        .for_each(|(t, y)| {
            assert!(((-t).exp() - y).abs() < TOLERANCE || ((-t).exp() / y - 1.0).abs() < TOLERANCE)
        });
}

#[test]
fn first_order_tensor_rank_0_nearby_evaluation_times() {
    let mut evaluation_times = zero_to_tau::<LENGTH>();
    evaluation_times[3] = evaluation_times[2] + 1e-10;
    evaluation_times[4] = evaluation_times[3] + 1e-10;
    evaluation_times[5] = evaluation_times[4] + 1e-10;
    let solution: TensorRank0List<LENGTH> = Ode1be {
        ..Default::default()
    }
    .integrate(
        |_: &TensorRank0, y: &TensorRank0| -y,
        |_: &TensorRank0, _: &TensorRank0| -1.0,
        0.0,
        1.0,
        &evaluation_times,
    )
    .expect("the unexpected");
    evaluation_times
        .iter()
        .zip(solution.iter())
        .for_each(|(t, y)| {
            assert!(((-t).exp() - y).abs() < TOLERANCE || ((-t).exp() / y - 1.0).abs() < TOLERANCE)
        });
}

#[test]
fn second_order_tensor_rank_0() {
    let evaluation_times = zero_to_tau::<LENGTH>();
    let solution: TensorRank1List<2, 1, LENGTH> = Ode1be {
        ..Default::default()
    }
    .integrate(
        |_: &TensorRank0, y: &TensorRank1<2, 1>| TensorRank1::new([y[1], -y[0]]),
        |_: &TensorRank0, _: &TensorRank1<2, 1>| TensorRank2::new([[0.0, -1.0], [1.0, 0.0]]),
        0.0,
        TensorRank1::new([0.0, 1.0]),
        &evaluation_times,
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
    let solution: TensorRank1List<3, 1, LENGTH> = Ode1be {
        ..Default::default()
    }
    .integrate(
        |_: &TensorRank0, y: &TensorRank1<3, 1>| TensorRank1::new([y[1], y[2], -y[1]]),
        |_: &TensorRank0, _: &TensorRank1<3, 1>| {
            TensorRank2::new([[0.0, 0.0, 0.0], [1.0, 0.0, -1.0], [0.0, 1.0, 0.0]])
        },
        0.0,
        TensorRank1::new([0.0, 1.0, 0.0]),
        &evaluation_times,
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
    let solution: TensorRank1List<4, 1, LENGTH> = Ode1be {
        ..Default::default()
    }
    .integrate(
        |_: &TensorRank0, y: &TensorRank1<4, 1>| TensorRank1::new([y[1], y[2], y[3], y[0]]),
        |_: &TensorRank0, _: &TensorRank1<4, 1>| {
            TensorRank2::new([
                [0.0, 0.0, 0.0, 1.0],
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
            ])
        },
        0.0,
        TensorRank1::new([0.0, 1.0, 0.0, -1.0]),
        &evaluation_times,
    )
    .expect("the unexpected");
    evaluation_times
        .iter()
        .zip(solution.iter())
        .for_each(|(t, y)| {
            assert!((t.sin() - y[0]).abs() < TOLERANCE || (t.sin() / y[0] - 1.0).abs() < TOLERANCE)
        });
}
