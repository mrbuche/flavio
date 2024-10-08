use super::{
    super::rank_1::list::TensorRank1ListTrait, list_2d::TensorRank2List2DTrait, TensorRank0,
    TensorRank1, TensorRank1List, TensorRank1Trait, TensorRank2, TensorRank2List2D,
    TensorRank2Trait,
};
use crate::test::assert_eq_within_tols;
use std::cmp::Ordering;

fn get_array_dim_2() -> [[TensorRank0; 2]; 2] {
    [[1.0, 2.0], [3.0, 4.0]]
}

fn get_array_dim_3() -> [[TensorRank0; 3]; 3] {
    [[1.0, 4.0, 6.0], [7.0, 2.0, 5.0], [9.0, 8.0, 3.0]]
}

fn get_array_dim_4() -> [[TensorRank0; 4]; 4] {
    [
        [1.0, 4.0, 6.0, 6.0],
        [1.0, 5.0, 1.0, 0.0],
        [1.0, 3.0, 5.0, 0.0],
        [1.0, 4.0, 6.0, 0.0],
    ]
}

fn get_array_dim_9() -> [[TensorRank0; 9]; 9] {
    [
        [2.0, 2.0, 4.0, 0.0, 0.0, 1.0, 1.0, 3.0, 3.0],
        [0.0, 3.0, 1.0, 0.0, 0.0, 1.0, 4.0, 2.0, 1.0],
        [3.0, 0.0, 1.0, 2.0, 0.0, 3.0, 4.0, 4.0, 2.0],
        [4.0, 4.0, 0.0, 2.0, 1.0, 1.0, 0.0, 0.0, 4.0],
        [0.0, 1.0, 0.0, 1.0, 1.0, 3.0, 0.0, 1.0, 1.0],
        [4.0, 2.0, 3.0, 4.0, 2.0, 4.0, 3.0, 0.0, 4.0],
        [1.0, 3.0, 2.0, 0.0, 0.0, 0.0, 2.0, 4.0, 2.0],
        [2.0, 2.0, 2.0, 4.0, 1.0, 2.0, 4.0, 2.0, 2.0],
        [1.0, 2.0, 3.0, 4.0, 0.0, 1.0, 4.0, 2.0, 1.0],
    ]
}

fn get_tensor_rank_1_a() -> TensorRank1<4, 1> {
    TensorRank1::new([1.0, 2.0, 3.0, 4.0])
}

fn get_tensor_rank_1_b() -> TensorRank1<4, 1> {
    TensorRank1::new([5.0, 7.0, 6.0, 8.0])
}

fn get_tensor_rank_2_dim_2() -> TensorRank2<2, 1, 1> {
    TensorRank2::new(get_array_dim_2())
}

fn get_tensor_rank_2_dim_3() -> TensorRank2<3, 1, 1> {
    TensorRank2::new(get_array_dim_3())
}

fn get_tensor_rank_2_dim_4() -> TensorRank2<4, 1, 1> {
    TensorRank2::new(get_array_dim_4())
}

fn get_tensor_rank_2_dim_9() -> TensorRank2<9, 1, 1> {
    TensorRank2::new(get_array_dim_9())
}

fn get_other_tensor_rank_2_dim_2() -> TensorRank2<2, 1, 1> {
    TensorRank2::new([[5.0, 6.0], [7.0, 8.0]])
}

fn get_other_tensor_rank_2_dim_3() -> TensorRank2<3, 1, 1> {
    TensorRank2::new([[3.0, 2.0, 3.0], [6.0, 5.0, 2.0], [4.0, 5.0, 0.0]])
}

fn get_other_tensor_rank_2_dim_4() -> TensorRank2<4, 1, 1> {
    TensorRank2::new([
        [3.0, 2.0, 3.0, 5.0],
        [6.0, 5.0, 2.0, 4.0],
        [4.0, 5.0, 0.0, 4.0],
        [4.0, 4.0, 1.0, 6.0],
    ])
}

fn get_diagonal_tensor_rank_2_dim_4() -> TensorRank2<4, 1, 1> {
    TensorRank2::new([
        [3.0, 0.0, 0.0, 0.0],
        [0.0, 5.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 6.0],
    ])
}

fn get_other_tensor_rank_2_dim_9() -> TensorRank2<9, 1, 1> {
    TensorRank2::new([
        [0.0, 4.0, 2.0, 0.0, 1.0, 4.0, 2.0, 4.0, 1.0],
        [1.0, 2.0, 2.0, 1.0, 0.0, 3.0, 0.0, 2.0, 0.0],
        [3.0, 0.0, 2.0, 3.0, 3.0, 0.0, 0.0, 0.0, 2.0],
        [2.0, 3.0, 0.0, 0.0, 1.0, 3.0, 3.0, 4.0, 2.0],
        [0.0, 4.0, 1.0, 3.0, 1.0, 1.0, 1.0, 2.0, 1.0],
        [1.0, 3.0, 0.0, 3.0, 3.0, 2.0, 1.0, 3.0, 4.0],
        [0.0, 0.0, 0.0, 1.0, 0.0, 3.0, 1.0, 3.0, 4.0],
        [2.0, 0.0, 4.0, 3.0, 1.0, 2.0, 0.0, 3.0, 4.0],
        [4.0, 2.0, 0.0, 0.0, 4.0, 0.0, 4.0, 2.0, 2.0],
    ])
}

fn get_other_tensor_rank_2_mul_tensor_rank_1_dim_4() -> TensorRank1<4, 1> {
    TensorRank1::new([51.0, 14.0, 22.0, 27.0])
}

fn get_other_tensor_rank_2_add_tensor_rank_2_dim_4() -> TensorRank2<4, 1, 1> {
    TensorRank2::new([
        [4.0, 6.0, 9.0, 11.0],
        [7.0, 10.0, 3.0, 4.0],
        [5.0, 8.0, 5.0, 4.0],
        [5.0, 8.0, 7.0, 6.0],
    ])
}

fn get_other_tensor_rank_2_sub_tensor_rank_2_dim_4() -> TensorRank2<4, 1, 1> {
    TensorRank2::new([
        [-2.0, 2.0, 3.0, 1.0],
        [-5.0, 0.0, -1.0, -4.0],
        [-3.0, -2.0, 5.0, -4.0],
        [-3.0, 0.0, 5.0, -6.0],
    ])
}

fn get_other_tensor_rank_2_mul_tensor_rank_2_dim_4() -> TensorRank2<4, 1, 1> {
    TensorRank2::new([
        [75.0, 76.0, 17.0, 81.0],
        [37.0, 32.0, 13.0, 29.0],
        [41.0, 42.0, 9.0, 37.0],
        [51.0, 52.0, 11.0, 45.0],
    ])
}

fn get_tensor_rank_1_list() -> TensorRank1List<3, 1, 8> {
    TensorRank1List::new([
        [5.0, 0.0, 0.0],
        [5.0, 5.0, 6.0],
        [3.0, 1.0, 4.0],
        [3.0, 4.0, 2.0],
        [1.0, 0.0, 3.0],
        [1.0, 3.0, 1.0],
        [1.0, 6.0, 0.0],
        [1.0, 1.0, 1.0],
    ])
}

fn get_tensor_rank_2_list_2d() -> TensorRank2List2D<3, 1, 1, 2, 2> {
    TensorRank2List2D::new([
        [
            [[1.0, 4.0, 6.0], [7.0, 2.0, 5.0], [9.0, 8.0, 3.0]],
            [[3.0, 2.0, 3.0], [6.0, 5.0, 2.0], [4.0, 5.0, 0.0]],
        ],
        [
            [[5.0, 2.0, 9.0], [2.0, 4.0, 5.0], [1.0, 3.0, 8.0]],
            [[4.0, 3.0, 2.0], [2.0, 5.0, 4.0], [1.0, 7.0, 1.0]],
        ],
    ])
}

fn get_tensor_rank_2_mul_tensor_rank_2_list_2d() -> TensorRank2List2D<3, 1, 1, 2, 2> {
    TensorRank2List2D::new([
        [
            [[83.0, 60.0, 44.0], [66.0, 72.0, 67.0], [92.0, 76.0, 103.0]],
            [[51.0, 52.0, 11.0], [53.0, 49.0, 25.0], [87.0, 73.0, 43.0]],
        ],
        [
            [[19.0, 36.0, 77.0], [44.0, 37.0, 113.0], [64.0, 59.0, 145.0]],
            [[18.0, 65.0, 24.0], [37.0, 66.0, 27.0], [55.0, 88.0, 53.0]],
        ],
    ])
}

#[test]
fn add_tensor_rank_2_to_self() {
    (get_tensor_rank_2_dim_4() + get_other_tensor_rank_2_dim_4())
        .iter()
        .zip(get_other_tensor_rank_2_add_tensor_rank_2_dim_4().iter())
        .for_each(|(tensor_rank_2_i, res_tensor_rank_2_i)| {
            tensor_rank_2_i
                .iter()
                .zip(res_tensor_rank_2_i.iter())
                .for_each(|(tensor_rank_2_ij, res_tensor_rank_2_ij)| {
                    assert_eq!(tensor_rank_2_ij, res_tensor_rank_2_ij)
                })
        });
}

#[test]
fn add_tensor_rank_2_ref_to_self() {
    (get_tensor_rank_2_dim_4() + &get_other_tensor_rank_2_dim_4())
        .iter()
        .zip(get_other_tensor_rank_2_add_tensor_rank_2_dim_4().iter())
        .for_each(|(tensor_rank_2_i, res_tensor_rank_2_i)| {
            tensor_rank_2_i
                .iter()
                .zip(res_tensor_rank_2_i.iter())
                .for_each(|(tensor_rank_2_ij, res_tensor_rank_2_ij)| {
                    assert_eq!(tensor_rank_2_ij, res_tensor_rank_2_ij)
                })
        });
}

#[test]
fn add_tensor_rank_2_to_self_ref() {
    (&get_tensor_rank_2_dim_4() + get_other_tensor_rank_2_dim_4())
        .iter()
        .zip(get_other_tensor_rank_2_add_tensor_rank_2_dim_4().iter())
        .for_each(|(tensor_rank_2_i, res_tensor_rank_2_i)| {
            tensor_rank_2_i
                .iter()
                .zip(res_tensor_rank_2_i.iter())
                .for_each(|(tensor_rank_2_ij, res_tensor_rank_2_ij)| {
                    assert_eq!(tensor_rank_2_ij, res_tensor_rank_2_ij)
                })
        });
}

#[test]
fn add_assign_tensor_rank_2() {
    let mut tensor_rank_2 = get_tensor_rank_2_dim_4();
    tensor_rank_2 += get_other_tensor_rank_2_dim_4();
    tensor_rank_2
        .iter()
        .zip(get_other_tensor_rank_2_add_tensor_rank_2_dim_4().iter())
        .for_each(|(tensor_rank_2_i, res_tensor_rank_2_i)| {
            tensor_rank_2_i
                .iter()
                .zip(res_tensor_rank_2_i.iter())
                .for_each(|(tensor_rank_2_ij, res_tensor_rank_2_ij)| {
                    assert_eq!(tensor_rank_2_ij, res_tensor_rank_2_ij)
                })
        });
}

#[test]
fn add_assign_tensor_rank_2_ref() {
    let mut tensor_rank_2 = get_tensor_rank_2_dim_4();
    tensor_rank_2 += &get_other_tensor_rank_2_dim_4();
    tensor_rank_2
        .iter()
        .zip(get_other_tensor_rank_2_add_tensor_rank_2_dim_4().iter())
        .for_each(|(tensor_rank_2_i, res_tensor_rank_2_i)| {
            tensor_rank_2_i
                .iter()
                .zip(res_tensor_rank_2_i.iter())
                .for_each(|(tensor_rank_2_ij, res_tensor_rank_2_ij)| {
                    assert_eq!(tensor_rank_2_ij, res_tensor_rank_2_ij)
                })
        });
}

#[test]
fn as_array_dim_2() {
    get_tensor_rank_2_dim_2()
        .as_array()
        .iter()
        .zip(get_array_dim_2().iter())
        .for_each(|(tensor_rank_2_as_array_i, array_i)| {
            tensor_rank_2_as_array_i
                .iter()
                .zip(array_i.iter())
                .for_each(|(tensor_rank_2_as_array_ij, array_ij)| {
                    assert_eq!(tensor_rank_2_as_array_ij, array_ij)
                })
        });
}

#[test]
fn as_array_dim_3() {
    get_tensor_rank_2_dim_3()
        .as_array()
        .iter()
        .zip(get_array_dim_3().iter())
        .for_each(|(tensor_rank_2_as_array_i, array_i)| {
            tensor_rank_2_as_array_i
                .iter()
                .zip(array_i.iter())
                .for_each(|(tensor_rank_2_as_array_ij, array_ij)| {
                    assert_eq!(tensor_rank_2_as_array_ij, array_ij)
                })
        });
}

#[test]
fn as_array_dim_4() {
    get_tensor_rank_2_dim_4()
        .as_array()
        .iter()
        .zip(get_array_dim_4().iter())
        .for_each(|(tensor_rank_2_as_array_i, array_i)| {
            tensor_rank_2_as_array_i
                .iter()
                .zip(array_i.iter())
                .for_each(|(tensor_rank_2_as_array_ij, array_ij)| {
                    assert_eq!(tensor_rank_2_as_array_ij, array_ij)
                })
        });
}

#[test]
fn as_array_dim_9() {
    get_tensor_rank_2_dim_9()
        .as_array()
        .iter()
        .zip(get_array_dim_9().iter())
        .for_each(|(tensor_rank_2_as_array_i, array_i)| {
            tensor_rank_2_as_array_i
                .iter()
                .zip(array_i.iter())
                .for_each(|(tensor_rank_2_as_array_ij, array_ij)| {
                    assert_eq!(tensor_rank_2_as_array_ij, array_ij)
                })
        });
}

#[test]
fn div_tensor_rank_0_to_self() {
    (get_tensor_rank_2_dim_4() / 3.3)
        .iter()
        .zip(get_array_dim_4().iter())
        .for_each(|(tensor_rank_2_i, array_i)| {
            tensor_rank_2_i
                .iter()
                .zip(array_i.iter())
                .for_each(|(tensor_rank_2_ij, array_ij)| {
                    assert_eq!(tensor_rank_2_ij, &(array_ij / 3.3))
                })
        });
}

#[test]
fn div_tensor_rank_0_to_self_ref() {
    (&get_tensor_rank_2_dim_4() / 3.3)
        .iter()
        .zip(get_array_dim_4().iter())
        .for_each(|(tensor_rank_2_i, array_i)| {
            tensor_rank_2_i
                .iter()
                .zip(array_i.iter())
                .for_each(|(tensor_rank_2_ij, array_ij)| {
                    assert_eq!(tensor_rank_2_ij, &(array_ij / 3.3))
                })
        });
}

#[test]
#[allow(clippy::op_ref)]
fn div_tensor_rank_0_ref_to_self() {
    (get_tensor_rank_2_dim_4() / &3.3)
        .iter()
        .zip(get_array_dim_4().iter())
        .for_each(|(tensor_rank_2_i, array_i)| {
            tensor_rank_2_i
                .iter()
                .zip(array_i.iter())
                .for_each(|(tensor_rank_2_ij, array_ij)| {
                    assert_eq!(tensor_rank_2_ij, &(array_ij / 3.3))
                })
        });
}

#[test]
#[allow(clippy::op_ref)]
fn div_tensor_rank_0_ref_to_self_ref() {
    (&get_tensor_rank_2_dim_4() / &3.3)
        .iter()
        .zip(get_array_dim_4().iter())
        .for_each(|(tensor_rank_2_i, array_i)| {
            tensor_rank_2_i
                .iter()
                .zip(array_i.iter())
                .for_each(|(tensor_rank_2_ij, array_ij)| {
                    assert_eq!(tensor_rank_2_ij, &(array_ij / 3.3))
                })
        });
}

#[test]
fn div_assign_tensor_rank_0() {
    let mut tensor_rank_2 = get_tensor_rank_2_dim_4();
    tensor_rank_2 /= 3.3;
    tensor_rank_2
        .iter()
        .zip(get_array_dim_4().iter())
        .for_each(|(tensor_rank_2_i, array_i)| {
            tensor_rank_2_i
                .iter()
                .zip(array_i.iter())
                .for_each(|(tensor_rank_2_ij, array_ij)| {
                    assert_eq!(tensor_rank_2_ij, &(array_ij / 3.3))
                })
        });
}

#[test]
fn div_assign_tensor_rank_0_ref() {
    let mut tensor_rank_2 = get_tensor_rank_2_dim_4();
    tensor_rank_2 /= &3.3;
    tensor_rank_2
        .iter()
        .zip(get_array_dim_4().iter())
        .for_each(|(tensor_rank_2_i, array_i)| {
            tensor_rank_2_i
                .iter()
                .zip(array_i.iter())
                .for_each(|(tensor_rank_2_ij, array_ij)| {
                    assert_eq!(tensor_rank_2_ij, &(array_ij / 3.3))
                })
        });
}

#[test]
fn determinant_dim_2() {
    assert_eq_within_tols(&get_tensor_rank_2_dim_2().determinant(), &(-2.0));
}

#[test]
fn determinant_dim_3() {
    assert_eq_within_tols(&get_tensor_rank_2_dim_3().determinant(), &290.0);
}

#[test]
fn determinant_dim_4() {
    assert_eq_within_tols(&get_tensor_rank_2_dim_4().determinant(), &36.0);
}

#[test]
fn determinant_dim_9() {
    assert_eq_within_tols(&get_tensor_rank_2_dim_9().determinant(), &2398.0);
}

#[test]
fn deviatoric_dim_2() {
    let tensor_rank_2 = get_tensor_rank_2_dim_2();
    let trace = tensor_rank_2.trace();
    let deviatoric_tensor_rank_2 = tensor_rank_2.deviatoric();
    assert_eq!(deviatoric_tensor_rank_2.trace(), 0.0);
    deviatoric_tensor_rank_2
        .iter()
        .enumerate()
        .zip(tensor_rank_2.iter())
        .for_each(|((i, deviatoric_tensor_rank_2_i), tensor_rank_2_i)| {
            deviatoric_tensor_rank_2_i
                .iter()
                .enumerate()
                .zip(tensor_rank_2_i.iter())
                .for_each(|((j, deviatoric_tensor_rank_2_ij), tensor_rank_2_ij)| {
                    assert_eq!(
                        deviatoric_tensor_rank_2_ij,
                        &(tensor_rank_2_ij - (((i == j) as u8) as TensorRank0) * trace / 2.0)
                    )
                })
        });
}

#[test]
fn deviatoric_dim_3() {
    let tensor_rank_2 = get_tensor_rank_2_dim_3();
    let trace = tensor_rank_2.trace();
    let deviatoric_tensor_rank_2 = tensor_rank_2.deviatoric();
    assert_eq!(deviatoric_tensor_rank_2.trace(), 0.0);
    deviatoric_tensor_rank_2
        .iter()
        .enumerate()
        .zip(tensor_rank_2.iter())
        .for_each(|((i, deviatoric_tensor_rank_2_i), tensor_rank_2_i)| {
            deviatoric_tensor_rank_2_i
                .iter()
                .enumerate()
                .zip(tensor_rank_2_i.iter())
                .for_each(|((j, deviatoric_tensor_rank_2_ij), tensor_rank_2_ij)| {
                    assert_eq!(
                        deviatoric_tensor_rank_2_ij,
                        &(tensor_rank_2_ij - (((i == j) as u8) as TensorRank0) * trace / 3.0)
                    )
                })
        });
}

#[test]
fn deviatoric_dim_4() {
    let tensor_rank_2 = get_tensor_rank_2_dim_4();
    let trace = tensor_rank_2.trace();
    let deviatoric_tensor_rank_2 = tensor_rank_2.deviatoric();
    assert_eq!(deviatoric_tensor_rank_2.trace(), 0.0);
    deviatoric_tensor_rank_2
        .iter()
        .enumerate()
        .zip(tensor_rank_2.iter())
        .for_each(|((i, deviatoric_tensor_rank_2_i), tensor_rank_2_i)| {
            deviatoric_tensor_rank_2_i
                .iter()
                .enumerate()
                .zip(tensor_rank_2_i.iter())
                .for_each(|((j, deviatoric_tensor_rank_2_ij), tensor_rank_2_ij)| {
                    assert_eq!(
                        deviatoric_tensor_rank_2_ij,
                        &(tensor_rank_2_ij - (((i == j) as u8) as TensorRank0) * trace / 4.0)
                    )
                })
        });
}

#[test]
fn deviatoric_dim_9() {
    let tensor_rank_2 = get_tensor_rank_2_dim_9();
    let trace = tensor_rank_2.trace();
    let deviatoric_tensor_rank_2 = tensor_rank_2.deviatoric();
    assert_eq!(deviatoric_tensor_rank_2.trace(), 0.0);
    deviatoric_tensor_rank_2
        .iter()
        .enumerate()
        .zip(tensor_rank_2.iter())
        .for_each(|((i, deviatoric_tensor_rank_2_i), tensor_rank_2_i)| {
            deviatoric_tensor_rank_2_i
                .iter()
                .enumerate()
                .zip(tensor_rank_2_i.iter())
                .for_each(|((j, deviatoric_tensor_rank_2_ij), tensor_rank_2_ij)| {
                    assert_eq!(
                        deviatoric_tensor_rank_2_ij,
                        &(tensor_rank_2_ij - (((i == j) as u8) as TensorRank0) * trace / 9.0)
                    )
                })
        });
}

#[test]
fn deviatoric_and_trace_dim_2() {
    let tensor_rank_2 = get_tensor_rank_2_dim_2();
    let (deviatoric, trace) = tensor_rank_2.deviatoric_and_trace();
    assert_eq!(trace, tensor_rank_2.trace());
    deviatoric
        .iter()
        .zip(tensor_rank_2.deviatoric().iter())
        .for_each(|(deviatoric_i, tensor_rank_2_deviatoric_i)| {
            deviatoric_i
                .iter()
                .zip(tensor_rank_2_deviatoric_i.iter())
                .for_each(|(deviatoric_ij, tensor_rank_2_deviatoric_ij)| {
                    assert_eq!(deviatoric_ij, tensor_rank_2_deviatoric_ij)
                })
        });
}

#[test]
fn deviatoric_and_trace_dim_3() {
    let tensor_rank_2 = get_tensor_rank_2_dim_3();
    let (deviatoric, trace) = tensor_rank_2.deviatoric_and_trace();
    assert_eq!(trace, tensor_rank_2.trace());
    deviatoric
        .iter()
        .zip(tensor_rank_2.deviatoric().iter())
        .for_each(|(deviatoric_i, tensor_rank_2_deviatoric_i)| {
            deviatoric_i
                .iter()
                .zip(tensor_rank_2_deviatoric_i.iter())
                .for_each(|(deviatoric_ij, tensor_rank_2_deviatoric_ij)| {
                    assert_eq!(deviatoric_ij, tensor_rank_2_deviatoric_ij)
                })
        });
}

#[test]
fn deviatoric_and_trace_dim_4() {
    let tensor_rank_2 = get_tensor_rank_2_dim_4();
    let (deviatoric, trace) = tensor_rank_2.deviatoric_and_trace();
    assert_eq!(trace, tensor_rank_2.trace());
    deviatoric
        .iter()
        .zip(tensor_rank_2.deviatoric().iter())
        .for_each(|(deviatoric_i, tensor_rank_2_deviatoric_i)| {
            deviatoric_i
                .iter()
                .zip(tensor_rank_2_deviatoric_i.iter())
                .for_each(|(deviatoric_ij, tensor_rank_2_deviatoric_ij)| {
                    assert_eq!(deviatoric_ij, tensor_rank_2_deviatoric_ij)
                })
        });
}

#[test]
fn deviatoric_and_trace_dim_9() {
    let tensor_rank_2 = get_tensor_rank_2_dim_9();
    let (deviatoric, trace) = tensor_rank_2.deviatoric_and_trace();
    assert_eq!(trace, tensor_rank_2.trace());
    deviatoric
        .iter()
        .zip(tensor_rank_2.deviatoric().iter())
        .for_each(|(deviatoric_i, tensor_rank_2_deviatoric_i)| {
            deviatoric_i
                .iter()
                .zip(tensor_rank_2_deviatoric_i.iter())
                .for_each(|(deviatoric_ij, tensor_rank_2_deviatoric_ij)| {
                    assert_eq!(deviatoric_ij, tensor_rank_2_deviatoric_ij)
                })
        });
}

#[test]
fn dyad() {
    let tensor_rank_1_a = get_tensor_rank_1_a();
    let tensor_rank_1_b = get_tensor_rank_1_b();
    let tensor_rank_2 = TensorRank2::dyad(&tensor_rank_1_a, &tensor_rank_1_b);
    tensor_rank_2.iter().zip(tensor_rank_1_a.iter()).for_each(
        |(tensor_rank_2_i, tensor_rank_1_a_i)| {
            tensor_rank_2_i.iter().zip(tensor_rank_1_b.iter()).for_each(
                |(tensor_rank_2_ij, tensor_rank_1_b_j)| {
                    assert_eq!(tensor_rank_2_ij, &(tensor_rank_1_a_i * tensor_rank_1_b_j))
                },
            )
        },
    );
}

#[test]
fn from_iter() {
    let into_iterator = get_tensor_rank_2_dim_4().0.into_iter();
    let tensor_rank_2 = TensorRank2::<4, 1, 1>::from_iter(get_tensor_rank_2_dim_4().0);
    tensor_rank_2
        .iter()
        .zip(into_iterator)
        .for_each(|(tensor_rank_2_i, value_i)| {
            tensor_rank_2_i
                .iter()
                .zip(value_i.iter())
                .for_each(|(tensor_rank_2_ij, value_ij)| assert_eq!(tensor_rank_2_ij, value_ij))
        });
}

#[test]
fn full_contraction_dim_2() {
    assert_eq_within_tols(
        &get_tensor_rank_2_dim_2().full_contraction(&get_other_tensor_rank_2_dim_2()),
        &70.0,
    );
}

#[test]
fn full_contraction_dim_3() {
    assert_eq_within_tols(
        &get_tensor_rank_2_dim_3().full_contraction(&get_other_tensor_rank_2_dim_3()),
        &167.0,
    );
}

#[test]
fn full_contraction_dim_4() {
    assert_eq_within_tols(
        &get_tensor_rank_2_dim_4().full_contraction(&get_other_tensor_rank_2_dim_4()),
        &137.0,
    );
}

#[test]
fn full_contraction_dim_9() {
    assert_eq_within_tols(
        &get_tensor_rank_2_dim_9().full_contraction(&get_other_tensor_rank_2_dim_9()),
        &269.0,
    );
}

#[test]
fn identity() {
    TensorRank2::<9, 1, 1>::identity()
        .iter()
        .enumerate()
        .for_each(|(i, tensor_rank_2_i)| {
            tensor_rank_2_i
                .iter()
                .enumerate()
                .for_each(|(j, tensor_rank_2_ij)| {
                    if i == j {
                        assert_eq!(tensor_rank_2_ij, &1.0)
                    } else {
                        assert_eq!(tensor_rank_2_ij, &0.0)
                    }
                })
        });
}

#[test]
fn inverse_dim_2() {
    (get_tensor_rank_2_dim_2() * get_tensor_rank_2_dim_2().inverse())
        .iter()
        .enumerate()
        .for_each(|(i, tensor_rank_2_i)| {
            tensor_rank_2_i
                .iter()
                .enumerate()
                .for_each(|(j, tensor_rank_2_ij)| {
                    if i == j {
                        assert_eq_within_tols(tensor_rank_2_ij, &1.0)
                    } else {
                        assert_eq_within_tols(tensor_rank_2_ij, &0.0)
                    }
                })
        });
}

#[test]
fn inverse_dim_3() {
    (get_tensor_rank_2_dim_3() * get_tensor_rank_2_dim_3().inverse())
        .iter()
        .enumerate()
        .for_each(|(i, tensor_rank_2_i)| {
            tensor_rank_2_i
                .iter()
                .enumerate()
                .for_each(|(j, tensor_rank_2_ij)| {
                    if i == j {
                        assert_eq_within_tols(tensor_rank_2_ij, &1.0)
                    } else {
                        assert_eq_within_tols(tensor_rank_2_ij, &0.0)
                    }
                })
        });
}

#[test]
fn inverse_dim_4() {
    (get_tensor_rank_2_dim_4() * get_tensor_rank_2_dim_4().inverse())
        .iter()
        .enumerate()
        .for_each(|(i, tensor_rank_2_i)| {
            tensor_rank_2_i
                .iter()
                .enumerate()
                .for_each(|(j, tensor_rank_2_ij)| {
                    if i == j {
                        assert_eq_within_tols(tensor_rank_2_ij, &1.0)
                    } else {
                        assert_eq_within_tols(tensor_rank_2_ij, &0.0)
                    }
                })
        });
}

#[test]
fn inverse_dim_9() {
    (get_tensor_rank_2_dim_9() * get_tensor_rank_2_dim_9().inverse())
        .iter()
        .enumerate()
        .for_each(|(i, tensor_rank_2_i)| {
            tensor_rank_2_i
                .iter()
                .enumerate()
                .for_each(|(j, tensor_rank_2_ij)| {
                    if i == j {
                        assert_eq_within_tols(tensor_rank_2_ij, &1.0)
                    } else {
                        assert_eq_within_tols(tensor_rank_2_ij, &0.0)
                    }
                })
        });
}

#[test]
fn inverse_and_determinant_dim_2() {
    let tensor_rank_2 = get_tensor_rank_2_dim_2();
    let (inverse, determinant) = tensor_rank_2.inverse_and_determinant();
    assert_eq!(determinant, tensor_rank_2.determinant());
    inverse.iter().zip(tensor_rank_2.inverse().iter()).for_each(
        |(inverse_i, tensor_rank_2_inverse_i)| {
            inverse_i
                .iter()
                .zip(tensor_rank_2_inverse_i.iter())
                .for_each(|(inverse_ij, tensor_rank_2_inverse_ij)| {
                    assert_eq!(inverse_ij, tensor_rank_2_inverse_ij)
                })
        },
    );
}

#[test]
fn inverse_and_determinant_dim_3() {
    let tensor_rank_2 = get_tensor_rank_2_dim_3();
    let (inverse, determinant) = tensor_rank_2.inverse_and_determinant();
    assert_eq!(determinant, tensor_rank_2.determinant());
    inverse.iter().zip(tensor_rank_2.inverse().iter()).for_each(
        |(inverse_i, tensor_rank_2_inverse_i)| {
            inverse_i
                .iter()
                .zip(tensor_rank_2_inverse_i.iter())
                .for_each(|(inverse_ij, tensor_rank_2_inverse_ij)| {
                    assert_eq!(inverse_ij, tensor_rank_2_inverse_ij)
                })
        },
    );
}

#[test]
fn inverse_and_determinant_dim_4() {
    let tensor_rank_2 = get_tensor_rank_2_dim_4();
    let (inverse, determinant) = tensor_rank_2.inverse_and_determinant();
    assert_eq!(determinant, tensor_rank_2.determinant());
    inverse.iter().zip(tensor_rank_2.inverse().iter()).for_each(
        |(inverse_i, tensor_rank_2_inverse_i)| {
            inverse_i
                .iter()
                .zip(tensor_rank_2_inverse_i.iter())
                .for_each(|(inverse_ij, tensor_rank_2_inverse_ij)| {
                    assert_eq!(inverse_ij, tensor_rank_2_inverse_ij)
                })
        },
    );
}

#[test]
fn inverse_transpose_dim_2() {
    (get_tensor_rank_2_dim_2().transpose() * get_tensor_rank_2_dim_2().inverse_transpose())
        .iter()
        .enumerate()
        .for_each(|(i, tensor_rank_2_i)| {
            tensor_rank_2_i
                .iter()
                .enumerate()
                .for_each(|(j, tensor_rank_2_ij)| {
                    if i == j {
                        assert_eq!(tensor_rank_2_ij, &1.0)
                    } else {
                        assert_eq!(tensor_rank_2_ij, &0.0)
                    }
                })
        });
}

#[test]
fn inverse_transpose_dim_3() {
    (get_tensor_rank_2_dim_3().transpose() * get_tensor_rank_2_dim_3().inverse_transpose())
        .iter()
        .enumerate()
        .for_each(|(i, tensor_rank_2_i)| {
            tensor_rank_2_i
                .iter()
                .enumerate()
                .for_each(|(j, tensor_rank_2_ij)| {
                    if i == j {
                        assert_eq_within_tols(tensor_rank_2_ij, &1.0)
                    } else {
                        assert_eq_within_tols(tensor_rank_2_ij, &0.0)
                    }
                })
        });
}

#[test]
fn inverse_transpose_dim_4() {
    (get_tensor_rank_2_dim_4().transpose() * get_tensor_rank_2_dim_4().inverse_transpose())
        .iter()
        .enumerate()
        .for_each(|(i, tensor_rank_2_i)| {
            tensor_rank_2_i
                .iter()
                .enumerate()
                .for_each(|(j, tensor_rank_2_ij)| {
                    if i == j {
                        assert_eq_within_tols(tensor_rank_2_ij, &1.0)
                    } else {
                        assert_eq_within_tols(tensor_rank_2_ij, &0.0)
                    }
                })
        });
}

#[test]
fn inverse_transpose_9() {
    (get_tensor_rank_2_dim_9().transpose() * get_tensor_rank_2_dim_9().inverse_transpose())
        .iter()
        .enumerate()
        .for_each(|(i, tensor_rank_2_i)| {
            tensor_rank_2_i
                .iter()
                .enumerate()
                .for_each(|(j, tensor_rank_2_ij)| {
                    if i == j {
                        assert_eq_within_tols(tensor_rank_2_ij, &1.0)
                    } else {
                        assert_eq_within_tols(tensor_rank_2_ij, &0.0)
                    }
                })
        });
}

#[test]
fn inverse_transpose_and_determinant_dim_2() {
    let tensor_rank_2 = get_tensor_rank_2_dim_2();
    let (inverse_transpose, determinant) = tensor_rank_2.inverse_transpose_and_determinant();
    assert_eq!(determinant, tensor_rank_2.determinant());
    inverse_transpose
        .iter()
        .zip(tensor_rank_2.inverse_transpose().iter())
        .for_each(|(inverse_transpose_i, tensor_rank_2_inverse_transpose_i)| {
            inverse_transpose_i
                .iter()
                .zip(tensor_rank_2_inverse_transpose_i.iter())
                .for_each(
                    |(inverse_transpose_ij, tensor_rank_2_inverse_transpose_ij)| {
                        assert_eq!(inverse_transpose_ij, tensor_rank_2_inverse_transpose_ij)
                    },
                )
        });
}

#[test]
fn inverse_transpose_and_determinant_dim_3() {
    let tensor_rank_2 = get_tensor_rank_2_dim_3();
    let (inverse_transpose, determinant) = tensor_rank_2.inverse_transpose_and_determinant();
    assert_eq!(determinant, tensor_rank_2.determinant());
    inverse_transpose
        .iter()
        .zip(tensor_rank_2.inverse_transpose().iter())
        .for_each(|(inverse_transpose_i, tensor_rank_2_inverse_transpose_i)| {
            inverse_transpose_i
                .iter()
                .zip(tensor_rank_2_inverse_transpose_i.iter())
                .for_each(
                    |(inverse_transpose_ij, tensor_rank_2_inverse_transpose_ij)| {
                        assert_eq!(inverse_transpose_ij, tensor_rank_2_inverse_transpose_ij)
                    },
                )
        });
}

#[test]
fn inverse_transpose_and_determinant_dim_4() {
    let tensor_rank_2 = get_tensor_rank_2_dim_4();
    let (inverse_transpose, determinant) = tensor_rank_2.inverse_transpose_and_determinant();
    assert_eq!(determinant, tensor_rank_2.determinant());
    inverse_transpose
        .iter()
        .zip(tensor_rank_2.inverse_transpose().iter())
        .for_each(|(inverse_transpose_i, tensor_rank_2_inverse_transpose_i)| {
            inverse_transpose_i
                .iter()
                .zip(tensor_rank_2_inverse_transpose_i.iter())
                .for_each(
                    |(inverse_transpose_ij, tensor_rank_2_inverse_transpose_ij)| {
                        assert_eq!(inverse_transpose_ij, tensor_rank_2_inverse_transpose_ij)
                    },
                )
        });
}

#[test]
fn is_diagonal() {
    assert!(get_diagonal_tensor_rank_2_dim_4().is_diagonal())
}

#[test]
fn is_not_diagonal() {
    assert!(!get_other_tensor_rank_2_dim_4().is_diagonal())
}

#[test]
fn is_diagonal_identity() {
    assert!(TensorRank2::<3, 0, 0>::identity().is_diagonal())
}

#[test]
fn is_diagonal_zero() {
    assert!(TensorRank2::<4, 1, 1>::zero().is_diagonal())
}

#[test]
fn is_identity_dim_3() {
    assert!(TensorRank2::<3, 0, 0>::identity().is_identity())
}

#[test]
fn is_not_identity_dim_3() {
    assert!(!TensorRank2::<3, 0, 0>::zero().is_identity())
}

#[test]
fn is_identity_dim_4() {
    assert!(TensorRank2::<4, 1, 1>::identity().is_identity())
}

#[test]
fn is_not_identity_dim_4() {
    assert!(!get_diagonal_tensor_rank_2_dim_4().is_identity())
}

#[test]
fn is_zero_dim_3() {
    assert!(TensorRank2::<3, 0, 0>::zero().is_zero())
}

#[test]
fn is_not_zero_dim_3() {
    assert!(!TensorRank2::<3, 0, 0>::identity().is_zero())
}

#[test]
fn is_zero_dim_4() {
    assert!(TensorRank2::<4, 1, 1>::zero().is_zero())
}

#[test]
fn is_not_zero_dim_4() {
    assert!(!get_other_tensor_rank_2_dim_4().is_zero())
}

#[test]
fn iter() {
    get_tensor_rank_2_dim_4()
        .iter()
        .zip(get_array_dim_4().iter())
        .for_each(|(tensor_rank_2_i, array_i)| {
            tensor_rank_2_i
                .iter()
                .zip(array_i.iter())
                .for_each(|(tensor_rank_2_ij, array_ij)| assert_eq!(tensor_rank_2_ij, array_ij))
        });
}

#[test]
fn iter_mut() {
    get_tensor_rank_2_dim_4()
        .iter_mut()
        .zip(get_array_dim_4().iter_mut())
        .for_each(|(tensor_rank_2_i, array_i)| {
            tensor_rank_2_i
                .iter_mut()
                .zip(array_i.iter_mut())
                .for_each(|(tensor_rank_2_ij, array_ij)| assert_eq!(tensor_rank_2_ij, array_ij))
        });
}

#[test]
fn lu_decomposition() {
    let (tensor_l, tensor_u) = get_tensor_rank_2_dim_9().lu_decomposition();
    tensor_l
        .iter()
        .enumerate()
        .zip(tensor_u.iter())
        .for_each(|((i, tensor_l_i), tensor_u_i)| {
            tensor_l_i
                .iter()
                .enumerate()
                .zip(tensor_u_i.iter())
                .for_each(|((j, tensor_l_ij), tensor_u_ij)| match i.cmp(&j) {
                    Ordering::Greater => assert_eq!(tensor_u_ij, &0.0),
                    Ordering::Less => assert_eq!(tensor_l_ij, &0.0),
                    _ => (),
                })
        });
}

#[test]
fn mul_tensor_rank_0_to_self() {
    (get_tensor_rank_2_dim_4() * 3.3)
        .iter()
        .zip(get_array_dim_4().iter())
        .for_each(|(tensor_rank_2_i, array_i)| {
            tensor_rank_2_i
                .iter()
                .zip(array_i.iter())
                .for_each(|(tensor_rank_2_ij, array_ij)| {
                    assert_eq!(tensor_rank_2_ij, &(array_ij * 3.3))
                })
        });
}

#[test]
fn mul_tensor_rank_0_to_self_ref() {
    (&get_tensor_rank_2_dim_4() * 3.3)
        .iter()
        .zip(get_array_dim_4().iter())
        .for_each(|(tensor_rank_2_i, array_i)| {
            tensor_rank_2_i
                .iter()
                .zip(array_i.iter())
                .for_each(|(tensor_rank_2_ij, array_ij)| {
                    assert_eq!(tensor_rank_2_ij, &(array_ij * 3.3))
                })
        });
}

#[test]
#[allow(clippy::op_ref)]
fn mul_tensor_rank_0_ref_to_self() {
    (get_tensor_rank_2_dim_4() * &3.3)
        .iter()
        .zip(get_array_dim_4().iter())
        .for_each(|(tensor_rank_2_i, array_i)| {
            tensor_rank_2_i
                .iter()
                .zip(array_i.iter())
                .for_each(|(tensor_rank_2_ij, array_ij)| {
                    assert_eq!(tensor_rank_2_ij, &(array_ij * 3.3))
                })
        });
}

#[test]
#[allow(clippy::op_ref)]
fn mul_tensor_rank_0_ref_to_self_ref() {
    (&get_tensor_rank_2_dim_4() * &3.3)
        .iter()
        .zip(get_array_dim_4().iter())
        .for_each(|(tensor_rank_2_i, array_i)| {
            tensor_rank_2_i
                .iter()
                .zip(array_i.iter())
                .for_each(|(tensor_rank_2_ij, array_ij)| {
                    assert_eq!(tensor_rank_2_ij, &(array_ij * 3.3))
                })
        });
}

#[test]
fn mul_assign_tensor_rank_0() {
    let mut tensor_rank_2 = get_tensor_rank_2_dim_4();
    tensor_rank_2 *= 3.3;
    tensor_rank_2
        .iter()
        .zip(get_array_dim_4().iter())
        .for_each(|(tensor_rank_2_i, array_i)| {
            tensor_rank_2_i
                .iter()
                .zip(array_i.iter())
                .for_each(|(tensor_rank_2_ij, array_ij)| {
                    assert_eq!(tensor_rank_2_ij, &(array_ij * 3.3))
                })
        });
}

#[test]
fn mul_assign_tensor_rank_0_ref() {
    let mut tensor_rank_2 = get_tensor_rank_2_dim_4();
    tensor_rank_2 *= &3.3;
    tensor_rank_2
        .iter()
        .zip(get_array_dim_4().iter())
        .for_each(|(tensor_rank_2_i, array_i)| {
            tensor_rank_2_i
                .iter()
                .zip(array_i.iter())
                .for_each(|(tensor_rank_2_ij, array_ij)| {
                    assert_eq!(tensor_rank_2_ij, &(array_ij * 3.3))
                })
        });
}

#[test]
fn mul_tensor_rank_1_to_self() {
    (get_tensor_rank_2_dim_4() * get_tensor_rank_1_a())
        .iter()
        .zip(get_other_tensor_rank_2_mul_tensor_rank_1_dim_4().iter())
        .for_each(|(tensor_rank_1_i, res_tensor_rank_1_i)| {
            assert_eq!(tensor_rank_1_i, res_tensor_rank_1_i)
        });
}

#[test]
fn mul_tensor_rank_1_ref_to_self() {
    (get_tensor_rank_2_dim_4() * &get_tensor_rank_1_a())
        .iter()
        .zip(get_other_tensor_rank_2_mul_tensor_rank_1_dim_4().iter())
        .for_each(|(tensor_rank_1_i, res_tensor_rank_1_i)| {
            assert_eq!(tensor_rank_1_i, res_tensor_rank_1_i)
        });
}

#[test]
fn mul_tensor_rank_1_to_self_ref() {
    (&get_tensor_rank_2_dim_4() * get_tensor_rank_1_a())
        .iter()
        .zip(get_other_tensor_rank_2_mul_tensor_rank_1_dim_4().iter())
        .for_each(|(tensor_rank_1_i, res_tensor_rank_1_i)| {
            assert_eq!(tensor_rank_1_i, res_tensor_rank_1_i)
        });
}

#[test]
fn mul_tensor_rank_1_ref_to_self_ref() {
    (&get_tensor_rank_2_dim_4() * &get_tensor_rank_1_a())
        .iter()
        .zip(get_other_tensor_rank_2_mul_tensor_rank_1_dim_4().iter())
        .for_each(|(tensor_rank_1_i, res_tensor_rank_1_i)| {
            assert_eq!(tensor_rank_1_i, res_tensor_rank_1_i)
        });
}

#[test]
fn mul_tensor_rank_2_to_self() {
    (get_tensor_rank_2_dim_4() * get_other_tensor_rank_2_dim_4())
        .iter()
        .zip(get_other_tensor_rank_2_mul_tensor_rank_2_dim_4().iter())
        .for_each(|(tensor_rank_2_i, res_tensor_rank_2_i)| {
            tensor_rank_2_i
                .iter()
                .zip(res_tensor_rank_2_i.iter())
                .for_each(|(tensor_rank_2_ij, res_tensor_rank_2_ij)| {
                    assert_eq!(tensor_rank_2_ij, res_tensor_rank_2_ij)
                })
        });
}

#[test]
fn mul_tensor_rank_2_ref_to_self() {
    (get_tensor_rank_2_dim_4() * &get_other_tensor_rank_2_dim_4())
        .iter()
        .zip(get_other_tensor_rank_2_mul_tensor_rank_2_dim_4().iter())
        .for_each(|(tensor_rank_2_i, res_tensor_rank_2_i)| {
            tensor_rank_2_i
                .iter()
                .zip(res_tensor_rank_2_i.iter())
                .for_each(|(tensor_rank_2_ij, res_tensor_rank_2_ij)| {
                    assert_eq!(tensor_rank_2_ij, res_tensor_rank_2_ij)
                })
        });
}

#[test]
fn mul_tensor_rank_2_to_self_ref() {
    (&get_tensor_rank_2_dim_4() * get_other_tensor_rank_2_dim_4())
        .iter()
        .zip(get_other_tensor_rank_2_mul_tensor_rank_2_dim_4().iter())
        .for_each(|(tensor_rank_2_i, res_tensor_rank_2_i)| {
            tensor_rank_2_i
                .iter()
                .zip(res_tensor_rank_2_i.iter())
                .for_each(|(tensor_rank_2_ij, res_tensor_rank_2_ij)| {
                    assert_eq!(tensor_rank_2_ij, res_tensor_rank_2_ij)
                })
        });
}

#[test]
fn mul_tensor_rank_2_ref_to_self_ref() {
    (&get_tensor_rank_2_dim_4() * &get_other_tensor_rank_2_dim_4())
        .iter()
        .zip(get_other_tensor_rank_2_mul_tensor_rank_2_dim_4().iter())
        .for_each(|(tensor_rank_2_i, res_tensor_rank_2_i)| {
            tensor_rank_2_i
                .iter()
                .zip(res_tensor_rank_2_i.iter())
                .for_each(|(tensor_rank_2_ij, res_tensor_rank_2_ij)| {
                    assert_eq!(tensor_rank_2_ij, res_tensor_rank_2_ij)
                })
        });
}

#[test]
fn mul_tensor_rank_1_list_to_self() {
    (get_tensor_rank_2_dim_3() * get_tensor_rank_1_list())
        .iter()
        .zip(get_tensor_rank_1_list().iter())
        .for_each(|(res_tensor_rank_1, tensor_rank_1)| {
            res_tensor_rank_1
                .iter()
                .zip((get_tensor_rank_2_dim_3() * tensor_rank_1).iter())
                .for_each(|(res_tensor_rank_1_i, tensor_rank_1_i)| {
                    assert_eq!(res_tensor_rank_1_i, tensor_rank_1_i)
                })
        })
}

#[test]
fn mul_tensor_rank_1_list_ref_to_self() {
    (get_tensor_rank_2_dim_3() * &get_tensor_rank_1_list())
        .iter()
        .zip(get_tensor_rank_1_list().iter())
        .for_each(|(res_tensor_rank_1, tensor_rank_1)| {
            res_tensor_rank_1
                .iter()
                .zip((get_tensor_rank_2_dim_3() * tensor_rank_1).iter())
                .for_each(|(res_tensor_rank_1_i, tensor_rank_1_i)| {
                    assert_eq!(res_tensor_rank_1_i, tensor_rank_1_i)
                })
        })
}

#[test]
fn mul_tensor_rank_1_list_to_self_ref() {
    (&get_tensor_rank_2_dim_3() * get_tensor_rank_1_list())
        .iter()
        .zip(get_tensor_rank_1_list().iter())
        .for_each(|(res_tensor_rank_1, tensor_rank_1)| {
            res_tensor_rank_1
                .iter()
                .zip((get_tensor_rank_2_dim_3() * tensor_rank_1).iter())
                .for_each(|(res_tensor_rank_1_i, tensor_rank_1_i)| {
                    assert_eq!(res_tensor_rank_1_i, tensor_rank_1_i)
                })
        })
}

#[test]
fn mul_tensor_rank_1_list_ref_to_self_ref() {
    (&get_tensor_rank_2_dim_3() * &get_tensor_rank_1_list())
        .iter()
        .zip(get_tensor_rank_1_list().iter())
        .for_each(|(res_tensor_rank_1, tensor_rank_1)| {
            res_tensor_rank_1
                .iter()
                .zip((get_tensor_rank_2_dim_3() * tensor_rank_1).iter())
                .for_each(|(res_tensor_rank_1_i, tensor_rank_1_i)| {
                    assert_eq!(res_tensor_rank_1_i, tensor_rank_1_i)
                })
        })
}

#[test]
fn mul_tensor_rank_2_list_2d_to_self() {
    (get_tensor_rank_2_dim_3() * get_tensor_rank_2_list_2d())
        .iter()
        .zip(get_tensor_rank_2_mul_tensor_rank_2_list_2d().iter())
        .for_each(|(tensor_rank_2_list_2d_entry, res_entry)| {
            tensor_rank_2_list_2d_entry
                .iter()
                .zip(res_entry.iter())
                .for_each(|(tensor_rank_2, res)| {
                    tensor_rank_2
                        .iter()
                        .zip(res.iter())
                        .for_each(|(tensor_rank_2_i, res_i)| {
                            tensor_rank_2_i.iter().zip(res_i.iter()).for_each(
                                |(tensor_rank_2_ij, res_ij)| assert_eq!(tensor_rank_2_ij, res_ij),
                            )
                        })
                })
        });
}

#[test]
fn mul_tensor_rank_2_list_2d_to_self_ref() {
    (&get_tensor_rank_2_dim_3() * get_tensor_rank_2_list_2d())
        .iter()
        .zip(get_tensor_rank_2_mul_tensor_rank_2_list_2d().iter())
        .for_each(|(tensor_rank_2_list_2d_entry, res_entry)| {
            tensor_rank_2_list_2d_entry
                .iter()
                .zip(res_entry.iter())
                .for_each(|(tensor_rank_2, res)| {
                    tensor_rank_2
                        .iter()
                        .zip(res.iter())
                        .for_each(|(tensor_rank_2_i, res_i)| {
                            tensor_rank_2_i.iter().zip(res_i.iter()).for_each(
                                |(tensor_rank_2_ij, res_ij)| assert_eq!(tensor_rank_2_ij, res_ij),
                            )
                        })
                })
        });
}

#[test]
fn new() {
    get_tensor_rank_2_dim_4()
        .iter()
        .zip(get_array_dim_4().iter())
        .for_each(|(tensor_rank_2_i, array_i)| {
            tensor_rank_2_i
                .iter()
                .zip(array_i.iter())
                .for_each(|(tensor_rank_2_ij, array_ij)| assert_eq!(tensor_rank_2_ij, array_ij))
        });
}

#[test]
fn norm_dim_2() {
    assert_eq!(get_tensor_rank_2_dim_2().norm(), 5.477_225_575_051_661);
}

#[test]
fn norm_dim_3() {
    assert_eq!(get_tensor_rank_2_dim_3().norm(), 16.881_943_016_134_134);
}

#[test]
fn norm_dim_4() {
    assert_eq!(get_tensor_rank_2_dim_4().norm(), 14.282_856_857_085_7);
}

#[test]
fn norm_dim_9() {
    assert_eq!(get_tensor_rank_2_dim_9().norm(), 20.976_176_963_403_03);
}

#[test]
fn size() {
    assert_eq!(
        std::mem::size_of::<TensorRank2::<3, 1, 1>>(),
        std::mem::size_of::<[TensorRank1::<3, 1>; 3]>()
    )
}

#[test]
fn second_invariant() {
    assert_eq!(get_tensor_rank_2_dim_4().second_invariant(), 16.0);
}

#[test]
fn squared_trace_dim_2() {
    assert_eq_within_tols(&get_tensor_rank_2_dim_2().squared_trace(), &29.0);
}

#[test]
fn squared_trace_dim_3() {
    assert_eq_within_tols(&get_tensor_rank_2_dim_3().squared_trace(), &258.0);
}

#[test]
fn squared_trace_dim_4() {
    assert_eq_within_tols(&get_tensor_rank_2_dim_4().squared_trace(), &89.0);
}

#[test]
fn squared_trace_dim_9() {
    assert_eq_within_tols(&get_tensor_rank_2_dim_9().squared_trace(), &318.0);
}

#[test]
fn sub_tensor_rank_2_to_self() {
    (get_tensor_rank_2_dim_4() - get_other_tensor_rank_2_dim_4())
        .iter()
        .zip(get_other_tensor_rank_2_sub_tensor_rank_2_dim_4().iter())
        .for_each(|(tensor_rank_2_i, res_tensor_rank_2_i)| {
            tensor_rank_2_i
                .iter()
                .zip(res_tensor_rank_2_i.iter())
                .for_each(|(tensor_rank_2_ij, res_tensor_rank_2_ij)| {
                    assert_eq!(tensor_rank_2_ij, res_tensor_rank_2_ij)
                })
        });
}

#[test]
fn sub_tensor_rank_2_ref_to_self() {
    (get_tensor_rank_2_dim_4() - &get_other_tensor_rank_2_dim_4())
        .iter()
        .zip(get_other_tensor_rank_2_sub_tensor_rank_2_dim_4().iter())
        .for_each(|(tensor_rank_2_i, res_tensor_rank_2_i)| {
            tensor_rank_2_i
                .iter()
                .zip(res_tensor_rank_2_i.iter())
                .for_each(|(tensor_rank_2_ij, res_tensor_rank_2_ij)| {
                    assert_eq!(tensor_rank_2_ij, res_tensor_rank_2_ij)
                })
        });
}

#[test]
fn sub_assign_tensor_rank_2() {
    let mut tensor_rank_2 = get_tensor_rank_2_dim_4();
    tensor_rank_2 -= get_other_tensor_rank_2_dim_4();
    tensor_rank_2
        .iter()
        .zip(get_other_tensor_rank_2_sub_tensor_rank_2_dim_4().iter())
        .for_each(|(tensor_rank_2_i, res_tensor_rank_2_i)| {
            tensor_rank_2_i
                .iter()
                .zip(res_tensor_rank_2_i.iter())
                .for_each(|(tensor_rank_2_ij, res_tensor_rank_2_ij)| {
                    assert_eq!(tensor_rank_2_ij, res_tensor_rank_2_ij)
                })
        });
}

#[test]
fn sub_assign_tensor_rank_2_ref() {
    let mut tensor_rank_2 = get_tensor_rank_2_dim_4();
    tensor_rank_2 -= &get_other_tensor_rank_2_dim_4();
    tensor_rank_2
        .iter()
        .zip(get_other_tensor_rank_2_sub_tensor_rank_2_dim_4().iter())
        .for_each(|(tensor_rank_2_i, res_tensor_rank_2_i)| {
            tensor_rank_2_i
                .iter()
                .zip(res_tensor_rank_2_i.iter())
                .for_each(|(tensor_rank_2_ij, res_tensor_rank_2_ij)| {
                    assert_eq!(tensor_rank_2_ij, res_tensor_rank_2_ij)
                })
        });
}

#[test]
fn trace_dim_2() {
    assert_eq!(get_tensor_rank_2_dim_2().trace(), 5.0);
}

#[test]
fn trace_dim_3() {
    assert_eq!(get_tensor_rank_2_dim_3().trace(), 6.0);
}

#[test]
fn trace_dim_4() {
    assert_eq!(get_tensor_rank_2_dim_4().trace(), 11.0);
}

#[test]
fn trace_dim_9() {
    assert_eq!(get_tensor_rank_2_dim_9().trace(), 18.0);
}

#[test]
fn transpose() {
    let tensor_rank_2 = get_tensor_rank_2_dim_4();
    let tensor_rank_2_transpose = tensor_rank_2.transpose();
    tensor_rank_2
        .iter()
        .enumerate()
        .for_each(|(i, tensor_rank_2_i)| {
            tensor_rank_2_i
                .iter()
                .enumerate()
                .for_each(|(j, tensor_rank_2_ij)| {
                    assert_eq!(tensor_rank_2_ij, &tensor_rank_2_transpose[j][i])
                })
        });
    tensor_rank_2_transpose
        .iter()
        .enumerate()
        .for_each(|(i, tensor_rank_2_transpose_i)| {
            tensor_rank_2_transpose_i.iter().enumerate().for_each(
                |(j, tensor_rank_2_transpose_ij)| {
                    assert_eq!(tensor_rank_2_transpose_ij, &tensor_rank_2[j][i])
                },
            )
        });
}

#[test]
fn zero_dim_2() {
    TensorRank2::<2, 1, 1>::zero()
        .iter()
        .for_each(|tensor_rank_2_i| {
            tensor_rank_2_i
                .iter()
                .for_each(|tensor_rank_2_ij| assert_eq!(tensor_rank_2_ij, &0.0))
        });
}

#[test]
fn zero_dim_3() {
    TensorRank2::<3, 1, 1>::zero()
        .iter()
        .for_each(|tensor_rank_2_i| {
            tensor_rank_2_i
                .iter()
                .for_each(|tensor_rank_2_ij| assert_eq!(tensor_rank_2_ij, &0.0))
        });
}

#[test]
fn zero_dim_4() {
    TensorRank2::<4, 1, 1>::zero()
        .iter()
        .for_each(|tensor_rank_2_i| {
            tensor_rank_2_i
                .iter()
                .for_each(|tensor_rank_2_ij| assert_eq!(tensor_rank_2_ij, &0.0))
        });
}

#[test]
fn zero_dim_9() {
    TensorRank2::<9, 1, 1>::zero()
        .iter()
        .for_each(|tensor_rank_2_i| {
            tensor_rank_2_i
                .iter()
                .for_each(|tensor_rank_2_ij| assert_eq!(tensor_rank_2_ij, &0.0))
        });
}
