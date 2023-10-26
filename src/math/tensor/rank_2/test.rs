use crate::test::assert_eq_within_tols;
use super::
{
    TensorRank0,
    TensorRank1,
    TensorRank1Trait,
    TensorRank2,
    TensorRank2Trait
};

fn get_array_dim_2() -> [[TensorRank0; 2]; 2]
{
    [[1.0, 2.0],
     [3.0, 4.0]]
}

fn get_array_dim_3() -> [[TensorRank0; 3]; 3]
{
    [[1.0, 4.0, 6.0],
     [7.0, 2.0, 5.0],
     [9.0, 8.0, 3.0]]
}

fn get_array_dim_4() -> [[TensorRank0; 4]; 4]
{
    [[1.0, 4.0, 6.0, 6.0],
     [1.0, 5.0, 1.0, 0.0],
     [1.0, 3.0, 5.0, 0.0],
     [1.0, 4.0, 6.0, 0.0]]
}

fn get_array_dim_9() -> [[TensorRank0; 9]; 9]
{
    [[2.0, 2.0, 4.0, 0.0, 0.0, 1.0, 1.0, 3.0, 3.0],
     [0.0, 3.0, 1.0, 0.0, 0.0, 1.0, 4.0, 2.0, 1.0],
     [3.0, 0.0, 1.0, 2.0, 0.0, 3.0, 4.0, 4.0, 2.0],
     [4.0, 4.0, 0.0, 2.0, 1.0, 1.0, 0.0, 0.0, 4.0],
     [0.0, 1.0, 0.0, 1.0, 1.0, 3.0, 0.0, 1.0, 1.0],
     [4.0, 2.0, 3.0, 4.0, 2.0, 4.0, 3.0, 0.0, 4.0],
     [1.0, 3.0, 2.0, 0.0, 0.0, 0.0, 2.0, 4.0, 2.0],
     [2.0, 2.0, 2.0, 4.0, 1.0, 2.0, 4.0, 2.0, 2.0],
     [1.0, 2.0, 3.0, 4.0, 0.0, 1.0, 4.0, 2.0, 1.0]]
}

fn get_tensor_rank_1_a() -> TensorRank1<4, 1>
{
    TensorRank1::new([1.0, 2.0, 3.0, 4.0])
}

fn get_tensor_rank_1_b() -> TensorRank1<4, 1>
{
    TensorRank1::new([5.0, 7.0, 6.0, 8.0])
}

fn get_tensor_rank_2_dim_2() -> TensorRank2<2, 1, 1>
{
    TensorRank2::new(get_array_dim_2())
}

fn get_tensor_rank_2_dim_3() -> TensorRank2<3, 1, 1>
{
    TensorRank2::new(get_array_dim_3())
}

fn get_tensor_rank_2_dim_4() -> TensorRank2<4, 1, 1>
{
    TensorRank2::new(get_array_dim_4())
}

fn get_tensor_rank_2_dim_9() -> TensorRank2<9, 1, 1>
{
    TensorRank2::new(get_array_dim_9())
}

fn get_other_tensor_rank_2_dim_2() -> TensorRank2<2, 1, 1>
{
    TensorRank2::new([
        [5.0, 6.0],
        [7.0, 8.0]
    ])
}

fn get_other_tensor_rank_2_dim_3() -> TensorRank2<3, 1, 1>
{
    TensorRank2::new([
        [3.0, 2.0, 3.0],
        [6.0, 5.0, 2.0],
        [4.0, 5.0, 0.0]
    ])
}

fn get_other_tensor_rank_2_dim_4() -> TensorRank2<4, 1, 1>
{
    TensorRank2::new([
        [3.0, 2.0, 3.0, 5.0],
        [6.0, 5.0, 2.0, 4.0],
        [4.0, 5.0, 0.0, 4.0],
        [4.0, 4.0, 1.0, 6.0]
    ])
}

fn get_other_tensor_rank_2_dim_9() -> TensorRank2<9, 1, 1>
{
    TensorRank2::new([
        [0.0, 4.0, 2.0, 0.0, 1.0, 4.0, 2.0, 4.0, 1.0],
        [1.0, 2.0, 2.0, 1.0, 0.0, 3.0, 0.0, 2.0, 0.0],
        [3.0, 0.0, 2.0, 3.0, 3.0, 0.0, 0.0, 0.0, 2.0],
        [2.0, 3.0, 0.0, 0.0, 1.0, 3.0, 3.0, 4.0, 2.0],
        [0.0, 4.0, 1.0, 3.0, 1.0, 1.0, 1.0, 2.0, 1.0],
        [1.0, 3.0, 0.0, 3.0, 3.0, 2.0, 1.0, 3.0, 4.0],
        [0.0, 0.0, 0.0, 1.0, 0.0, 3.0, 1.0, 3.0, 4.0],
        [2.0, 0.0, 4.0, 3.0, 1.0, 2.0, 0.0, 3.0, 4.0],
        [4.0, 2.0, 0.0, 0.0, 4.0, 0.0, 4.0, 2.0, 2.0]
    ])
}

fn get_other_tensor_rank_2_mul_tensor_rank_1_dim_4() -> TensorRank1<4, 1>
{
    TensorRank1::new([51.0, 14.0, 22.0, 27.0])
}

fn get_other_tensor_rank_2_add_tensor_rank_2_dim_4() -> TensorRank2<4, 1, 1>
{
    TensorRank2::new([
        [4.0,  6.0, 9.0, 11.0],
        [7.0, 10.0, 3.0,  4.0],
        [5.0,  8.0, 5.0,  4.0],
        [5.0,  8.0, 7.0,  6.0]
    ])
}

fn get_other_tensor_rank_2_sub_tensor_rank_2_dim_4() -> TensorRank2<4, 1, 1>
{
    TensorRank2::new([
        [-2.0,  2.0,  3.0,  1.0],
        [-5.0,  0.0, -1.0, -4.0],
        [-3.0, -2.0,  5.0, -4.0],
        [-3.0,  0.0,  5.0, -6.0]
    ])
}

fn get_other_tensor_rank_2_mul_tensor_rank_2_dim_4() -> TensorRank2<4, 1, 1>
{
    TensorRank2::new([
        [75.0, 76.0, 17.0, 81.0],
        [37.0, 32.0, 13.0, 29.0],
        [41.0, 42.0,  9.0, 37.0],
        [51.0, 52.0, 11.0, 45.0]
    ])
}

#[test]
fn need_to_test_configuration_switches_and_compatibility_enforcement()
{
    todo!("how to test the enforcement of not contracting different coordinate systems?")
}

#[test]
fn add_tensor_rank_2_to_self()
{
    (get_tensor_rank_2_dim_4() + get_other_tensor_rank_2_dim_4()).iter()
    .zip(get_other_tensor_rank_2_add_tensor_rank_2_dim_4().iter())
    .for_each(|(tensor_rank_2_i, res_tensor_rank_2_i)|
        tensor_rank_2_i.iter()
        .zip(res_tensor_rank_2_i.iter())
        .for_each(|(tensor_rank_2_ij, res_tensor_rank_2_ij)|
            assert_eq!(tensor_rank_2_ij, res_tensor_rank_2_ij)
        )
    );
}

#[test]
fn add_tensor_rank_2_ref_to_self()
{
    (get_tensor_rank_2_dim_4() + &get_other_tensor_rank_2_dim_4()).iter()
    .zip(get_other_tensor_rank_2_add_tensor_rank_2_dim_4().iter())
    .for_each(|(tensor_rank_2_i, res_tensor_rank_2_i)|
        tensor_rank_2_i.iter()
        .zip(res_tensor_rank_2_i.iter())
        .for_each(|(tensor_rank_2_ij, res_tensor_rank_2_ij)|
            assert_eq!(tensor_rank_2_ij, res_tensor_rank_2_ij)
        )
    );
}

#[test]
fn add_tensor_rank_2_to_self_ref()
{
    (&get_tensor_rank_2_dim_4() + get_other_tensor_rank_2_dim_4()).iter()
    .zip(get_other_tensor_rank_2_add_tensor_rank_2_dim_4().iter())
    .for_each(|(tensor_rank_2_i, res_tensor_rank_2_i)|
        tensor_rank_2_i.iter()
        .zip(res_tensor_rank_2_i.iter())
        .for_each(|(tensor_rank_2_ij, res_tensor_rank_2_ij)|
            assert_eq!(tensor_rank_2_ij, res_tensor_rank_2_ij)
        )
    );
}

#[test]
fn add_assign_tensor_rank_2()
{
    let mut tensor_rank_2 = get_tensor_rank_2_dim_4();
    tensor_rank_2 += get_other_tensor_rank_2_dim_4();
    tensor_rank_2.iter()
    .zip(get_other_tensor_rank_2_add_tensor_rank_2_dim_4().iter())
    .for_each(|(tensor_rank_2_i, res_tensor_rank_2_i)|
        tensor_rank_2_i.iter()
        .zip(res_tensor_rank_2_i.iter())
        .for_each(|(tensor_rank_2_ij, res_tensor_rank_2_ij)|
            assert_eq!(tensor_rank_2_ij, res_tensor_rank_2_ij)
        )
    );
}

#[test]
fn add_assign_tensor_rank_2_ref()
{
    let mut tensor_rank_2 = get_tensor_rank_2_dim_4();
    tensor_rank_2 += &get_other_tensor_rank_2_dim_4();
    tensor_rank_2.iter()
    .zip(get_other_tensor_rank_2_add_tensor_rank_2_dim_4().iter())
    .for_each(|(tensor_rank_2_i, res_tensor_rank_2_i)|
        tensor_rank_2_i.iter()
        .zip(res_tensor_rank_2_i.iter())
        .for_each(|(tensor_rank_2_ij, res_tensor_rank_2_ij)|
            assert_eq!(tensor_rank_2_ij, res_tensor_rank_2_ij)
        )
    );
}

#[test]
fn div_tensor_rank_0_to_self()
{
    (get_tensor_rank_2_dim_4() / 3.3).iter()
    .zip(get_array_dim_4().iter())
    .for_each(|(tensor_rank_2_i, array_i)|
        tensor_rank_2_i.iter()
        .zip(array_i.iter())
        .for_each(|(tensor_rank_2_ij, array_ij)|
            assert_eq!(tensor_rank_2_ij, &(array_ij / 3.3))
        )
    );
}

#[test]
fn div_tensor_rank_0_to_self_ref()
{
    (&get_tensor_rank_2_dim_4() / 3.3).iter()
    .zip(get_array_dim_4().iter())
    .for_each(|(tensor_rank_2_i, array_i)|
        tensor_rank_2_i.iter()
        .zip(array_i.iter())
        .for_each(|(tensor_rank_2_ij, array_ij)|
            assert_eq!(tensor_rank_2_ij, &(array_ij / 3.3))
        )
    );
}

#[test]
fn div_tensor_rank_0_ref_to_self()
{
    (get_tensor_rank_2_dim_4() / &3.3).iter()
    .zip(get_array_dim_4().iter())
    .for_each(|(tensor_rank_2_i, array_i)|
        tensor_rank_2_i.iter()
        .zip(array_i.iter())
        .for_each(|(tensor_rank_2_ij, array_ij)|
            assert_eq!(tensor_rank_2_ij, &(array_ij / 3.3))
        )
    );
}

#[test]
fn div_tensor_rank_0_ref_to_self_ref()
{
    (&get_tensor_rank_2_dim_4() / &3.3).iter()
    .zip(get_array_dim_4().iter())
    .for_each(|(tensor_rank_2_i, array_i)|
        tensor_rank_2_i.iter()
        .zip(array_i.iter())
        .for_each(|(tensor_rank_2_ij, array_ij)|
            assert_eq!(tensor_rank_2_ij, &(array_ij / 3.3))
        )
    );
}

#[test]
fn div_assign_tensor_rank_0()
{
    let mut tensor_rank_2 = get_tensor_rank_2_dim_4();
    tensor_rank_2 /= 3.3;
    tensor_rank_2.iter()
    .zip(get_array_dim_4().iter())
    .for_each(|(tensor_rank_2_i, array_i)|
        tensor_rank_2_i.iter()
        .zip(array_i.iter())
        .for_each(|(tensor_rank_2_ij, array_ij)|
            assert_eq!(tensor_rank_2_ij, &(array_ij / 3.3))
        )
    );
}

#[test]
fn div_assign_tensor_rank_0_ref()
{
    let mut tensor_rank_2 = get_tensor_rank_2_dim_4();
    tensor_rank_2 /= &3.3;
    tensor_rank_2.iter()
    .zip(get_array_dim_4().iter())
    .for_each(|(tensor_rank_2_i, array_i)|
        tensor_rank_2_i.iter()
        .zip(array_i.iter())
        .for_each(|(tensor_rank_2_ij, array_ij)|
            assert_eq!(tensor_rank_2_ij, &(array_ij / 3.3))
        )
    );
}

#[test]
fn determinant_dim_2()
{
    assert_eq_within_tols(&get_tensor_rank_2_dim_2().determinant(), &(-2.0));
}

#[test]
fn determinant_dim_3()
{
    assert_eq_within_tols(&get_tensor_rank_2_dim_3().determinant(), &290.0);
}

#[test]
fn determinant_dim_4()
{
    assert_eq_within_tols(&get_tensor_rank_2_dim_4().determinant(), &36.0);
} 

#[test]
fn determinant_dim_9()
{
    assert_eq_within_tols(&get_tensor_rank_2_dim_9().determinant(), &2398.0);
}

#[test]
fn deviatoric_dim_2()
{
    let tensor_rank_2 = get_tensor_rank_2_dim_2();
    let trace = tensor_rank_2.trace();
    let deviatoric_tensor_rank_2 = tensor_rank_2.deviatoric();
    assert_eq!(deviatoric_tensor_rank_2.trace(), 0.0);
    deviatoric_tensor_rank_2.iter()
    .enumerate()
    .zip(tensor_rank_2.iter())
    .for_each(|((i, deviatoric_tensor_rank_2_i), tensor_rank_2_i)|
        deviatoric_tensor_rank_2_i.iter()
        .enumerate()
        .zip(tensor_rank_2_i.iter())
        .for_each(|((j, deviatoric_tensor_rank_2_ij), tensor_rank_2_ij)|
            assert_eq!(deviatoric_tensor_rank_2_ij, &(tensor_rank_2_ij - (((i == j) as u8) as TensorRank0) * trace / 2.0))
        )
    );
}

#[test]
fn deviatoric_dim_3()
{
    let tensor_rank_2 = get_tensor_rank_2_dim_3();
    let trace = tensor_rank_2.trace();
    let deviatoric_tensor_rank_2 = tensor_rank_2.deviatoric();
    assert_eq!(deviatoric_tensor_rank_2.trace(), 0.0);
    deviatoric_tensor_rank_2.iter()
    .enumerate()
    .zip(tensor_rank_2.iter())
    .for_each(|((i, deviatoric_tensor_rank_2_i), tensor_rank_2_i)|
        deviatoric_tensor_rank_2_i.iter()
        .enumerate()
        .zip(tensor_rank_2_i.iter())
        .for_each(|((j, deviatoric_tensor_rank_2_ij), tensor_rank_2_ij)|
            assert_eq!(deviatoric_tensor_rank_2_ij, &(tensor_rank_2_ij - (((i == j) as u8) as TensorRank0) * trace / 3.0))
        )
    );
}

#[test]
fn deviatoric_dim_4()
{
    let tensor_rank_2 = get_tensor_rank_2_dim_4();
    let trace = tensor_rank_2.trace();
    let deviatoric_tensor_rank_2 = tensor_rank_2.deviatoric();
    assert_eq!(deviatoric_tensor_rank_2.trace(), 0.0);
    deviatoric_tensor_rank_2.iter()
    .enumerate()
    .zip(tensor_rank_2.iter())
    .for_each(|((i, deviatoric_tensor_rank_2_i), tensor_rank_2_i)|
        deviatoric_tensor_rank_2_i.iter()
        .enumerate()
        .zip(tensor_rank_2_i.iter())
        .for_each(|((j, deviatoric_tensor_rank_2_ij), tensor_rank_2_ij)|
            assert_eq!(deviatoric_tensor_rank_2_ij, &(tensor_rank_2_ij - (((i == j) as u8) as TensorRank0) * trace / 4.0))
        )
    );
}

#[test]
fn deviatoric_dim_9()
{
    let tensor_rank_2 = get_tensor_rank_2_dim_9();
    let trace = tensor_rank_2.trace();
    let deviatoric_tensor_rank_2 = tensor_rank_2.deviatoric();
    assert_eq!(deviatoric_tensor_rank_2.trace(), 0.0);
    deviatoric_tensor_rank_2.iter()
    .enumerate()
    .zip(tensor_rank_2.iter())
    .for_each(|((i, deviatoric_tensor_rank_2_i), tensor_rank_2_i)|
        deviatoric_tensor_rank_2_i.iter()
        .enumerate()
        .zip(tensor_rank_2_i.iter())
        .for_each(|((j, deviatoric_tensor_rank_2_ij), tensor_rank_2_ij)|
            assert_eq!(deviatoric_tensor_rank_2_ij, &(tensor_rank_2_ij - (((i == j) as u8) as TensorRank0) * trace / 9.0))
        )
    );
}

#[test]
fn deviatoric_and_trace_dim_2()
{
    let tensor_rank_2 = get_tensor_rank_2_dim_2();
    let (deviatoric, trace) = tensor_rank_2.deviatoric_and_trace();
    assert_eq!(trace, tensor_rank_2.trace());
    deviatoric.iter()
    .zip(tensor_rank_2.deviatoric().iter())
    .for_each(|(deviatoric_i, tensor_rank_2_deviatoric_i)|
        deviatoric_i.iter()
        .zip(tensor_rank_2_deviatoric_i.iter())
        .for_each(|(deviatoric_ij, tensor_rank_2_deviatoric_ij)|
        assert_eq!(deviatoric_ij, tensor_rank_2_deviatoric_ij)
        )
    );
}

#[test]
fn deviatoric_and_trace_dim_3()
{
    let tensor_rank_2 = get_tensor_rank_2_dim_3();
    let (deviatoric, trace) = tensor_rank_2.deviatoric_and_trace();
    assert_eq!(trace, tensor_rank_2.trace());
    deviatoric.iter()
    .zip(tensor_rank_2.deviatoric().iter())
    .for_each(|(deviatoric_i, tensor_rank_2_deviatoric_i)|
        deviatoric_i.iter()
        .zip(tensor_rank_2_deviatoric_i.iter())
        .for_each(|(deviatoric_ij, tensor_rank_2_deviatoric_ij)|
        assert_eq!(deviatoric_ij, tensor_rank_2_deviatoric_ij)
        )
    );
}

#[test]
fn deviatoric_and_trace_dim_4()
{
    let tensor_rank_2 = get_tensor_rank_2_dim_4();
    let (deviatoric, trace) = tensor_rank_2.deviatoric_and_trace();
    assert_eq!(trace, tensor_rank_2.trace());
    deviatoric.iter()
    .zip(tensor_rank_2.deviatoric().iter())
    .for_each(|(deviatoric_i, tensor_rank_2_deviatoric_i)|
        deviatoric_i.iter()
        .zip(tensor_rank_2_deviatoric_i.iter())
        .for_each(|(deviatoric_ij, tensor_rank_2_deviatoric_ij)|
        assert_eq!(deviatoric_ij, tensor_rank_2_deviatoric_ij)
        )
    );
}

#[test]
fn deviatoric_and_trace_dim_9()
{
    let tensor_rank_2 = get_tensor_rank_2_dim_9();
    let (deviatoric, trace) = tensor_rank_2.deviatoric_and_trace();
    assert_eq!(trace, tensor_rank_2.trace());
    deviatoric.iter()
    .zip(tensor_rank_2.deviatoric().iter())
    .for_each(|(deviatoric_i, tensor_rank_2_deviatoric_i)|
        deviatoric_i.iter()
        .zip(tensor_rank_2_deviatoric_i.iter())
        .for_each(|(deviatoric_ij, tensor_rank_2_deviatoric_ij)|
        assert_eq!(deviatoric_ij, tensor_rank_2_deviatoric_ij)
        )
    );
}

#[test]
fn dyad()
{
    let tensor_rank_1_a = get_tensor_rank_1_a();
    let tensor_rank_1_b = get_tensor_rank_1_b();
    let tensor_rank_2 = TensorRank2::dyad(&tensor_rank_1_a, &tensor_rank_1_b);
    tensor_rank_2.iter()
    .zip(tensor_rank_1_a.iter())
    .for_each(|(tensor_rank_2_i, tensor_rank_1_a_i)|
        tensor_rank_2_i.iter()
        .zip(tensor_rank_1_b.iter())
        .for_each(|(tensor_rank_2_ij, tensor_rank_1_b_j)|
            assert_eq!(tensor_rank_2_ij, &(tensor_rank_1_a_i * tensor_rank_1_b_j))
        )
    );
}

#[test]
fn from_iter()
{
    let into_iterator = get_tensor_rank_2_dim_4().0.into_iter();
    let tensor_rank_2 = TensorRank2::<4, 1, 1>::from_iter(get_tensor_rank_2_dim_4().0.into_iter());
    tensor_rank_2.iter()
    .zip(into_iterator)
    .for_each(|(tensor_rank_2_i, value_i)|
        tensor_rank_2_i.iter()
        .zip(value_i.iter())
        .for_each(|(tensor_rank_2_ij, value_ij)|
            assert_eq!(tensor_rank_2_ij, value_ij)
        )
    );
}

#[test]
fn full_contraction_dim_2()
{
    assert_eq_within_tols(&get_tensor_rank_2_dim_2().full_contraction(&get_other_tensor_rank_2_dim_2()), &70.0);
}

#[test]
fn full_contraction_dim_3()
{
    assert_eq_within_tols(&get_tensor_rank_2_dim_3().full_contraction(&get_other_tensor_rank_2_dim_3()), &167.0);
}

#[test]
fn full_contraction_dim_4()
{
    assert_eq_within_tols(&get_tensor_rank_2_dim_4().full_contraction(&get_other_tensor_rank_2_dim_4()), &137.0);
}

#[test]
fn full_contraction_dim_9()
{
    assert_eq_within_tols(&get_tensor_rank_2_dim_9().full_contraction(&get_other_tensor_rank_2_dim_9()), &269.0);
}

#[test]
fn identity()
{
    TensorRank2::<9, 1, 1>::identity().iter()
    .enumerate()
    .for_each(|(i, tensor_rank_2_i)|
        tensor_rank_2_i.iter()
        .enumerate()
        .for_each(|(j, tensor_rank_2_ij)|
            if i == j
            {
                assert_eq!(tensor_rank_2_ij, &1.0)
            }
            else
            {
                assert_eq!(tensor_rank_2_ij, &0.0)
            }
        )
    );
}

#[test]
fn inverse_dim_2()
{
    (get_tensor_rank_2_dim_2() * get_tensor_rank_2_dim_2().inverse()).iter()
    .enumerate()
    .for_each(|(i, tensor_rank_2_i)|
        tensor_rank_2_i.iter()
        .enumerate()
        .for_each(|(j, tensor_rank_2_ij)|
            if i == j
            {
                assert_eq_within_tols(tensor_rank_2_ij, &1.0)
            }
            else
            {
                assert_eq_within_tols(tensor_rank_2_ij, &0.0)
            }
        )
    );
}

#[test]
fn inverse_dim_3()
{
    (get_tensor_rank_2_dim_3() * get_tensor_rank_2_dim_3().inverse()).iter()
    .enumerate()
    .for_each(|(i, tensor_rank_2_i)|
        tensor_rank_2_i.iter()
        .enumerate()
        .for_each(|(j, tensor_rank_2_ij)|
            if i == j
            {
                assert_eq_within_tols(tensor_rank_2_ij, &1.0)
            }
            else
            {
                assert_eq_within_tols(tensor_rank_2_ij, &0.0)
            }
        )
    );
}

#[test]
fn inverse_dim_4()
{
    (get_tensor_rank_2_dim_4() * get_tensor_rank_2_dim_4().inverse()).iter()
    .enumerate()
    .for_each(|(i, tensor_rank_2_i)|
        tensor_rank_2_i.iter()
        .enumerate()
        .for_each(|(j, tensor_rank_2_ij)|
            if i == j
            {
                assert_eq_within_tols(tensor_rank_2_ij, &1.0)
            }
            else
            {
                assert_eq_within_tols(tensor_rank_2_ij, &0.0)
            }
        )
    );
}

#[test]
fn inverse_dim_9()
{
    (get_tensor_rank_2_dim_9() * get_tensor_rank_2_dim_9().inverse()).iter()
    .enumerate()
    .for_each(|(i, tensor_rank_2_i)|
        tensor_rank_2_i.iter()
        .enumerate()
        .for_each(|(j, tensor_rank_2_ij)|
            if i == j
            {
                assert_eq_within_tols(tensor_rank_2_ij, &1.0)
            }
            else
            {
                assert_eq_within_tols(tensor_rank_2_ij, &0.0)
            }
        )
    );
}

#[test]
fn inverse_and_determinant_dim_2()
{
    let tensor_rank_2 = get_tensor_rank_2_dim_2();
    let (inverse, determinant) = tensor_rank_2.inverse_and_determinant();
    assert_eq!(determinant, tensor_rank_2.determinant());
    inverse.iter()
    .zip(tensor_rank_2.inverse().iter())
    .for_each(|(inverse_i, tensor_rank_2_inverse_i)|
        inverse_i.iter()
        .zip(tensor_rank_2_inverse_i.iter())
        .for_each(|(inverse_ij, tensor_rank_2_inverse_ij)|
        assert_eq!(inverse_ij, tensor_rank_2_inverse_ij)
        )
    );
}

#[test]
fn inverse_and_determinant_dim_3()
{
    let tensor_rank_2 = get_tensor_rank_2_dim_3();
    let (inverse, determinant) = tensor_rank_2.inverse_and_determinant();
    assert_eq!(determinant, tensor_rank_2.determinant());
    inverse.iter()
    .zip(tensor_rank_2.inverse().iter())
    .for_each(|(inverse_i, tensor_rank_2_inverse_i)|
        inverse_i.iter()
        .zip(tensor_rank_2_inverse_i.iter())
        .for_each(|(inverse_ij, tensor_rank_2_inverse_ij)|
        assert_eq!(inverse_ij, tensor_rank_2_inverse_ij)
        )
    );
}

#[test]
fn inverse_and_determinant_dim_4()
{
    let tensor_rank_2 = get_tensor_rank_2_dim_4();
    let (inverse, determinant) = tensor_rank_2.inverse_and_determinant();
    assert_eq!(determinant, tensor_rank_2.determinant());
    inverse.iter()
    .zip(tensor_rank_2.inverse().iter())
    .for_each(|(inverse_i, tensor_rank_2_inverse_i)|
        inverse_i.iter()
        .zip(tensor_rank_2_inverse_i.iter())
        .for_each(|(inverse_ij, tensor_rank_2_inverse_ij)|
        assert_eq!(inverse_ij, tensor_rank_2_inverse_ij)
        )
    );
}

#[test]
fn inverse_transpose_dim_2()
{
    (get_tensor_rank_2_dim_2().transpose() * get_tensor_rank_2_dim_2().inverse_transpose()).iter()
    .enumerate()
    .for_each(|(i, tensor_rank_2_i)|
        tensor_rank_2_i.iter()
        .enumerate()
        .for_each(|(j, tensor_rank_2_ij)|
            if i == j
            {
                assert_eq!(tensor_rank_2_ij, &1.0)
            }
            else
            {
                assert_eq!(tensor_rank_2_ij, &0.0)
            }
        )
    );
}

#[test]
fn inverse_transpose_dim_3()
{
    (get_tensor_rank_2_dim_3().transpose() * get_tensor_rank_2_dim_3().inverse_transpose()).iter()
    .enumerate()
    .for_each(|(i, tensor_rank_2_i)|
        tensor_rank_2_i.iter()
        .enumerate()
        .for_each(|(j, tensor_rank_2_ij)|
            if i == j
            {
                assert_eq_within_tols(tensor_rank_2_ij, &1.0)
            }
            else
            {
                assert_eq_within_tols(tensor_rank_2_ij, &0.0)
            }
        )
    );
}

#[test]
fn inverse_transpose_dim_4()
{
    (get_tensor_rank_2_dim_4().transpose() * get_tensor_rank_2_dim_4().inverse_transpose()).iter()
    .enumerate()
    .for_each(|(i, tensor_rank_2_i)|
        tensor_rank_2_i.iter()
        .enumerate()
        .for_each(|(j, tensor_rank_2_ij)|
            if i == j
            {
                assert_eq_within_tols(tensor_rank_2_ij, &1.0)
            }
            else
            {
                assert_eq_within_tols(tensor_rank_2_ij, &0.0)
            }
        )
    );
}

#[test]
fn inverse_transpose_9()
{
    (get_tensor_rank_2_dim_9().transpose() * get_tensor_rank_2_dim_9().inverse_transpose()).iter()
    .enumerate()
    .for_each(|(i, tensor_rank_2_i)|
        tensor_rank_2_i.iter()
        .enumerate()
        .for_each(|(j, tensor_rank_2_ij)|
            if i == j
            {
                assert_eq_within_tols(tensor_rank_2_ij, &1.0)
            }
            else
            {
                assert_eq_within_tols(tensor_rank_2_ij, &0.0)
            }
        )
    );
}

#[test]
fn inverse_transpose_and_determinant_dim_2()
{
    let tensor_rank_2 = get_tensor_rank_2_dim_2();
    let (inverse_transpose, determinant) = tensor_rank_2.inverse_transpose_and_determinant();
    assert_eq!(determinant, tensor_rank_2.determinant());
    inverse_transpose.iter()
    .zip(tensor_rank_2.inverse_transpose().iter())
    .for_each(|(inverse_transpose_i, tensor_rank_2_inverse_transpose_i)|
    inverse_transpose_i.iter()
        .zip(tensor_rank_2_inverse_transpose_i.iter())
        .for_each(|(inverse_transpose_ij, tensor_rank_2_inverse_transpose_ij)|
            assert_eq!(inverse_transpose_ij, tensor_rank_2_inverse_transpose_ij)
        )
    );
}

#[test]
fn inverse_transpose_and_determinant_dim_3()
{
    let tensor_rank_2 = get_tensor_rank_2_dim_3();
    let (inverse_transpose, determinant) = tensor_rank_2.inverse_transpose_and_determinant();
    assert_eq!(determinant, tensor_rank_2.determinant());
    inverse_transpose.iter()
    .zip(tensor_rank_2.inverse_transpose().iter())
    .for_each(|(inverse_transpose_i, tensor_rank_2_inverse_transpose_i)|
    inverse_transpose_i.iter()
        .zip(tensor_rank_2_inverse_transpose_i.iter())
        .for_each(|(inverse_transpose_ij, tensor_rank_2_inverse_transpose_ij)|
            assert_eq!(inverse_transpose_ij, tensor_rank_2_inverse_transpose_ij)
        )
    );
}

#[test]
fn inverse_transpose_and_determinant_dim_4()
{
    let tensor_rank_2 = get_tensor_rank_2_dim_4();
    let (inverse_transpose, determinant) = tensor_rank_2.inverse_transpose_and_determinant();
    assert_eq!(determinant, tensor_rank_2.determinant());
    inverse_transpose.iter()
    .zip(tensor_rank_2.inverse_transpose().iter())
    .for_each(|(inverse_transpose_i, tensor_rank_2_inverse_transpose_i)|
    inverse_transpose_i.iter()
        .zip(tensor_rank_2_inverse_transpose_i.iter())
        .for_each(|(inverse_transpose_ij, tensor_rank_2_inverse_transpose_ij)|
            assert_eq!(inverse_transpose_ij, tensor_rank_2_inverse_transpose_ij)
        )
    );
}

#[test]
fn iter()
{
    get_tensor_rank_2_dim_4().iter()
    .zip(get_array_dim_4().iter())
    .for_each(|(tensor_rank_2_i, array_i)|
        tensor_rank_2_i.iter()
        .zip(array_i.iter())
        .for_each(|(tensor_rank_2_ij, array_ij)|
            assert_eq!(tensor_rank_2_ij, array_ij)
        )
    );
}

#[test]
fn iter_mut()
{
    get_tensor_rank_2_dim_4().iter_mut()
    .zip(get_array_dim_4().iter_mut())
    .for_each(|(tensor_rank_2_i, array_i)|
        tensor_rank_2_i.iter_mut()
        .zip(array_i.iter_mut())
        .for_each(|(tensor_rank_2_ij, array_ij)|
            assert_eq!(tensor_rank_2_ij, array_ij)
        )
    );
}

#[test]
fn lu_decomposition()
{
    let (tensor_l, tensor_u) = get_tensor_rank_2_dim_9().lu_decomposition();
    tensor_l.iter()
    .enumerate()
    .zip(tensor_u.iter())
    .for_each(|((i, tensor_l_i), tensor_u_i)|
        tensor_l_i.iter()
        .enumerate()
        .zip(tensor_u_i.iter())
        .for_each(|((j, tensor_l_ij), tensor_u_ij)|
            if i > j
            {
                assert_eq!(tensor_u_ij, &0.0);
            }
            else if i < j
            {
                assert_eq!(tensor_l_ij, &0.0);
            }
        )
    );
}

#[test]
fn mul_tensor_rank_0_to_self()
{
    (get_tensor_rank_2_dim_4() * 3.3).iter()
    .zip(get_array_dim_4().iter())
    .for_each(|(tensor_rank_2_i, array_i)|
        tensor_rank_2_i.iter()
        .zip(array_i.iter())
        .for_each(|(tensor_rank_2_ij, array_ij)|
            assert_eq!(tensor_rank_2_ij, &(array_ij * 3.3))
        )
    );
}

#[test]
fn mul_tensor_rank_0_to_self_ref()
{
    (&get_tensor_rank_2_dim_4() * 3.3).iter()
    .zip(get_array_dim_4().iter())
    .for_each(|(tensor_rank_2_i, array_i)|
        tensor_rank_2_i.iter()
        .zip(array_i.iter())
        .for_each(|(tensor_rank_2_ij, array_ij)|
            assert_eq!(tensor_rank_2_ij, &(array_ij * 3.3))
        )
    );
}

#[test]
fn mul_tensor_rank_0_ref_to_self()
{
    (get_tensor_rank_2_dim_4() * &3.3).iter()
    .zip(get_array_dim_4().iter())
    .for_each(|(tensor_rank_2_i, array_i)|
        tensor_rank_2_i.iter()
        .zip(array_i.iter())
        .for_each(|(tensor_rank_2_ij, array_ij)|
            assert_eq!(tensor_rank_2_ij, &(array_ij * 3.3))
        )
    );
}

#[test]
fn mul_tensor_rank_0_ref_to_self_ref()
{
    (&get_tensor_rank_2_dim_4() * &3.3).iter()
    .zip(get_array_dim_4().iter())
    .for_each(|(tensor_rank_2_i, array_i)|
        tensor_rank_2_i.iter()
        .zip(array_i.iter())
        .for_each(|(tensor_rank_2_ij, array_ij)|
            assert_eq!(tensor_rank_2_ij, &(array_ij * 3.3))
        )
    );
}

#[test]
fn mul_assign_tensor_rank_0()
{
    let mut tensor_rank_2 = get_tensor_rank_2_dim_4();
    tensor_rank_2 *= 3.3;
    tensor_rank_2.iter()
    .zip(get_array_dim_4().iter())
    .for_each(|(tensor_rank_2_i, array_i)|
        tensor_rank_2_i.iter()
        .zip(array_i.iter())
        .for_each(|(tensor_rank_2_ij, array_ij)|
            assert_eq!(tensor_rank_2_ij, &(array_ij * 3.3))
        )
    );
}

#[test]
fn mul_assign_tensor_rank_0_ref()
{
    let mut tensor_rank_2 = get_tensor_rank_2_dim_4();
    tensor_rank_2 *= &3.3;
    tensor_rank_2.iter()
    .zip(get_array_dim_4().iter())
    .for_each(|(tensor_rank_2_i, array_i)|
        tensor_rank_2_i.iter()
        .zip(array_i.iter())
        .for_each(|(tensor_rank_2_ij, array_ij)|
            assert_eq!(tensor_rank_2_ij, &(array_ij * 3.3))
        )
    );
}

#[test]
fn mul_tensor_rank_1_to_self()
{
    (get_tensor_rank_2_dim_4() * get_tensor_rank_1_a()).iter()
    .zip(get_other_tensor_rank_2_mul_tensor_rank_1_dim_4().iter())
    .for_each(|(tensor_rank_1_i, res_tensor_rank_1_i)|
        assert_eq!(tensor_rank_1_i, res_tensor_rank_1_i)
    );
}

#[test]
fn mul_tensor_rank_1_ref_to_self()
{
    (get_tensor_rank_2_dim_4() * &get_tensor_rank_1_a()).iter()
    .zip(get_other_tensor_rank_2_mul_tensor_rank_1_dim_4().iter())
    .for_each(|(tensor_rank_1_i, res_tensor_rank_1_i)|
        assert_eq!(tensor_rank_1_i, res_tensor_rank_1_i)
    );
}

#[test]
fn mul_tensor_rank_2_to_self()
{
    (get_tensor_rank_2_dim_4() * get_other_tensor_rank_2_dim_4()).iter()
    .zip(get_other_tensor_rank_2_mul_tensor_rank_2_dim_4().iter())
    .for_each(|(tensor_rank_2_i, res_tensor_rank_2_i)|
        tensor_rank_2_i.iter()
        .zip(res_tensor_rank_2_i.iter())
        .for_each(|(tensor_rank_2_ij, res_tensor_rank_2_ij)|
            assert_eq!(tensor_rank_2_ij, res_tensor_rank_2_ij)
        )
    );
}

#[test]
fn mul_tensor_rank_2_ref_to_self()
{
    (get_tensor_rank_2_dim_4() * &get_other_tensor_rank_2_dim_4()).iter()
    .zip(get_other_tensor_rank_2_mul_tensor_rank_2_dim_4().iter())
    .for_each(|(tensor_rank_2_i, res_tensor_rank_2_i)|
        tensor_rank_2_i.iter()
        .zip(res_tensor_rank_2_i.iter())
        .for_each(|(tensor_rank_2_ij, res_tensor_rank_2_ij)|
            assert_eq!(tensor_rank_2_ij, res_tensor_rank_2_ij)
        )
    );
}

#[test]
fn mul_tensor_rank_2_to_self_ref()
{
    (&get_tensor_rank_2_dim_4() * get_other_tensor_rank_2_dim_4()).iter()
    .zip(get_other_tensor_rank_2_mul_tensor_rank_2_dim_4().iter())
    .for_each(|(tensor_rank_2_i, res_tensor_rank_2_i)|
        tensor_rank_2_i.iter()
        .zip(res_tensor_rank_2_i.iter())
        .for_each(|(tensor_rank_2_ij, res_tensor_rank_2_ij)|
            assert_eq!(tensor_rank_2_ij, res_tensor_rank_2_ij)
        )
    );
}

#[test]
fn mul_tensor_rank_2_ref_to_self_ref()
{
    (&get_tensor_rank_2_dim_4() * &get_other_tensor_rank_2_dim_4()).iter()
    .zip(get_other_tensor_rank_2_mul_tensor_rank_2_dim_4().iter())
    .for_each(|(tensor_rank_2_i, res_tensor_rank_2_i)|
        tensor_rank_2_i.iter()
        .zip(res_tensor_rank_2_i.iter())
        .for_each(|(tensor_rank_2_ij, res_tensor_rank_2_ij)|
            assert_eq!(tensor_rank_2_ij, res_tensor_rank_2_ij)
        )
    );
}

#[test]
fn new()
{
    get_tensor_rank_2_dim_4().iter()
    .zip(get_array_dim_4().iter())
    .for_each(|(tensor_rank_2_i, array_i)|
        tensor_rank_2_i.iter()
        .zip(array_i.iter())
        .for_each(|(tensor_rank_2_ij, array_ij)|
            assert_eq!(tensor_rank_2_ij, array_ij)
        )
    );
}

#[test]
fn norm_dim_2()
{
    assert_eq!(get_tensor_rank_2_dim_2().norm(), 3.872_983_346_207_417);
}

#[test]
fn norm_dim_3()
{
    assert_eq!(get_tensor_rank_2_dim_3().norm(), 11.937_336_386_313_323);
}

#[test]
fn norm_dim_4()
{
    assert_eq!(get_tensor_rank_2_dim_4().norm(), 10.099_504_938_362_077);
}

#[test]
fn norm_dim_9()
{
    assert_eq!(get_tensor_rank_2_dim_9().norm(), 14.832_396_974_191_326);
}

#[test]
fn second_invariant()
{
    assert_eq!(get_tensor_rank_2_dim_4().second_invariant(), 16.0);
}

#[test]
fn squared_trace_dim_2()
{
    assert_eq_within_tols(&get_tensor_rank_2_dim_2().squared_trace(), &29.0);
}

#[test]
fn squared_trace_dim_3()
{
    assert_eq_within_tols(&get_tensor_rank_2_dim_3().squared_trace(), &258.0);
}

#[test]
fn squared_trace_dim_4()
{
    assert_eq_within_tols(&get_tensor_rank_2_dim_4().squared_trace(), &89.0);
}

#[test]
fn squared_trace_dim_9()
{
    assert_eq_within_tols(&get_tensor_rank_2_dim_9().squared_trace(), &318.0);
}

#[test]
fn sub_tensor_rank_2_to_self()
{
    (get_tensor_rank_2_dim_4() - get_other_tensor_rank_2_dim_4()).iter()
    .zip(get_other_tensor_rank_2_sub_tensor_rank_2_dim_4().iter())
    .for_each(|(tensor_rank_2_i, res_tensor_rank_2_i)|
        tensor_rank_2_i.iter()
        .zip(res_tensor_rank_2_i.iter())
        .for_each(|(tensor_rank_2_ij, res_tensor_rank_2_ij)|
            assert_eq!(tensor_rank_2_ij, res_tensor_rank_2_ij)
        )
    );
}

#[test]
fn sub_tensor_rank_2_ref_to_self()
{
    (get_tensor_rank_2_dim_4() - &get_other_tensor_rank_2_dim_4()).iter()
    .zip(get_other_tensor_rank_2_sub_tensor_rank_2_dim_4().iter())
    .for_each(|(tensor_rank_2_i, res_tensor_rank_2_i)|
        tensor_rank_2_i.iter()
        .zip(res_tensor_rank_2_i.iter())
        .for_each(|(tensor_rank_2_ij, res_tensor_rank_2_ij)|
            assert_eq!(tensor_rank_2_ij, res_tensor_rank_2_ij)
        )
    );
}

#[test]
fn sub_assign_tensor_rank_2()
{
    let mut tensor_rank_2 = get_tensor_rank_2_dim_4();
    tensor_rank_2 -= get_other_tensor_rank_2_dim_4();
    tensor_rank_2.iter()
    .zip(get_other_tensor_rank_2_sub_tensor_rank_2_dim_4().iter())
    .for_each(|(tensor_rank_2_i, res_tensor_rank_2_i)|
        tensor_rank_2_i.iter()
        .zip(res_tensor_rank_2_i.iter())
        .for_each(|(tensor_rank_2_ij, res_tensor_rank_2_ij)|
            assert_eq!(tensor_rank_2_ij, res_tensor_rank_2_ij)
        )
    );
}

#[test]
fn sub_assign_tensor_rank_2_ref()
{
    let mut tensor_rank_2 = get_tensor_rank_2_dim_4();
    tensor_rank_2 -= &get_other_tensor_rank_2_dim_4();
    tensor_rank_2.iter()
    .zip(get_other_tensor_rank_2_sub_tensor_rank_2_dim_4().iter())
    .for_each(|(tensor_rank_2_i, res_tensor_rank_2_i)|
        tensor_rank_2_i.iter()
        .zip(res_tensor_rank_2_i.iter())
        .for_each(|(tensor_rank_2_ij, res_tensor_rank_2_ij)|
            assert_eq!(tensor_rank_2_ij, res_tensor_rank_2_ij)
        )
    );
}

#[test]
fn trace_dim_2()
{
    assert_eq!(get_tensor_rank_2_dim_2().trace(), 5.0);
}

#[test]
fn trace_dim_3()
{
    assert_eq!(get_tensor_rank_2_dim_3().trace(), 6.0);
}

#[test]
fn trace_dim_4()
{
    assert_eq!(get_tensor_rank_2_dim_4().trace(), 11.0);
}

#[test]
fn trace_dim_9()
{
    assert_eq!(get_tensor_rank_2_dim_9().trace(), 18.0);
}

#[test]
fn transpose()
{
    let tensor_rank_2 = get_tensor_rank_2_dim_4();
    let tensor_rank_2_transpose = tensor_rank_2.transpose();
    tensor_rank_2.iter()
    .enumerate()
    .for_each(|(i, tensor_rank_2_i)|
        tensor_rank_2_i.iter()
        .enumerate()
        .for_each(|(j, tensor_rank_2_ij)|
            assert_eq!(tensor_rank_2_ij, &tensor_rank_2_transpose[j][i])
        )
    );
    tensor_rank_2_transpose.iter()
    .enumerate()
    .for_each(|(i, tensor_rank_2_transpose_i)|
        tensor_rank_2_transpose_i.iter()
        .enumerate()
        .for_each(|(j, tensor_rank_2_transpose_ij)|
            assert_eq!(tensor_rank_2_transpose_ij, &tensor_rank_2[j][i])
        )
    );
}

#[test]
fn zero_dim_2()
{
    TensorRank2::<2, 1, 1>::zero().iter()
    .for_each(|tensor_rank_2_i|
        tensor_rank_2_i.iter()
        .for_each(|tensor_rank_2_ij|
            assert_eq!(tensor_rank_2_ij, &0.0)
        )
    );
}

#[test]
fn zero_dim_3()
{
    TensorRank2::<3, 1, 1>::zero().iter()
    .for_each(|tensor_rank_2_i|
        tensor_rank_2_i.iter()
        .for_each(|tensor_rank_2_ij|
            assert_eq!(tensor_rank_2_ij, &0.0)
        )
    );
}

#[test]
fn zero_dim_4()
{
    TensorRank2::<4, 1, 1>::zero().iter()
    .for_each(|tensor_rank_2_i|
        tensor_rank_2_i.iter()
        .for_each(|tensor_rank_2_ij|
            assert_eq!(tensor_rank_2_ij, &0.0)
        )
    );
}

#[test]
fn zero_dim_9()
{
    TensorRank2::<9, 1, 1>::zero().iter()
    .for_each(|tensor_rank_2_i|
        tensor_rank_2_i.iter()
        .for_each(|tensor_rank_2_ij|
            assert_eq!(tensor_rank_2_ij, &0.0)
        )
    );
}