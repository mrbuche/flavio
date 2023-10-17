use crate::test::assert_eq_within_tols;
use super::
{
    TensorRank0,
    TensorRank1,
    TensorRank1Traits,
    TensorRank2,
    TensorRank2Traits
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

fn get_tensor_rank_1_a() -> TensorRank1<4>
{
    TensorRank1::new([1.0, 2.0, 3.0, 4.0])
}

fn get_tensor_rank_1_b() -> TensorRank1<4>
{
    TensorRank1::new([5.0, 7.0, 6.0, 8.0])
}

fn get_tensor_rank_2_dim_2() -> TensorRank2<2>
{
    TensorRank2::new(get_array_dim_2())
}

fn get_tensor_rank_2_dim_3() -> TensorRank2<3>
{
    TensorRank2::new(get_array_dim_3())
}

fn get_tensor_rank_2_dim_4() -> TensorRank2<4>
{
    TensorRank2::new(get_array_dim_4())
}

fn get_tensor_rank_2_dim_9() -> TensorRank2<9>
{
    TensorRank2::new(get_array_dim_9())
}

fn get_tensor_rank_2_dim_4_squared() -> TensorRank2<4>
{
    TensorRank2::new([
        [17.0, 66.0, 76.0, 6.0],
        [ 7.0, 32.0, 16.0, 6.0],
        [ 9.0, 34.0, 34.0, 6.0],
        [11.0, 42.0, 40.0, 6.0]
    ])
}

#[test]
fn todo()
{
    todo!("Test ops and refs variants with ranks 0, 1, 2.");
}

// #[test]
// fn copy()
// {
//     let tensor_rank_2_a = get_tensor_rank_2_dim_4();
//     let tensor_rank_2_b = tensor_rank_2_a;
//     tensor_rank_2_a.iter()
//     .zip(tensor_rank_2_b.iter())
//     .for_each(|(tensor_rank_2_a_i, tensor_rank_2_b_i)|
//         tensor_rank_2_a_i.iter()
//         .zip(tensor_rank_2_b_i.iter())
//         .for_each(|(tensor_rank_2_a_ij, tensor_rank_2_b_ij)|
//             assert_eq!(tensor_rank_2_a_ij, tensor_rank_2_b_ij)
//         )
//     );
// }

#[test]
fn determinant_dim_2()
{
    assert_eq!(get_tensor_rank_2_dim_2().determinant(), -2.0);
}

#[test]
fn determinant_dim_3()
{
    assert_eq!(get_tensor_rank_2_dim_3().determinant(), 290.0);
}

#[test]
fn determinant_dim_4()
{
    assert_eq!(get_tensor_rank_2_dim_4().determinant(), 36.0);
} 

#[test]
fn deviatoric()
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
    let tensor_rank_2 = TensorRank2::<4>::from_iter(get_tensor_rank_2_dim_4().0.into_iter());
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
fn full_contraction()
{
    assert_eq_within_tols(
        &(get_tensor_rank_2_dim_4().full_contraction(&(get_tensor_rank_2_dim_4() / 3.3))),
        &(get_tensor_rank_2_dim_4().norm().powi(2) * 2.0 / 3.3)
    )
}

#[test]
fn identity()
{
    TensorRank2::<9>::identity().iter()
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
fn index()
{
    get_tensor_rank_2_dim_4().iter()
    .zip(get_array_dim_4().iter())
    .for_each(|(tensor_rank_2_i, array_i)|
        assert_eq!(tensor_rank_2_i.0, *array_i)
    );
}

#[test]
fn index_mut()
{
    get_tensor_rank_2_dim_4().iter_mut()
    .zip(get_array_dim_4().iter_mut())
    .for_each(|(tensor_rank_2_i, array_i)|
        assert_eq!(tensor_rank_2_i.0, *array_i)
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
fn inverse_lower_triangular()
{
    let tensor_l = get_tensor_rank_2_dim_9().lu_decomposition().0;
    let inverse_tensor_l = get_tensor_rank_2_dim_9().lu_decomposition().0.inverse_lower_triangular();
    (&inverse_tensor_l * &tensor_l).iter()
    .enumerate()
    .zip(inverse_tensor_l.iter())
    .for_each(|((i, tensor_rank_2_i), inverse_tensor_l_i)|
        tensor_rank_2_i.iter()
        .enumerate()
        .zip(inverse_tensor_l_i.iter())
        .for_each(|((j, tensor_rank_2_ij), inverse_tensor_l_ij)|
            if i == j
            {
                assert_eq_within_tols(tensor_rank_2_ij, &1.0)
            }
            else if i < j
            {
                assert_eq!(inverse_tensor_l_ij, &0.0);
                assert_eq_within_tols(tensor_rank_2_ij, &0.0)
            }
            else
            {
                assert_eq_within_tols(tensor_rank_2_ij, &0.0)
            }
        )
    );
}

#[test]
fn inverse_upper_triangular()
{
    let tensor_u = get_tensor_rank_2_dim_9().lu_decomposition().1;
    let inverse_tensor_u = get_tensor_rank_2_dim_9().lu_decomposition().1.inverse_upper_triangular();
    (&inverse_tensor_u * &tensor_u).iter()
    .enumerate()
    .zip(inverse_tensor_u.iter())
    .for_each(|((i, tensor_rank_2_i), inverse_tensor_u_i)|
        tensor_rank_2_i.iter()
        .enumerate()
        .zip(inverse_tensor_u_i.iter())
        .for_each(|((j, tensor_rank_2_ij), inverse_tensor_u_ij)|
            if i == j
            {
                assert_eq_within_tols(tensor_rank_2_ij, &1.0)
            }
            else if i > j
            {
                assert_eq!(inverse_tensor_u_ij, &0.0);
                assert_eq_within_tols(tensor_rank_2_ij, &0.0)
            }
            else
            {
                assert_eq_within_tols(tensor_rank_2_ij, &0.0)
            }
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
fn norm()
{
    assert_eq!(get_tensor_rank_2_dim_4().norm(), 10.099_504_938_362_077);
}

#[test]
fn second_invariant()
{
    assert_eq!(get_tensor_rank_2_dim_4().second_invariant(), 16.0);
}

#[test]
fn squared()
{
    get_tensor_rank_2_dim_4().squared().iter()
    .zip(get_tensor_rank_2_dim_4_squared().iter())
    .for_each(|(tensor_rank_2_i, squared_tensor_rank_2_i)|
        tensor_rank_2_i.iter()
        .zip(squared_tensor_rank_2_i.iter())
        .for_each(|(tensor_rank_2_ij, squared_tensor_rank_2_ij)|
            assert_eq!(tensor_rank_2_ij, squared_tensor_rank_2_ij)
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
fn transpose_dim_4()
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
fn zero()
{
    TensorRank2::<9>::zero().iter()
    .for_each(|tensor_rank_2_i|
        tensor_rank_2_i.iter()
        .for_each(|tensor_rank_2_ij|
            assert_eq!(tensor_rank_2_ij, &0.0)
        )
    );
}