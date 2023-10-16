use super::
{
    TensorRank0,
    TensorRank2,
    TensorRank2Traits
};

fn get_dummy_array() -> [[TensorRank0; 4]; 4]
{
    [[ 1.0,  2.0,  3.0,  4.0],
     [ 5.0,  6.0,  7.0,  8.0],
     [ 9.0, 10.0, 11.0, 12.0],
     [13.0, 14.0, 15.0, 16.0]]
}

fn get_dummy_tensor_rank_2() -> TensorRank2<4>
{
    TensorRank2::new(get_dummy_array())
}

#[test]
fn todo()
{
    todo!("Make comprehensive tests, and traits for <D> and <3> instead of base impls in order to use for mechanics.");
}

#[test]
fn copy()
{
    let tensor_rank_2_a = get_dummy_tensor_rank_2();
    let tensor_rank_2_b = tensor_rank_2_a;
    tensor_rank_2_a.iter().zip(tensor_rank_2_b.iter()).for_each(|(tensor_rank_2_a_i, tensor_rank_2_b_i)|
        tensor_rank_2_a_i.iter().zip(tensor_rank_2_b_i.iter()).for_each(|(tensor_rank_2_a_ij, tensor_rank_2_b_ij)|
            assert_eq!(tensor_rank_2_a_ij, tensor_rank_2_b_ij)
        )
    );
}

#[test]
fn dyad()
{
    todo!();
}

#[test]
fn from_iter()
{
    todo!();
}

#[test]
fn identity()
{
    TensorRank2::<8>::identity().iter().enumerate().for_each(|(i, tensor_rank_2_i)|
        tensor_rank_2_i.iter().enumerate().for_each(|(j, tensor_rank_2_ij)|
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
fn inverse_2()
{
    todo!();
}

#[test]
fn inverse_3()
{
    todo!();
}

#[test]
fn inverse_4()
{
    todo!();
}

#[test]
fn inverse_5()
{
    todo!();
}

#[test]
fn inverse_9()
{
    todo!();
}

#[test]
fn new()
{
    get_dummy_tensor_rank_2().iter().zip(get_dummy_array().iter()).for_each(|(tensor_rank_2_i, array_i)|
        tensor_rank_2_i.iter().zip(array_i.iter()).for_each(|(tensor_rank_2_ij, array_ij)|
            assert_eq!(tensor_rank_2_ij, array_ij)
        )
    );
}

#[test]
fn norm()
{
    todo!();
}

#[test]
fn transpose()
{
    let tensor_rank_2 = get_dummy_tensor_rank_2();
    let tensor_rank_2_transpose = tensor_rank_2.transpose();
    tensor_rank_2.iter().enumerate().for_each(|(i, tensor_rank_2_i)|
        tensor_rank_2_i.iter().enumerate().for_each(|(j, tensor_rank_2_ij)|
            assert_eq!(tensor_rank_2_ij, &tensor_rank_2_transpose[j][i])
        )
    );
    tensor_rank_2_transpose.iter().enumerate().for_each(|(i, tensor_rank_2_transpose_i)|
        tensor_rank_2_transpose_i.iter().enumerate().for_each(|(j, tensor_rank_2_transpose_ij)|
            assert_eq!(tensor_rank_2_transpose_ij, &tensor_rank_2[j][i])
        )
    );
}

#[test]
fn zero()
{
    TensorRank2::<8>::zero().iter().for_each(|tensor_rank_2_i|
        tensor_rank_2_i.iter().for_each(|tensor_rank_2_ij|
            assert_eq!(tensor_rank_2_ij, &0.0)
        )
    );
}