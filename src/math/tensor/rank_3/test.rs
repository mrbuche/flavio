use super::
{
    TensorRank0,
    TensorRank3,
    TensorRank3Traits
};

fn get_array() -> [[[TensorRank0; 4]; 4]; 4]
{
    [[
        [1.0, 2.0, 1.0, 2.0],
        [0.0, 3.0, 1.0, 3.0],
        [3.0, 1.0, 2.0, 0.0],
        [2.0, 0.0, 3.0, 2.0]
        ], [
        [2.0, 0.0, 3.0, 3.0],
        [0.0, 2.0, 3.0, 1.0],
        [3.0, 3.0, 2.0, 2.0],
        [3.0, 2.0, 3.0, 1.0]
        ], [
        [2.0, 2.0, 0.0, 1.0],
        [1.0, 0.0, 0.0, 1.0],
        [0.0, 2.0, 1.0, 3.0],
        [1.0, 2.0, 3.0, 0.0]
        ], [
        [2.0, 2.0, 1.0, 3.0],
        [0.0, 0.0, 0.0, 3.0],
        [2.0, 2.0, 3.0, 1.0],
        [1.0, 0.0, 3.0, 0.0]
    ]]
}

fn get_tensor_rank_3() -> TensorRank3<4>
{
    TensorRank3::new(get_array())
}

fn get_other_tensor_rank_3() -> TensorRank3<4>
{
    TensorRank3::new([[
        [2.0, 1.0, 1.0, 2.0],
        [3.0, 1.0, 1.0, 2.0],
        [1.0, 3.0, 2.0, 1.0],
        [3.0, 0.0, 0.0, 1.0]
        ], [
        [1.0, 3.0, 3.0, 1.0],
        [3.0, 2.0, 3.0, 2.0],
        [0.0, 3.0, 2.0, 0.0],
        [2.0, 3.0, 3.0, 3.0]
        ], [
        [1.0, 2.0, 1.0, 0.0],
        [2.0, 1.0, 2.0, 3.0],
        [1.0, 3.0, 2.0, 1.0],
        [1.0, 2.0, 3.0, 3.0]
        ], [
        [3.0, 1.0, 1.0, 3.0],
        [0.0, 0.0, 3.0, 2.0],
        [0.0, 2.0, 0.0, 2.0],
        [0.0, 1.0, 2.0, 1.0]
    ]])
}

fn get_other_tensor_rank_3_add_tensor_rank_3() -> TensorRank3<4>
{
    TensorRank3::new([[
        [3.0, 3.0, 2.0, 4.0],
        [3.0, 4.0, 2.0, 5.0],
        [4.0, 4.0, 4.0, 1.0],
        [5.0, 0.0, 3.0, 3.0]
        ], [
        [3.0, 3.0, 6.0, 4.0],
        [3.0, 4.0, 6.0, 3.0],
        [3.0, 6.0, 4.0, 2.0],
        [5.0, 5.0, 6.0, 4.0]
        ], [
        [3.0, 4.0, 1.0, 1.0],
        [3.0, 1.0, 2.0, 4.0],
        [1.0, 5.0, 3.0, 4.0],
        [2.0, 4.0, 6.0, 3.0]
        ], [
        [5.0, 3.0, 2.0, 6.0],
        [0.0, 0.0, 3.0, 5.0],
        [2.0, 4.0, 3.0, 3.0],
        [1.0, 1.0, 5.0, 1.0]
    ]])
}

fn get_other_tensor_rank_3_sub_tensor_rank_3() -> TensorRank3<4>
{
    TensorRank3::new([[
        [-1.0,  1.0,  0.0,  0.0],
        [-3.0,  2.0,  0.0,  1.0],
        [ 2.0, -2.0,  0.0, -1.0],
        [-1.0,  0.0,  3.0,  1.0]
        ], [
        [ 1.0, -3.0,  0.0,  2.0],
        [-3.0,  0.0,  0.0, -1.0],
        [ 3.0,  0.0,  0.0,  2.0],
        [ 1.0, -1.0,  0.0, -2.0]
        ], [
        [ 1.0,  0.0, -1.0,  1.0],
        [-1.0, -1.0, -2.0, -2.0],
        [-1.0, -1.0, -1.0,  2.0],
        [ 0.0,  0.0,  0.0, -3.0]
        ], [
        [-1.0,  1.0,  0.0,  0.0],
        [ 0.0,  0.0, -3.0,  1.0],
        [ 2.0,  0.0,  3.0, -1.0],
        [ 1.0, -1.0,  1.0, -1.0]
    ]])
}

#[test]
fn add_tensor_rank_3_to_self()
{
    (get_tensor_rank_3() + get_other_tensor_rank_3()).iter()
    .zip(get_other_tensor_rank_3_add_tensor_rank_3().iter())
    .for_each(|(tensor_rank_3_i, res_tensor_rank_3_i)|
        tensor_rank_3_i.iter()
        .zip(res_tensor_rank_3_i.iter())
        .for_each(|(tensor_rank_3_ij, res_tensor_rank_3_ij)|
            tensor_rank_3_ij.iter()
            .zip(res_tensor_rank_3_ij.iter())
            .for_each(|(tensor_rank_3_ijk, res_tensor_rank_3_ijk)|
                assert_eq!(tensor_rank_3_ijk, res_tensor_rank_3_ijk)
            )
        )
    );
}

#[test]
fn add_tensor_rank_3_ref_to_self()
{
    (get_tensor_rank_3() + &get_other_tensor_rank_3()).iter()
    .zip(get_other_tensor_rank_3_add_tensor_rank_3().iter())
    .for_each(|(tensor_rank_3_i, res_tensor_rank_3_i)|
        tensor_rank_3_i.iter()
        .zip(res_tensor_rank_3_i.iter())
        .for_each(|(tensor_rank_3_ij, res_tensor_rank_3_ij)|
            tensor_rank_3_ij.iter()
            .zip(res_tensor_rank_3_ij.iter())
            .for_each(|(tensor_rank_3_ijk, res_tensor_rank_3_ijk)|
                assert_eq!(tensor_rank_3_ijk, res_tensor_rank_3_ijk)
            )
        )
    );
}

#[test]
fn add_tensor_rank_3_to_self_ref()
{
    (&get_tensor_rank_3() + get_other_tensor_rank_3()).iter()
    .zip(get_other_tensor_rank_3_add_tensor_rank_3().iter())
    .for_each(|(tensor_rank_3_i, res_tensor_rank_3_i)|
        tensor_rank_3_i.iter()
        .zip(res_tensor_rank_3_i.iter())
        .for_each(|(tensor_rank_3_ij, res_tensor_rank_3_ij)|
            tensor_rank_3_ij.iter()
            .zip(res_tensor_rank_3_ij.iter())
            .for_each(|(tensor_rank_3_ijk, res_tensor_rank_3_ijk)|
                assert_eq!(tensor_rank_3_ijk, res_tensor_rank_3_ijk)
            )
        )
    );
}

#[test]
fn add_assign_tensor_rank_3()
{
    let mut tensor_rank_3 = get_tensor_rank_3();
    tensor_rank_3 += get_other_tensor_rank_3();
    tensor_rank_3.iter()
    .zip(get_other_tensor_rank_3_add_tensor_rank_3().iter())
    .for_each(|(tensor_rank_3_i, res_tensor_rank_3_i)|
        tensor_rank_3_i.iter()
        .zip(res_tensor_rank_3_i.iter())
        .for_each(|(tensor_rank_3_ij, res_tensor_rank_3_ij)|
            tensor_rank_3_ij.iter()
            .zip(res_tensor_rank_3_ij.iter())
            .for_each(|(tensor_rank_3_ijk, res_tensor_rank_3_ijk)|
                assert_eq!(tensor_rank_3_ijk, res_tensor_rank_3_ijk)
            )
        )
    );
}

#[test]
fn add_assign_tensor_rank_3_ref()
{
    let mut tensor_rank_3 = get_tensor_rank_3();
    tensor_rank_3 += &get_other_tensor_rank_3();
    tensor_rank_3.iter()
    .zip(get_other_tensor_rank_3_add_tensor_rank_3().iter())
    .for_each(|(tensor_rank_3_i, res_tensor_rank_3_i)|
        tensor_rank_3_i.iter()
        .zip(res_tensor_rank_3_i.iter())
        .for_each(|(tensor_rank_3_ij, res_tensor_rank_3_ij)|
            tensor_rank_3_ij.iter()
            .zip(res_tensor_rank_3_ij.iter())
            .for_each(|(tensor_rank_3_ijk, res_tensor_rank_3_ijk)|
                assert_eq!(tensor_rank_3_ijk, res_tensor_rank_3_ijk)
            )
        )
    );
}

#[test]
fn div_tensor_rank_0_to_self()
{
    (get_tensor_rank_3() / 3.3).iter()
    .zip(get_array().iter())
    .for_each(|(tensor_rank_3_i, array_i)|
        tensor_rank_3_i.iter()
        .zip(array_i.iter())
        .for_each(|(tensor_rank_3_ij, array_ij)|
            tensor_rank_3_ij.iter()
            .zip(array_ij.iter())
            .for_each(|(tensor_rank_3_ijk, array_ijk)|
                assert_eq!(tensor_rank_3_ijk, &(array_ijk / 3.3))
            )
        )
    );
}

#[test]
fn div_tensor_rank_0_ref_to_self()
{
    (get_tensor_rank_3() / &3.3).iter()
    .zip(get_array().iter())
    .for_each(|(tensor_rank_3_i, array_i)|
        tensor_rank_3_i.iter()
        .zip(array_i.iter())
        .for_each(|(tensor_rank_3_ij, array_ij)|
            tensor_rank_3_ij.iter()
            .zip(array_ij.iter())
            .for_each(|(tensor_rank_3_ijk, array_ijk)|
                assert_eq!(tensor_rank_3_ijk, &(array_ijk / 3.3))
            )
        )
    );
}

#[test]
fn div_assign_tensor_rank_0()
{
    let mut tensor_rank_3 = get_tensor_rank_3();
    tensor_rank_3 /= 3.3;
    tensor_rank_3.iter()
    .zip(get_array().iter())
    .for_each(|(tensor_rank_3_i, array_i)|
        tensor_rank_3_i.iter()
        .zip(array_i.iter())
        .for_each(|(tensor_rank_3_ij, array_ij)|
            tensor_rank_3_ij.iter()
            .zip(array_ij.iter())
            .for_each(|(tensor_rank_3_ijk, array_ijk)|
                assert_eq!(tensor_rank_3_ijk, &(array_ijk / 3.3))
            )
        )
    );
}

#[test]
fn div_assign_tensor_rank_0_ref()
{
    let mut tensor_rank_3 = get_tensor_rank_3();
    tensor_rank_3 /= &3.3;
    tensor_rank_3.iter()
    .zip(get_array().iter())
    .for_each(|(tensor_rank_3_i, array_i)|
        tensor_rank_3_i.iter()
        .zip(array_i.iter())
        .for_each(|(tensor_rank_3_ij, array_ij)|
            tensor_rank_3_ij.iter()
            .zip(array_ij.iter())
            .for_each(|(tensor_rank_3_ijk, array_ijk)|
                assert_eq!(tensor_rank_3_ijk, &(array_ijk / 3.3))
            )
        )
    );
}

#[test]
fn from_iter()
{
    let into_iterator = get_tensor_rank_3().0.into_iter();
    let tensor_rank_3 = TensorRank3::<4>::from_iter(get_tensor_rank_3().0.into_iter());
    tensor_rank_3.iter()
    .zip(into_iterator)
    .for_each(|(tensor_rank_3_i, value_i)|
        tensor_rank_3_i.iter()
        .zip(value_i.iter())
        .for_each(|(tensor_rank_3_ij, value_ij)|
            tensor_rank_3_ij.iter()
            .zip(value_ij.iter())
            .for_each(|(tensor_rank_3_ijk, value_ijk)|
                assert_eq!(tensor_rank_3_ijk, value_ijk)
            )
        )
    );
}

#[test]
fn iter()
{
    get_tensor_rank_3().iter()
    .zip(get_array().iter())
    .for_each(|(tensor_rank_3_i, array_i)|
        tensor_rank_3_i.iter()
        .zip(array_i.iter())
        .for_each(|(tensor_rank_3_ij, array_ij)|
            assert_eq!(tensor_rank_3_ij.0, *array_ij)
        )
    );
}

#[test]
fn iter_mut()
{
    get_tensor_rank_3().iter_mut()
    .zip(get_array().iter_mut())
    .for_each(|(tensor_rank_3_i, array_i)|
        tensor_rank_3_i.iter_mut()
        .zip(array_i.iter_mut())
        .for_each(|(tensor_rank_3_ij, array_ij)|
            assert_eq!(tensor_rank_3_ij.0, *array_ij)
        )
    );
}

#[test]
fn mul_tensor_rank_0_to_self()
{
    (get_tensor_rank_3() * 3.3).iter()
    .zip(get_array().iter())
    .for_each(|(tensor_rank_3_i, array_i)|
        tensor_rank_3_i.iter()
        .zip(array_i.iter())
        .for_each(|(tensor_rank_3_ij, array_ij)|
            tensor_rank_3_ij.iter()
            .zip(array_ij.iter())
            .for_each(|(tensor_rank_3_ijk, array_ijk)|
                assert_eq!(tensor_rank_3_ijk, &(array_ijk * 3.3))
            )
        )
    );
}

#[test]
fn mul_tensor_rank_0_ref_to_self()
{
    (get_tensor_rank_3() * &3.3).iter()
    .zip(get_array().iter())
    .for_each(|(tensor_rank_3_i, array_i)|
        tensor_rank_3_i.iter()
        .zip(array_i.iter())
        .for_each(|(tensor_rank_3_ij, array_ij)|
            tensor_rank_3_ij.iter()
            .zip(array_ij.iter())
            .for_each(|(tensor_rank_3_ijk, array_ijk)|
                assert_eq!(tensor_rank_3_ijk, &(array_ijk * 3.3))
            )
        )
    );
}

#[test]
fn mul_assign_tensor_rank_0()
{
    let mut tensor_rank_3 = get_tensor_rank_3();
    tensor_rank_3*= 3.3;
    tensor_rank_3.iter()
    .zip(get_array().iter())
    .for_each(|(tensor_rank_3_i, array_i)|
        tensor_rank_3_i.iter()
        .zip(array_i.iter())
        .for_each(|(tensor_rank_3_ij, array_ij)|
            tensor_rank_3_ij.iter()
            .zip(array_ij.iter())
            .for_each(|(tensor_rank_3_ijk, array_ijk)|
                assert_eq!(tensor_rank_3_ijk, &(array_ijk * 3.3))
            )
        )
    );
}

#[test]
fn mul_assign_tensor_rank_0_ref()
{
    let mut tensor_rank_3 = get_tensor_rank_3();
    tensor_rank_3*= &3.3;
    tensor_rank_3.iter()
    .zip(get_array().iter())
    .for_each(|(tensor_rank_3_i, array_i)|
        tensor_rank_3_i.iter()
        .zip(array_i.iter())
        .for_each(|(tensor_rank_3_ij, array_ij)|
            tensor_rank_3_ij.iter()
            .zip(array_ij.iter())
            .for_each(|(tensor_rank_3_ijk, array_ijk)|
                assert_eq!(tensor_rank_3_ijk, &(array_ijk * 3.3))
            )
        )
    );
}

#[test]
fn new()
{
    get_tensor_rank_3().iter()
    .zip(get_array().iter())
    .for_each(|(tensor_rank_3_i, array_i)|
        tensor_rank_3_i.iter()
        .zip(array_i.iter())
        .for_each(|(tensor_rank_3_ij, array_ij)|
            tensor_rank_3_ij.iter()
            .zip(array_ij.iter())
            .for_each(|(tensor_rank_3_ijk, array_ijk)|
                assert_eq!(tensor_rank_3_ijk, array_ijk)
            )
        )
    );
}

#[test]
fn sub_tensor_rank_3_to_self()
{
    (get_tensor_rank_3() - get_other_tensor_rank_3()).iter()
    .zip(get_other_tensor_rank_3_sub_tensor_rank_3().iter())
    .for_each(|(tensor_rank_3_i, res_tensor_rank_3_i)|
        tensor_rank_3_i.iter()
        .zip(res_tensor_rank_3_i.iter())
        .for_each(|(tensor_rank_3_ij, res_tensor_rank_3_ij)|
            tensor_rank_3_ij.iter()
            .zip(res_tensor_rank_3_ij.iter())
            .for_each(|(tensor_rank_3_ijk, res_tensor_rank_3_ijk)|
                assert_eq!(tensor_rank_3_ijk, res_tensor_rank_3_ijk)
            )
        )
    );
}

#[test]
fn sub_tensor_rank_3_ref_to_self()
{
    (get_tensor_rank_3() - &get_other_tensor_rank_3()).iter()
    .zip(get_other_tensor_rank_3_sub_tensor_rank_3().iter())
    .for_each(|(tensor_rank_3_i, res_tensor_rank_3_i)|
        tensor_rank_3_i.iter()
        .zip(res_tensor_rank_3_i.iter())
        .for_each(|(tensor_rank_3_ij, res_tensor_rank_3_ij)|
            tensor_rank_3_ij.iter()
            .zip(res_tensor_rank_3_ij.iter())
            .for_each(|(tensor_rank_3_ijk, res_tensor_rank_3_ijk)|
                assert_eq!(tensor_rank_3_ijk, res_tensor_rank_3_ijk)
            )
        )
    );
}

#[test]
fn sub_assign_tensor_rank_3()
{
    let mut tensor_rank_3 = get_tensor_rank_3();
    tensor_rank_3 -= get_other_tensor_rank_3();
    tensor_rank_3.iter()
    .zip(get_other_tensor_rank_3_sub_tensor_rank_3().iter())
    .for_each(|(tensor_rank_3_i, res_tensor_rank_3_i)|
        tensor_rank_3_i.iter()
        .zip(res_tensor_rank_3_i.iter())
        .for_each(|(tensor_rank_3_ij, res_tensor_rank_3_ij)|
            tensor_rank_3_ij.iter()
            .zip(res_tensor_rank_3_ij.iter())
            .for_each(|(tensor_rank_3_ijk, res_tensor_rank_3_ijk)|
                assert_eq!(tensor_rank_3_ijk, res_tensor_rank_3_ijk)
            )
        )
    );
}

#[test]
fn sub_assign_tensor_rank_3_ref()
{
    let mut tensor_rank_3 = get_tensor_rank_3();
    tensor_rank_3 -= &get_other_tensor_rank_3();
    tensor_rank_3.iter()
    .zip(get_other_tensor_rank_3_sub_tensor_rank_3().iter())
    .for_each(|(tensor_rank_3_i, res_tensor_rank_3_i)|
        tensor_rank_3_i.iter()
        .zip(res_tensor_rank_3_i.iter())
        .for_each(|(tensor_rank_3_ij, res_tensor_rank_3_ij)|
            tensor_rank_3_ij.iter()
            .zip(res_tensor_rank_3_ij.iter())
            .for_each(|(tensor_rank_3_ijk, res_tensor_rank_3_ijk)|
                assert_eq!(tensor_rank_3_ijk, res_tensor_rank_3_ijk)
            )
        )
    );
}

#[test]
fn zero()
{
    TensorRank3::<4>::zero().iter()
    .for_each(|tensor_rank_3_i|
        tensor_rank_3_i.iter()
        .for_each(|tensor_rank_3_ij|
            tensor_rank_3_ij.iter()
            .for_each(|tensor_rank_3_ijk|
                assert_eq!(tensor_rank_3_ijk, &0.0)
            )
        )
    );
}