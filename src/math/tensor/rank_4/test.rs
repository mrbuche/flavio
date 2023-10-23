use super::
{
    TensorRank0,
    TensorRank2,
    TensorRank2Traits,
    TensorRank4,
    TensorRank4Traits
};

fn get_array() -> [[[[TensorRank0; 4]; 4]; 4]; 4]
{
    [[[
        [1.0, 0.0, 3.0, 3.0],
        [1.0, 1.0, 3.0, 2.0],
        [1.0, 2.0, 3.0, 3.0],
        [1.0, 1.0, 0.0, 2.0]
        ], [
        [0.0, 3.0, 1.0, 2.0],
        [0.0, 2.0, 3.0, 2.0],
        [0.0, 2.0, 0.0, 1.0],
        [2.0, 3.0, 0.0, 2.0]
        ], [
        [1.0, 0.0, 1.0, 2.0],
        [0.0, 0.0, 1.0, 3.0],
        [1.0, 0.0, 1.0, 3.0],
        [1.0, 3.0, 1.0, 1.0]
        ], [
        [2.0, 1.0, 0.0, 2.0],
        [1.0, 3.0, 2.0, 2.0],
        [2.0, 2.0, 1.0, 1.0],
        [2.0, 1.0, 3.0, 0.0]
    ]], [[
        [0.0, 3.0, 2.0, 3.0],
        [3.0, 1.0, 1.0, 0.0],
        [2.0, 3.0, 2.0, 1.0],
        [2.0, 1.0, 0.0, 0.0]
        ], [
        [3.0, 0.0, 0.0, 2.0],
        [2.0, 1.0, 1.0, 2.0],
        [3.0, 1.0, 3.0, 3.0],
        [0.0, 1.0, 2.0, 3.0]
        ], [
        [2.0, 2.0, 1.0, 0.0],
        [1.0, 0.0, 3.0, 0.0],
        [1.0, 1.0, 3.0, 2.0],
        [3.0, 3.0, 3.0, 1.0]
        ], [
        [2.0, 3.0, 3.0, 3.0],
        [2.0, 2.0, 3.0, 2.0],
        [2.0, 2.0, 0.0, 2.0],
        [2.0, 1.0, 2.0, 3.0]
    ]], [[
        [2.0, 1.0, 2.0, 0.0],
        [2.0, 1.0, 0.0, 2.0],
        [3.0, 2.0, 1.0, 0.0],
        [3.0, 2.0, 2.0, 3.0]
        ], [
        [1.0, 2.0, 3.0, 1.0],
        [0.0, 2.0, 2.0, 3.0],
        [0.0, 1.0, 2.0, 3.0],
        [2.0, 1.0, 0.0, 1.0]
        ], [
        [0.0, 0.0, 0.0, 1.0],
        [1.0, 3.0, 3.0, 0.0],
        [3.0, 0.0, 0.0, 3.0],
        [1.0, 1.0, 3.0, 2.0]
        ], [
        [1.0, 1.0, 2.0, 0.0],
        [0.0, 3.0, 1.0, 1.0],
        [2.0, 1.0, 2.0, 3.0],
        [3.0, 1.0, 3.0, 3.0]
    ]], [[
        [1.0, 2.0, 2.0, 0.0],
        [3.0, 1.0, 0.0, 1.0],
        [0.0, 3.0, 3.0, 1.0],
        [3.0, 3.0, 3.0, 0.0]
        ], [
        [1.0, 1.0, 2.0, 1.0],
        [1.0, 3.0, 2.0, 2.0],
        [0.0, 3.0, 0.0, 3.0],
        [1.0, 0.0, 1.0, 0.0]
        ], [
        [1.0, 3.0, 1.0, 3.0],
        [1.0, 0.0, 1.0, 3.0],
        [0.0, 0.0, 1.0, 3.0],
        [2.0, 2.0, 3.0, 1.0]
        ], [
        [2.0, 1.0, 2.0, 0.0],
        [1.0, 0.0, 2.0, 3.0],
        [3.0, 0.0, 2.0, 2.0],
        [1.0, 2.0, 2.0, 3.0]
    ]]]
}

fn get_tensor_rank_4() -> TensorRank4<4>
{
    TensorRank4::new(get_array())
}

fn get_other_tensor_rank_4() -> TensorRank4<4>
{
    TensorRank4::new([[[
        [3., 0., 0., 1.],
        [0., 3., 2., 1.],
        [2., 1., 3., 3.],
        [0., 1., 3., 3.]],

       [[2., 3., 3., 2.],
        [1., 2., 1., 1.],
        [0., 0., 2., 2.],
        [3., 0., 3., 0.]],

       [[0., 0., 2., 0.],
        [2., 2., 2., 0.],
        [2., 3., 2., 2.],
        [1., 0., 1., 3.]],

       [[3., 1., 3., 0.],
        [2., 1., 3., 3.],
        [3., 1., 2., 1.],
        [0., 2., 2., 0.]]],


      [[[1., 1., 0., 2.],
        [0., 3., 1., 3.],
        [1., 0., 0., 2.],
        [1., 1., 0., 3.]],

       [[3., 2., 2., 3.],
        [3., 1., 1., 1.],
        [3., 1., 0., 3.],
        [2., 0., 3., 0.]],

       [[1., 1., 3., 1.],
        [0., 2., 2., 2.],
        [0., 0., 1., 2.],
        [3., 1., 0., 0.]],

       [[0., 0., 2., 0.],
        [3., 3., 0., 1.],
        [0., 1., 3., 0.],
        [3., 3., 3., 2.]]],


      [[[2., 0., 1., 0.],
        [0., 0., 1., 3.],
        [3., 2., 2., 0.],
        [3., 1., 3., 3.]],

       [[3., 1., 2., 3.],
        [0., 2., 1., 3.],
        [3., 2., 0., 1.],
        [2., 0., 2., 1.]],

       [[2., 3., 2., 1.],
        [2., 1., 0., 2.],
        [1., 2., 2., 0.],
        [3., 2., 3., 0.]],

       [[3., 1., 2., 2.],
        [3., 2., 2., 0.],
        [1., 1., 2., 0.],
        [3., 2., 3., 3.]]],


      [[[3., 2., 0., 3.],
        [3., 3., 1., 0.],
        [3., 0., 2., 3.],
        [2., 2., 1., 0.]],

       [[3., 2., 3., 0.],
        [2., 0., 0., 0.],
        [0., 3., 3., 0.],
        [0., 3., 0., 2.]],

       [[3., 2., 3., 1.],
        [0., 0., 0., 2.],
        [3., 0., 2., 0.],
        [3., 3., 0., 3.]],

       [[0., 1., 1., 0.],
        [0., 1., 0., 3.],
        [3., 0., 2., 3.],
        [2., 3., 2., 1.]
    ]]])
}

fn get_other_tensor_rank_4_add_tensor_rank_4() -> TensorRank4<4>
{
    TensorRank4::new([[[[4., 0., 3., 4.],
        [1., 4., 5., 3.],
        [3., 3., 6., 6.],
        [1., 2., 3., 5.]],

       [[2., 6., 4., 4.],
        [1., 4., 4., 3.],
        [0., 2., 2., 3.],
        [5., 3., 3., 2.]],

       [[1., 0., 3., 2.],
        [2., 2., 3., 3.],
        [3., 3., 3., 5.],
        [2., 3., 2., 4.]],

       [[5., 2., 3., 2.],
        [3., 4., 5., 5.],
        [5., 3., 3., 2.],
        [2., 3., 5., 0.]]],


      [[[1., 4., 2., 5.],
        [3., 4., 2., 3.],
        [3., 3., 2., 3.],
        [3., 2., 0., 3.]],

       [[6., 2., 2., 5.],
        [5., 2., 2., 3.],
        [6., 2., 3., 6.],
        [2., 1., 5., 3.]],

       [[3., 3., 4., 1.],
        [1., 2., 5., 2.],
        [1., 1., 4., 4.],
        [6., 4., 3., 1.]],

       [[2., 3., 5., 3.],
        [5., 5., 3., 3.],
        [2., 3., 3., 2.],
        [5., 4., 5., 5.]]],


      [[[4., 1., 3., 0.],
        [2., 1., 1., 5.],
        [6., 4., 3., 0.],
        [6., 3., 5., 6.]],

       [[4., 3., 5., 4.],
        [0., 4., 3., 6.],
        [3., 3., 2., 4.],
        [4., 1., 2., 2.]],

       [[2., 3., 2., 2.],
        [3., 4., 3., 2.],
        [4., 2., 2., 3.],
        [4., 3., 6., 2.]],

       [[4., 2., 4., 2.],
        [3., 5., 3., 1.],
        [3., 2., 4., 3.],
        [6., 3., 6., 6.]]],


      [[[4., 4., 2., 3.],
        [6., 4., 1., 1.],
        [3., 3., 5., 4.],
        [5., 5., 4., 0.]],

       [[4., 3., 5., 1.],
        [3., 3., 2., 2.],
        [0., 6., 3., 3.],
        [1., 3., 1., 2.]],

       [[4., 5., 4., 4.],
        [1., 0., 1., 5.],
        [3., 0., 3., 3.],
        [5., 5., 3., 4.]],

       [[2., 2., 3., 0.],
        [1., 1., 2., 6.],
        [6., 0., 4., 5.],
        [3., 5., 4., 4.]]]])
}

fn get_other_tensor_rank_4_sub_tensor_rank_4() -> TensorRank4<4>
{
    TensorRank4::new([[[[-2.,  0.,  3.,  2.],
        [ 1., -2.,  1.,  1.],
        [-1.,  1.,  0.,  0.],
        [ 1.,  0., -3., -1.]],

       [[-2.,  0., -2.,  0.],
        [-1.,  0.,  2.,  1.],
        [ 0.,  2., -2., -1.],
        [-1.,  3., -3.,  2.]],

       [[ 1.,  0., -1.,  2.],
        [-2., -2., -1.,  3.],
        [-1., -3., -1.,  1.],
        [ 0.,  3.,  0., -2.]],

       [[-1.,  0., -3.,  2.],
        [-1.,  2., -1., -1.],
        [-1.,  1., -1.,  0.],
        [ 2., -1.,  1.,  0.]]],


      [[[-1.,  2.,  2.,  1.],
        [ 3., -2.,  0., -3.],
        [ 1.,  3.,  2., -1.],
        [ 1.,  0.,  0., -3.]],

       [[ 0., -2., -2., -1.],
        [-1.,  0.,  0.,  1.],
        [ 0.,  0.,  3.,  0.],
        [-2.,  1., -1.,  3.]],

       [[ 1.,  1., -2., -1.],
        [ 1., -2.,  1., -2.],
        [ 1.,  1.,  2.,  0.],
        [ 0.,  2.,  3.,  1.]],

       [[ 2.,  3.,  1.,  3.],
        [-1., -1.,  3.,  1.],
        [ 2.,  1., -3.,  2.],
        [-1., -2., -1.,  1.]]],


      [[[ 0.,  1.,  1.,  0.],
        [ 2.,  1., -1., -1.],
        [ 0.,  0., -1.,  0.],
        [ 0.,  1., -1.,  0.]],

       [[-2.,  1.,  1., -2.],
        [ 0.,  0.,  1.,  0.],
        [-3., -1.,  2.,  2.],
        [ 0.,  1., -2.,  0.]],

       [[-2., -3., -2.,  0.],
        [-1.,  2.,  3., -2.],
        [ 2., -2., -2.,  3.],
        [-2., -1.,  0.,  2.]],

       [[-2.,  0.,  0., -2.],
        [-3.,  1., -1.,  1.],
        [ 1.,  0.,  0.,  3.],
        [ 0., -1.,  0.,  0.]]],


      [[[-2.,  0.,  2., -3.],
        [ 0., -2., -1.,  1.],
        [-3.,  3.,  1., -2.],
        [ 1.,  1.,  2.,  0.]],

       [[-2., -1., -1.,  1.],
        [-1.,  3.,  2.,  2.],
        [ 0.,  0., -3.,  3.],
        [ 1., -3.,  1., -2.]],

       [[-2.,  1., -2.,  2.],
        [ 1.,  0.,  1.,  1.],
        [-3.,  0., -1.,  3.],
        [-1., -1.,  3., -2.]],

       [[ 2.,  0.,  1.,  0.],
        [ 1., -1.,  2.,  0.],
        [ 0.,  0.,  0., -1.],
        [-1., -1.,  0.,  2.]
    ]]])
}

fn get_tensor_rank_2() -> TensorRank2<4>
{
    TensorRank2::new([
        [1.0, 4.0, 6.0, 6.0],
        [1.0, 5.0, 1.0, 0.0],
        [1.0, 3.0, 5.0, 0.0],
        [1.0, 4.0, 6.0, 0.0]
    ])
}

fn get_other_tensor_rank_2() -> TensorRank2<4>
{
    TensorRank2::new([
        [3.0, 2.0, 3.0, 5.0],
        [6.0, 5.0, 2.0, 4.0],
        [4.0, 5.0, 0.0, 4.0],
        [4.0, 4.0, 1.0, 6.0]
    ])
}

#[test]
fn various_contractions()
{
    todo!()
}

#[test]
fn muls_with_tensor_rank_2()
{
    todo!()
}

#[test]
fn as_tensor_rank_2()
{
    todo!()
}

#[test]
fn add_tensor_rank_4_to_self()
{
    (get_tensor_rank_4() + get_other_tensor_rank_4()).iter()
    .zip(get_other_tensor_rank_4_add_tensor_rank_4().iter())
    .for_each(|(tensor_rank_4_i, res_tensor_rank_4_i)|
        tensor_rank_4_i.iter()
        .zip(res_tensor_rank_4_i.iter())
        .for_each(|(tensor_rank_4_ij, res_tensor_rank_4_ij)|
            tensor_rank_4_ij.iter()
            .zip(res_tensor_rank_4_ij.iter())
            .for_each(|(tensor_rank_4_ijk, res_tensor_rank_4_ijk)|
                tensor_rank_4_ijk.iter()
                .zip(res_tensor_rank_4_ijk.iter())
                .for_each(|(tensor_rank_4_ijkl, res_tensor_rank_4_ijkl)|
                    assert_eq!(tensor_rank_4_ijkl, res_tensor_rank_4_ijkl)
                )
            )
        )
    );
}

#[test]
fn add_tensor_rank_4_ref_to_self()
{
    (get_tensor_rank_4() + &get_other_tensor_rank_4()).iter()
    .zip(get_other_tensor_rank_4_add_tensor_rank_4().iter())
    .for_each(|(tensor_rank_4_i, res_tensor_rank_4_i)|
        tensor_rank_4_i.iter()
        .zip(res_tensor_rank_4_i.iter())
        .for_each(|(tensor_rank_4_ij, res_tensor_rank_4_ij)|
            tensor_rank_4_ij.iter()
            .zip(res_tensor_rank_4_ij.iter())
            .for_each(|(tensor_rank_4_ijk, res_tensor_rank_4_ijk)|
                tensor_rank_4_ijk.iter()
                .zip(res_tensor_rank_4_ijk.iter())
                .for_each(|(tensor_rank_4_ijkl, res_tensor_rank_4_ijkl)|
                    assert_eq!(tensor_rank_4_ijkl, res_tensor_rank_4_ijkl)
                )
            )
        )
    );
}

#[test]
fn add_tensor_rank_4_to_self_ref()
{
    (&get_tensor_rank_4() + get_other_tensor_rank_4()).iter()
    .zip(get_other_tensor_rank_4_add_tensor_rank_4().iter())
    .for_each(|(tensor_rank_4_i, res_tensor_rank_4_i)|
        tensor_rank_4_i.iter()
        .zip(res_tensor_rank_4_i.iter())
        .for_each(|(tensor_rank_4_ij, res_tensor_rank_4_ij)|
            tensor_rank_4_ij.iter()
            .zip(res_tensor_rank_4_ij.iter())
            .for_each(|(tensor_rank_4_ijk, res_tensor_rank_4_ijk)|
                tensor_rank_4_ijk.iter()
                .zip(res_tensor_rank_4_ijk.iter())
                .for_each(|(tensor_rank_4_ijkl, res_tensor_rank_4_ijkl)|
                    assert_eq!(tensor_rank_4_ijkl, res_tensor_rank_4_ijkl)
                )
            )
        )
    );
}

#[test]
fn add_assign_tensor_rank_4()
{
    let mut tensor_rank_4 = get_tensor_rank_4();
    tensor_rank_4 += get_other_tensor_rank_4();
    tensor_rank_4.iter()
    .zip(get_other_tensor_rank_4_add_tensor_rank_4().iter())
    .for_each(|(tensor_rank_4_i, res_tensor_rank_4_i)|
        tensor_rank_4_i.iter()
        .zip(res_tensor_rank_4_i.iter())
        .for_each(|(tensor_rank_4_ij, res_tensor_rank_4_ij)|
            tensor_rank_4_ij.iter()
            .zip(res_tensor_rank_4_ij.iter())
            .for_each(|(tensor_rank_4_ijk, res_tensor_rank_4_ijk)|
                tensor_rank_4_ijk.iter()
                .zip(res_tensor_rank_4_ijk.iter())
                .for_each(|(tensor_rank_4_ijkl, res_tensor_rank_4_ijkl)|
                    assert_eq!(tensor_rank_4_ijkl, res_tensor_rank_4_ijkl)
                )
            )
        )
    );
}

#[test]
fn add_assign_tensor_rank_4_ref()
{
    let mut tensor_rank_4 = get_tensor_rank_4();
    tensor_rank_4 += &get_other_tensor_rank_4();
    tensor_rank_4.iter()
    .zip(get_other_tensor_rank_4_add_tensor_rank_4().iter())
    .for_each(|(tensor_rank_4_i, res_tensor_rank_4_i)|
        tensor_rank_4_i.iter()
        .zip(res_tensor_rank_4_i.iter())
        .for_each(|(tensor_rank_4_ij, res_tensor_rank_4_ij)|
            tensor_rank_4_ij.iter()
            .zip(res_tensor_rank_4_ij.iter())
            .for_each(|(tensor_rank_4_ijk, res_tensor_rank_4_ijk)|
                tensor_rank_4_ijk.iter()
                .zip(res_tensor_rank_4_ijk.iter())
                .for_each(|(tensor_rank_4_ijkl, res_tensor_rank_4_ijkl)|
                    assert_eq!(tensor_rank_4_ijkl, res_tensor_rank_4_ijkl)
                )
            )
        )
    );
}

#[test]
fn div_tensor_rank_0_to_self()
{
    (get_tensor_rank_4() / 3.3).iter()
    .zip(get_array().iter())
    .for_each(|(tensor_rank_4_i, array_i)|
        tensor_rank_4_i.iter()
        .zip(array_i.iter())
        .for_each(|(tensor_rank_4_ij, array_ij)|
            tensor_rank_4_ij.iter()
            .zip(array_ij.iter())
            .for_each(|(tensor_rank_4_ijk, array_ijk)|
                tensor_rank_4_ijk.iter()
                .zip(array_ijk.iter())
                .for_each(|(tensor_rank_4_ijkl, array_ijkl)|
                    assert_eq!(tensor_rank_4_ijkl, &(array_ijkl / 3.3))
                )
            )
        )
    );
}

#[test]
fn div_tensor_rank_0_ref_to_self()
{
    (get_tensor_rank_4() / &3.3).iter()
    .zip(get_array().iter())
    .for_each(|(tensor_rank_4_i, array_i)|
        tensor_rank_4_i.iter()
        .zip(array_i.iter())
        .for_each(|(tensor_rank_4_ij, array_ij)|
            tensor_rank_4_ij.iter()
            .zip(array_ij.iter())
            .for_each(|(tensor_rank_4_ijk, array_ijk)|
                tensor_rank_4_ijk.iter()
                .zip(array_ijk.iter())
                .for_each(|(tensor_rank_4_ijkl, array_ijkl)|
                    assert_eq!(tensor_rank_4_ijkl, &(array_ijkl / 3.3))
                )
            )
        )
    );
}

#[test]
fn div_assign_tensor_rank_0()
{
    let mut tensor_rank_4 = get_tensor_rank_4();
    tensor_rank_4 /= 3.3;
    tensor_rank_4.iter()
    .zip(get_array().iter())
    .for_each(|(tensor_rank_4_i, array_i)|
        tensor_rank_4_i.iter()
        .zip(array_i.iter())
        .for_each(|(tensor_rank_4_ij, array_ij)|
            tensor_rank_4_ij.iter()
            .zip(array_ij.iter())
            .for_each(|(tensor_rank_4_ijk, array_ijk)|
                tensor_rank_4_ijk.iter()
                .zip(array_ijk.iter())
                .for_each(|(tensor_rank_4_ijkl, array_ijkl)|
                    assert_eq!(tensor_rank_4_ijkl, &(array_ijkl / 3.3))
                )
            )
        )
    );
}

#[test]
fn div_assign_tensor_rank_0_ref()
{
    let mut tensor_rank_4 = get_tensor_rank_4();
    tensor_rank_4 /= &3.3;
    tensor_rank_4.iter()
    .zip(get_array().iter())
    .for_each(|(tensor_rank_4_i, array_i)|
        tensor_rank_4_i.iter()
        .zip(array_i.iter())
        .for_each(|(tensor_rank_4_ij, array_ij)|
            tensor_rank_4_ij.iter()
            .zip(array_ij.iter())
            .for_each(|(tensor_rank_4_ijk, array_ijk)|
                tensor_rank_4_ijk.iter()
                .zip(array_ijk.iter())
                .for_each(|(tensor_rank_4_ijkl, array_ijkl)|
                    assert_eq!(tensor_rank_4_ijkl, &(array_ijkl / 3.3))
                )
            )
        )
    );
}

#[test]
fn dyad_ij_kl()
{
    let tensor_a = get_tensor_rank_2();
    let tensor_b = get_other_tensor_rank_2();
    TensorRank4::<4>::dyad_ij_kl(&tensor_a, &tensor_b).iter()
    .zip(tensor_a.iter())
    .for_each(|(tensor_rank_4_i, tensor_a_i)|
        tensor_rank_4_i.iter()
        .zip(tensor_a_i.iter())
        .for_each(|(tensor_rank_4_ij, tensor_a_ij)|
            tensor_rank_4_ij.iter()
            .zip(tensor_b.iter())
            .for_each(|(tensor_rank_4_ijk, tensor_b_k)|
                tensor_rank_4_ijk.iter()
                .zip(tensor_b_k.iter())
                .for_each(|(tensor_rank_4_ijkl, tensor_b_kl)|
                    assert_eq!(tensor_rank_4_ijkl, &(tensor_a_ij * tensor_b_kl))
                )
            )
        )
    );
}

#[test]
fn dyad_ik_jl()
{
    let tensor_a = get_tensor_rank_2();
    let tensor_b = get_other_tensor_rank_2();
    TensorRank4::<4>::dyad_ik_jl(&tensor_a, &tensor_b).iter()
    .zip(tensor_a.iter())
    .for_each(|(tensor_rank_4_i, tensor_a_i)|
        tensor_rank_4_i.iter()
        .zip(tensor_b.iter())
        .for_each(|(tensor_rank_4_ij, tensor_b_j)|
            tensor_rank_4_ij.iter()
            .zip(tensor_a_i.iter())
            .for_each(|(tensor_rank_4_ijk, tensor_a_ik)|
                tensor_rank_4_ijk.iter()
                .zip(tensor_b_j.iter())
                .for_each(|(tensor_rank_4_ijkl, tensor_b_jl)|
                    assert_eq!(tensor_rank_4_ijkl, &(tensor_a_ik * tensor_b_jl))
                )
            )
        )
    );
}

#[test]
fn dyad_il_jk()
{
    let tensor_a = get_tensor_rank_2();
    let tensor_b = get_other_tensor_rank_2();
    TensorRank4::<4>::dyad_il_jk(&tensor_a, &tensor_b).iter()
    .zip(tensor_a.iter())
    .for_each(|(tensor_rank_4_i, tensor_a_i)|
        tensor_rank_4_i.iter()
        .zip(tensor_b.iter())
        .for_each(|(tensor_rank_4_ij, tensor_b_j)|
            tensor_rank_4_ij.iter()
            .zip(tensor_b_j.iter())
            .for_each(|(tensor_rank_4_ijk, tensor_b_jk)|
                tensor_rank_4_ijk.iter()
                .zip(tensor_a_i.iter())
                .for_each(|(tensor_rank_4_ijkl, tensor_a_il)|
                    assert_eq!(tensor_rank_4_ijkl, &(tensor_a_il * tensor_b_jk))
                )
            )
        )
    );
}

#[test]
fn dyad_il_kj()
{
    let tensor_a = get_tensor_rank_2();
    let tensor_b = get_other_tensor_rank_2();
    TensorRank4::<4>::dyad_il_kj(&tensor_a, &tensor_b).iter()
    .zip(tensor_a.iter())
    .for_each(|(tensor_rank_4_i, tensor_a_i)|
        tensor_rank_4_i.iter()
        .zip(tensor_b.transpose().iter())
        .for_each(|(tensor_rank_4_ij, tensor_b_j)|
            tensor_rank_4_ij.iter()
            .zip(tensor_b_j.iter())
            .for_each(|(tensor_rank_4_ijk, tensor_b_jk)|
                tensor_rank_4_ijk.iter()
                .zip(tensor_a_i.iter())
                .for_each(|(tensor_rank_4_ijkl, tensor_a_il)|
                    assert_eq!(tensor_rank_4_ijkl, &(tensor_a_il * tensor_b_jk))
                )
            )
        )
    );
}

#[test]
fn from_iter()
{
    let into_iterator = get_tensor_rank_4().0.into_iter();
    let tensor_rank_4 = TensorRank4::<4>::from_iter(get_tensor_rank_4().0.into_iter());
    tensor_rank_4.iter()
    .zip(into_iterator)
    .for_each(|(tensor_rank_4_i, value_i)|
        tensor_rank_4_i.iter()
        .zip(value_i.iter())
        .for_each(|(tensor_rank_4_ij, value_ij)|
            tensor_rank_4_ij.iter()
            .zip(value_ij.iter())
            .for_each(|(tensor_rank_4_ijk, value_ijk)|
                tensor_rank_4_ijk.iter()
                .zip(value_ijk.iter())
                .for_each(|(tensor_rank_4_ijkl, value_ijkl)|
                    assert_eq!(tensor_rank_4_ijkl, value_ijkl)
                )
            )
        )
    );
}

#[test]
fn inverse()
{
    todo!()
}

#[test]
fn iter()
{
    get_tensor_rank_4().iter()
    .zip(get_array().iter())
    .for_each(|(tensor_rank_4_i, array_i)|
        tensor_rank_4_i.iter()
        .zip(array_i.iter())
        .for_each(|(tensor_rank_4_ij, array_ij)|
            tensor_rank_4_ij.iter()
            .zip(array_ij.iter())
            .for_each(|(tensor_rank_4_ijk, array_ijk)|
                assert_eq!(tensor_rank_4_ijk.0, *array_ijk)
            )
        )
    );
}

#[test]
fn iter_mut()
{
    get_tensor_rank_4().iter_mut()
    .zip(get_array().iter_mut())
    .for_each(|(tensor_rank_4_i, array_i)|
        tensor_rank_4_i.iter_mut()
        .zip(array_i.iter_mut())
        .for_each(|(tensor_rank_4_ij, array_ij)|
            tensor_rank_4_ij.iter_mut()
            .zip(array_ij.iter_mut())
            .for_each(|(tensor_rank_4_ijk, array_ijk)|
                assert_eq!(tensor_rank_4_ijk.0, *array_ijk)
            )
        )
    );
}

#[test]
fn mul_tensor_rank_0_to_self()
{
    (get_tensor_rank_4() * 3.3).iter()
    .zip(get_array().iter())
    .for_each(|(tensor_rank_4_i, array_i)|
        tensor_rank_4_i.iter()
        .zip(array_i.iter())
        .for_each(|(tensor_rank_4_ij, array_ij)|
            tensor_rank_4_ij.iter()
            .zip(array_ij.iter())
            .for_each(|(tensor_rank_4_ijk, array_ijk)|
                tensor_rank_4_ijk.iter()
                .zip(array_ijk.iter())
                .for_each(|(tensor_rank_4_ijkl, array_ijkl)|
                    assert_eq!(tensor_rank_4_ijkl, &(array_ijkl * 3.3))
                )
            )
        )
    );
}

#[test]
fn mul_tensor_rank_0_ref_to_self()
{
    (get_tensor_rank_4() * &3.3).iter()
    .zip(get_array().iter())
    .for_each(|(tensor_rank_4_i, array_i)|
        tensor_rank_4_i.iter()
        .zip(array_i.iter())
        .for_each(|(tensor_rank_4_ij, array_ij)|
            tensor_rank_4_ij.iter()
            .zip(array_ij.iter())
            .for_each(|(tensor_rank_4_ijk, array_ijk)|
                tensor_rank_4_ijk.iter()
                .zip(array_ijk.iter())
                .for_each(|(tensor_rank_4_ijkl, array_ijkl)|
                    assert_eq!(tensor_rank_4_ijkl, &(array_ijkl * 3.3))
                )
            )
        )
    );
}

#[test]
fn mul_assign_tensor_rank_0()
{
    let mut tensor_rank_4 = get_tensor_rank_4();
    tensor_rank_4 *= 3.3;
    tensor_rank_4.iter()
    .zip(get_array().iter())
    .for_each(|(tensor_rank_4_i, array_i)|
        tensor_rank_4_i.iter()
        .zip(array_i.iter())
        .for_each(|(tensor_rank_4_ij, array_ij)|
            tensor_rank_4_ij.iter()
            .zip(array_ij.iter())
            .for_each(|(tensor_rank_4_ijk, array_ijk)|
                tensor_rank_4_ijk.iter()
                .zip(array_ijk.iter())
                .for_each(|(tensor_rank_4_ijkl, array_ijkl)|
                    assert_eq!(tensor_rank_4_ijkl, &(array_ijkl * 3.3))
                )
            )
        )
    );
}

#[test]
fn mul_assign_tensor_rank_0_ref()
{
    let mut tensor_rank_4 = get_tensor_rank_4();
    tensor_rank_4 *= &3.3;
    tensor_rank_4.iter()
    .zip(get_array().iter())
    .for_each(|(tensor_rank_4_i, array_i)|
        tensor_rank_4_i.iter()
        .zip(array_i.iter())
        .for_each(|(tensor_rank_4_ij, array_ij)|
            tensor_rank_4_ij.iter()
            .zip(array_ij.iter())
            .for_each(|(tensor_rank_4_ijk, array_ijk)|
                tensor_rank_4_ijk.iter()
                .zip(array_ijk.iter())
                .for_each(|(tensor_rank_4_ijkl, array_ijkl)|
                    assert_eq!(tensor_rank_4_ijkl, &(array_ijkl * 3.3))
                )
            )
        )
    );
}

#[test]
fn new()
{
    get_tensor_rank_4().iter()
    .zip(get_array().iter())
    .for_each(|(get_tensor_rank_4_i, array_i)|
        get_tensor_rank_4_i.iter()
        .zip(array_i.iter())
        .for_each(|(get_tensor_rank_4_ij, array_ij)|
            get_tensor_rank_4_ij.iter()
            .zip(array_ij.iter())
            .for_each(|(get_tensor_rank_4_ijk, array_ijk)|
                get_tensor_rank_4_ijk.iter()
                .zip(array_ijk.iter())
                .for_each(|(get_tensor_rank_4_ijkl, array_ijkl)|
                    assert_eq!(get_tensor_rank_4_ijkl, array_ijkl)
                )
            )
        )
    );
}

#[test]
fn sub_tensor_rank_4_to_self()
{
    (get_tensor_rank_4() - get_other_tensor_rank_4()).iter()
    .zip(get_other_tensor_rank_4_sub_tensor_rank_4().iter())
    .for_each(|(tensor_rank_4_i, res_tensor_rank_4_i)|
        tensor_rank_4_i.iter()
        .zip(res_tensor_rank_4_i.iter())
        .for_each(|(tensor_rank_4_ij, res_tensor_rank_4_ij)|
            tensor_rank_4_ij.iter()
            .zip(res_tensor_rank_4_ij.iter())
            .for_each(|(tensor_rank_4_ijk, res_tensor_rank_4_ijk)|
                tensor_rank_4_ijk.iter()
                .zip(res_tensor_rank_4_ijk.iter())
                .for_each(|(tensor_rank_4_ijkl, res_tensor_rank_4_ijkl)|
                    assert_eq!(tensor_rank_4_ijkl, res_tensor_rank_4_ijkl)
                )
            )
        )
    );
}

#[test]
fn sub_tensor_rank_4_ref_to_self()
{
    (get_tensor_rank_4() - &get_other_tensor_rank_4()).iter()
    .zip(get_other_tensor_rank_4_sub_tensor_rank_4().iter())
    .for_each(|(tensor_rank_4_i, res_tensor_rank_4_i)|
        tensor_rank_4_i.iter()
        .zip(res_tensor_rank_4_i.iter())
        .for_each(|(tensor_rank_4_ij, res_tensor_rank_4_ij)|
            tensor_rank_4_ij.iter()
            .zip(res_tensor_rank_4_ij.iter())
            .for_each(|(tensor_rank_4_ijk, res_tensor_rank_4_ijk)|
                tensor_rank_4_ijk.iter()
                .zip(res_tensor_rank_4_ijk.iter())
                .for_each(|(tensor_rank_4_ijkl, res_tensor_rank_4_ijkl)|
                    assert_eq!(tensor_rank_4_ijkl, res_tensor_rank_4_ijkl)
                )
            )
        )
    );
}

#[test]
fn sub_assign_tensor_rank_4()
{
    let mut tensor_rank_4 = get_tensor_rank_4();
    tensor_rank_4 -= get_other_tensor_rank_4();
    tensor_rank_4.iter()
    .zip(get_other_tensor_rank_4_sub_tensor_rank_4().iter())
    .for_each(|(tensor_rank_4_i, res_tensor_rank_4_i)|
        tensor_rank_4_i.iter()
        .zip(res_tensor_rank_4_i.iter())
        .for_each(|(tensor_rank_4_ij, res_tensor_rank_4_ij)|
            tensor_rank_4_ij.iter()
            .zip(res_tensor_rank_4_ij.iter())
            .for_each(|(tensor_rank_4_ijk, res_tensor_rank_4_ijk)|
                tensor_rank_4_ijk.iter()
                .zip(res_tensor_rank_4_ijk.iter())
                .for_each(|(tensor_rank_4_ijkl, res_tensor_rank_4_ijkl)|
                    assert_eq!(tensor_rank_4_ijkl, res_tensor_rank_4_ijkl)
                )
            )
        )
    );
}

#[test]
fn sub_assign_tensor_rank_4_ref()
{
    let mut tensor_rank_4 = get_tensor_rank_4();
    tensor_rank_4 -= &get_other_tensor_rank_4();
    tensor_rank_4.iter()
    .zip(get_other_tensor_rank_4_sub_tensor_rank_4().iter())
    .for_each(|(tensor_rank_4_i, res_tensor_rank_4_i)|
        tensor_rank_4_i.iter()
        .zip(res_tensor_rank_4_i.iter())
        .for_each(|(tensor_rank_4_ij, res_tensor_rank_4_ij)|
            tensor_rank_4_ij.iter()
            .zip(res_tensor_rank_4_ij.iter())
            .for_each(|(tensor_rank_4_ijk, res_tensor_rank_4_ijk)|
                tensor_rank_4_ijk.iter()
                .zip(res_tensor_rank_4_ijk.iter())
                .for_each(|(tensor_rank_4_ijkl, res_tensor_rank_4_ijkl)|
                    assert_eq!(tensor_rank_4_ijkl, res_tensor_rank_4_ijkl)
                )
            )
        )
    );
}

#[test]
fn zero()
{
    TensorRank4::<4>::zero().iter()
    .for_each(|tensor_rank_4_i|
        tensor_rank_4_i.iter()
        .for_each(|tensor_rank_4_ij|
            tensor_rank_4_ij.iter()
            .for_each(|tensor_rank_4_ijk|
                tensor_rank_4_ijk.iter()
                .for_each(|tensor_rank_4_ijkl|
                    assert_eq!(tensor_rank_4_ijkl, &0.0)
                )
            )
        )
    );
}