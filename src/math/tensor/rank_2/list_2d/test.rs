use super::
{
    TensorRank0,
    TensorRank2List2D,
    TensorRank2List2DTrait
};

fn get_array() -> [[[[TensorRank0; 3]; 3]; 2]; 2]
{
    [[[
        [1.0, 4.0, 6.0],
        [7.0, 2.0, 5.0],
        [9.0, 8.0, 3.0]
    ],[
        [3.0, 2.0, 3.0],
        [6.0, 5.0, 2.0],
        [4.0, 5.0, 0.0]
    ]], [[
        [5.0, 2.0, 9.0],
        [2.0, 4.0, 5.0],
        [1.0, 3.0, 8.0]
    ],[
        [4.0, 3.0, 2.0],
        [2.0, 5.0, 4.0],
        [1.0, 7.0, 1.0]
    ]]]
}

fn get_tensor_rank_2_list_2d() -> TensorRank2List2D<3, 1, 1, 2>
{
    TensorRank2List2D::new(get_array())
}

#[test]
fn from_iter()
{
    let into_iterator = get_tensor_rank_2_list_2d().0.into_iter();
    let tensor_rank_2_list_2d = TensorRank2List2D::<3, 1, 1, 2>::from_iter(get_tensor_rank_2_list_2d().0.into_iter());
    tensor_rank_2_list_2d.iter()
    .zip(into_iterator)
    .for_each(|(tensor_rank_2_list_2d_i, array_i)|
        tensor_rank_2_list_2d_i.iter()
        .zip(array_i.iter())
        .for_each(|(tensor_rank_2_list_2d_ij, array_ij)|
            tensor_rank_2_list_2d_ij.iter()
            .zip(array_ij.iter())
            .for_each(|(tensor_rank_2_list_2d_ij_w, array_ij_w)|
                tensor_rank_2_list_2d_ij_w.iter()
                .zip(array_ij_w.iter())
                .for_each(|(tensor_rank_2_list_2d_ij_ww, array_ij_ww)|
                    assert_eq!(tensor_rank_2_list_2d_ij_ww, array_ij_ww)
                )
            )
        )
    );
}

#[test]
fn iter()
{
    get_tensor_rank_2_list_2d().iter()
    .zip(get_array().iter())
    .for_each(|(tensor_rank_2_list_2d_i, array_i)|
        tensor_rank_2_list_2d_i.iter()
        .zip(array_i.iter())
        .for_each(|(tensor_rank_2_list_2d_ij, array_ij)|
            tensor_rank_2_list_2d_ij.iter()
            .zip(array_ij.iter())
            .for_each(|(tensor_rank_2_list_2d_ij_w, array_ij_w)|
                tensor_rank_2_list_2d_ij_w.iter()
                .zip(array_ij_w.iter())
                .for_each(|(tensor_rank_2_list_2d_ij_ww, array_ij_ww)|
                    assert_eq!(tensor_rank_2_list_2d_ij_ww, array_ij_ww)
                )
            )
        )
    );
}

#[test]
fn iter_mut()
{
    get_tensor_rank_2_list_2d().iter_mut()
    .zip(get_array().iter_mut())
    .for_each(|(tensor_rank_2_list_2d_i, array_i)|
        tensor_rank_2_list_2d_i.iter_mut()
        .zip(array_i.iter_mut())
        .for_each(|(tensor_rank_2_list_2d_ij, array_ij)|
            tensor_rank_2_list_2d_ij.iter_mut()
            .zip(array_ij.iter_mut())
            .for_each(|(tensor_rank_2_list_2d_ij_w, array_ij_w)|
                tensor_rank_2_list_2d_ij_w.iter_mut()
                .zip(array_ij_w.iter_mut())
                .for_each(|(tensor_rank_2_list_2d_ij_ww, array_ij_ww)|
                    assert_eq!(tensor_rank_2_list_2d_ij_ww, array_ij_ww)
                )
            )
        )
    );
}

#[test]
fn new()
{
    get_tensor_rank_2_list_2d().iter()
    .zip(get_array().iter())
    .for_each(|(tensor_rank_2_list_2d_i, array_i)|
        tensor_rank_2_list_2d_i.iter()
        .zip(array_i.iter())
        .for_each(|(tensor_rank_2_list_2d_ij, array_ij)|
            tensor_rank_2_list_2d_ij.iter()
            .zip(array_ij.iter())
            .for_each(|(tensor_rank_2_list_2d_ij_w, array_ij_w)|
                tensor_rank_2_list_2d_ij_w.iter()
                .zip(array_ij_w.iter())
                .for_each(|(tensor_rank_2_list_2d_ij_ww, array_ij_ww)|
                    assert_eq!(tensor_rank_2_list_2d_ij_ww, array_ij_ww)
                )
            )
        )
    );
}

#[test]
fn zero()
{
    TensorRank2List2D::<3, 1, 1, 8>::zero().iter()
    .for_each(|tensor_rank_0_list_i|
        tensor_rank_0_list_i.iter()
        .for_each(|tensor_rank_0_list_ij|
            tensor_rank_0_list_ij.iter()
            .for_each(|tensor_rank_0_list_ij_w|
                tensor_rank_0_list_ij_w.iter()
                .for_each(|tensor_rank_0_list_ij_ww|
                    assert_eq!(tensor_rank_0_list_ij_ww, &0.0)
                )
            )
        )
    );
}
