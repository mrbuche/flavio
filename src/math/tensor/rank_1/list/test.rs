use super::
{
    TensorRank1List,
    TensorRank1ListTraits
};

#[test]
fn add_tensor_rank_1_list_to_self()
{
    todo!()
}

#[test]
fn add_tensor_rank_1_list_ref_to_self()
{
    todo!()
}

#[test]
fn add_tensor_rank_1_list_to_self_ref()
{
    todo!()
}

#[test]
fn add_assign_tensor_rank_1_list()
{
    todo!()
}

#[test]
fn add_assign_tensor_rank_1_list_ref()
{
    todo!()
}

#[test]
fn div_tensor_rank_0_to_self_ref()
{
    todo!()
}

#[test]
fn div_tensor_rank_0_ref_to_self()
{
    todo!()
}

#[test]
fn div_tensor_rank_0_ref_to_self_ref()
{
    todo!()
}

#[test]
fn div_assign_tensor_rank_0()
{
    todo!()
}

#[test]
fn div_assign_tensor_rank_0_ref()
{
    todo!()
}

#[test]
fn from_iter()
{
    todo!()
}

#[test]
fn iter()
{
    todo!()
}

#[test]
fn iter_mut()
{
    todo!()
}

#[test]
fn mul_tensor_rank_0_to_self()
{
    todo!()
}

#[test]
fn mul_tensor_rank_0_to_self_ref()
{
    todo!()
}

#[test]
fn mul_tensor_rank_0_ref_to_self()
{
    todo!()
}

#[test]
fn mul_tensor_rank_0_ref_to_self_ref()
{
    todo!()
}

#[test]
fn mul_assign_tensor_rank_0()
{
    todo!()
}

#[test]
fn mul_assign_tensor_rank_0_ref()
{
    todo!()
}

#[test]
fn mul_tensor_rank_1_list_to_self()
{
    todo!()
}

#[test]
fn mul_tensor_rank_1_list_ref_to_self()
{
    todo!()
}

#[test]
fn mul_tensor_rank_1_list_to_self_ref()
{
    todo!()
}

#[test]
fn mul_tensor_rank_1_list_ref_to_self_ref()
{
    todo!()
}

#[test]
fn new()
{
    todo!()
}

#[test]
fn sub_tensor_rank_1_list_to_self()
{
    todo!()
}

#[test]
fn sub_tensor_rank_1_list_ref_to_self()
{
    todo!()
}

#[test]
fn sub_assign_tensor_rank_1_list()
{
    todo!()
}

#[test]
fn sub_assign_tensor_rank_1_list_ref()
{
    todo!()
}

#[test]
fn zero()
{
    TensorRank1List::<3, 8>::zero().iter()
    .for_each(|tensor_rank_1_entry|
        tensor_rank_1_entry.iter()
        .for_each(|tensor_rank_1_entry_i|
            assert_eq!(tensor_rank_1_entry_i, &0.0)
        )
    );
}
