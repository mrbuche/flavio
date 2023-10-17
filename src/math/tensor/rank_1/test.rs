use super::
{
    TensorRank0,
    TensorRank1,
    TensorRank1Traits
};

fn get_dummy_array() -> [TensorRank0; 4]
{
    [1.0, 2.0, 3.0, 4.0]
}

fn get_dummy_tensor_rank_1() -> TensorRank1<4>
{
    TensorRank1::new(get_dummy_array())
}

#[test]
fn todo()
{
    todo!("Make comprehensive tests, and traits for <D> and <3> instead of base impls in order to use for mechanics.");
}

#[test]
fn add_tensor_rank_1()
{
    todo!("also the reference variants in other tests?");
}

#[test]
fn add_assign_tensor_rank_1()
{
    todo!();
}

#[test]
fn copy()
{
    let tensor_rank_1_a = get_dummy_tensor_rank_1();
    let tensor_rank_1_b = tensor_rank_1_a;
    tensor_rank_1_a.iter().zip(tensor_rank_1_b.iter()).for_each(|(tensor_rank_1_a_i, tensor_rank_1_b_i)|
        assert_eq!(tensor_rank_1_a_i, tensor_rank_1_b_i)
    );
}

#[test]
fn from_iter()
{
    let into_iterator = (0..8).map(|x| x as TensorRank0).into_iter();
    let tensor_rank_1 = TensorRank1::<8>::from_iter(into_iterator.clone());
    tensor_rank_1.iter().zip(into_iterator).for_each(|(tensor_rank_1_i, value_i)|
        assert_eq!(tensor_rank_1_i, &value_i)
    );
}

#[test]
fn new()
{
    get_dummy_tensor_rank_1().iter().zip(get_dummy_array().iter()).for_each(|(tensor_rank_1_i, array_i)|
        assert_eq!(tensor_rank_1_i, array_i)
    );
}

#[test]
fn norm()
{
    assert_eq!(get_dummy_tensor_rank_1().norm(), 30.0);
}

#[test]
fn zero()
{
    TensorRank1::<8>::zero().iter().for_each(|tensor_rank_1_i|
        assert_eq!(tensor_rank_1_i, &0.0)
    );
}
