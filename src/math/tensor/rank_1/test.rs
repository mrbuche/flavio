use super::
{
    TensorRank0,
    TensorRank1,
    TensorRank1Traits
};

fn get_array() -> [TensorRank0; 4]
{
    [1.0, 2.0, 3.0, 4.0]
}

fn get_tensor_rank_1() -> TensorRank1<4>
{
    TensorRank1::new(get_array())
}

fn get_tensor_rank_1_a() -> TensorRank1<4>
{
    TensorRank1::new([5.0, 7.0, 6.0, 8.0])
}

fn get_tensor_rank_1_add_tensor_rank_1_a() -> TensorRank1<4>
{
    TensorRank1::new([6.0, 9.0, 9.0, 12.0])
}

fn get_tensor_rank_1_mul_tensor_rank_1_a() -> TensorRank0
{
    69.0
}

fn get_tensor_rank_1_sub_tensor_rank_1_a() -> TensorRank1<4>
{
    TensorRank1::new([-4.0, -5.0, -3.0, -4.0])
}

#[test]
fn todo()
{
    todo!("Make comprehensive tests, and traits for <D> and <3> instead of base impls in order to use for mechanics.");
}

#[test]
fn add_tensor_rank_1_to_self()
{
    (get_tensor_rank_1() + get_tensor_rank_1_a()).iter()
    .zip(get_tensor_rank_1_add_tensor_rank_1_a().iter())
    .for_each(|(tensor_rank_1_i, add_tensor_rank_1_i)|
        assert_eq!(tensor_rank_1_i, add_tensor_rank_1_i)
    );
}

#[test]
fn add_tensor_rank_1_ref_to_self()
{
    (get_tensor_rank_1() + &get_tensor_rank_1_a()).iter()
    .zip(get_tensor_rank_1_add_tensor_rank_1_a().iter())
    .for_each(|(tensor_rank_1_i, add_tensor_rank_1_i)|
        assert_eq!(tensor_rank_1_i, add_tensor_rank_1_i)
    );
}

#[test]
fn add_tensor_rank_1_to_self_ref()
{
    (&get_tensor_rank_1() + get_tensor_rank_1_a()).iter()
    .zip(get_tensor_rank_1_add_tensor_rank_1_a().iter())
    .for_each(|(tensor_rank_1_i, add_tensor_rank_1_i)|
        assert_eq!(tensor_rank_1_i, add_tensor_rank_1_i)
    );
}

#[test]
fn add_assign_tensor_rank_1()
{
    let mut tensor_rank_1 = get_tensor_rank_1();
    tensor_rank_1 += get_tensor_rank_1_a();
    tensor_rank_1.iter()
    .zip(get_tensor_rank_1_add_tensor_rank_1_a().iter())
    .for_each(|(tensor_rank_1_i, add_tensor_rank_1_i)|
        assert_eq!(tensor_rank_1_i, add_tensor_rank_1_i)
    );
}

#[test]
fn add_assign_tensor_rank_1_ref()
{
    let mut tensor_rank_1 = get_tensor_rank_1();
    tensor_rank_1 += &get_tensor_rank_1_a();
    tensor_rank_1.iter()
    .zip(get_tensor_rank_1_add_tensor_rank_1_a().iter())
    .for_each(|(tensor_rank_1_i, add_tensor_rank_1_i)|
        assert_eq!(tensor_rank_1_i, add_tensor_rank_1_i)
    );
}

#[test]
fn copy()
{
    let tensor_rank_1_a = get_tensor_rank_1();
    let tensor_rank_1_b = tensor_rank_1_a;
    tensor_rank_1_a.iter()
    .zip(tensor_rank_1_b.iter())
    .for_each(|(tensor_rank_1_a_i, tensor_rank_1_b_i)|
        assert_eq!(tensor_rank_1_a_i, tensor_rank_1_b_i)
    );
}

#[test]
fn from_iter()
{
    let into_iterator = (0..8).map(|x| x as TensorRank0).into_iter();
    let tensor_rank_1 = TensorRank1::<8>::from_iter(into_iterator.clone());
    tensor_rank_1.iter()
    .zip(into_iterator)
    .for_each(|(tensor_rank_1_i, value_i)|
        assert_eq!(tensor_rank_1_i, &value_i)
    );
}

#[test]
fn mul_tensor_rank_1_to_self()
{
    assert_eq!(get_tensor_rank_1() * get_tensor_rank_1_a(), get_tensor_rank_1_mul_tensor_rank_1_a())
}

#[test]
fn mul_tensor_rank_1_ref_to_self()
{
    assert_eq!(get_tensor_rank_1() * &get_tensor_rank_1_a(), get_tensor_rank_1_mul_tensor_rank_1_a())
}

#[test]
fn mul_tensor_rank_1_to_self_ref()
{
    assert_eq!(&get_tensor_rank_1() * get_tensor_rank_1_a(), get_tensor_rank_1_mul_tensor_rank_1_a())
}

#[test]
fn mul_tensor_rank_1_ref_to_self_ref()
{
    assert_eq!(&get_tensor_rank_1() * &get_tensor_rank_1_a(), get_tensor_rank_1_mul_tensor_rank_1_a())
}

#[test]
fn new()
{
    get_tensor_rank_1().iter()
    .zip(get_array().iter())
    .for_each(|(tensor_rank_1_i, array_i)|
        assert_eq!(tensor_rank_1_i, array_i)
    );
}

#[test]
fn norm()
{
    assert_eq!(get_tensor_rank_1().norm(), 30.0);
}

#[test]
fn sub_tensor_rank_1_to_self()
{
    (get_tensor_rank_1() - get_tensor_rank_1_a()).iter()
    .zip(get_tensor_rank_1_sub_tensor_rank_1_a().iter())
    .for_each(|(tensor_rank_1_i, sub_tensor_rank_1_i)|
        assert_eq!(tensor_rank_1_i, sub_tensor_rank_1_i)
    );
}

#[test]
fn sub_tensor_rank_1_ref_to_self()
{
    (get_tensor_rank_1() - &get_tensor_rank_1_a()).iter()
    .zip(get_tensor_rank_1_sub_tensor_rank_1_a().iter())
    .for_each(|(tensor_rank_1_i, sub_tensor_rank_1_i)|
        assert_eq!(tensor_rank_1_i, sub_tensor_rank_1_i)
    );
}

#[test]
fn sub_tensor_rank_1_to_self_ref()
{
    (&get_tensor_rank_1() - get_tensor_rank_1_a()).iter()
    .zip(get_tensor_rank_1_sub_tensor_rank_1_a().iter())
    .for_each(|(tensor_rank_1_i, sub_tensor_rank_1_i)|
        assert_eq!(tensor_rank_1_i, sub_tensor_rank_1_i)
    );
}

#[test]
fn sub_assign_tensor_rank_1()
{
    let mut tensor_rank_1 = get_tensor_rank_1();
    tensor_rank_1 -= get_tensor_rank_1_a();
    tensor_rank_1.iter()
    .zip(get_tensor_rank_1_sub_tensor_rank_1_a().iter())
    .for_each(|(tensor_rank_1_i, sub_tensor_rank_1_i)|
        assert_eq!(tensor_rank_1_i, sub_tensor_rank_1_i)
    );
}

#[test]
fn sub_assign_tensor_rank_1_ref()
{
    let mut tensor_rank_1 = get_tensor_rank_1();
    tensor_rank_1 -= &get_tensor_rank_1_a();
    tensor_rank_1.iter()
    .zip(get_tensor_rank_1_sub_tensor_rank_1_a().iter())
    .for_each(|(tensor_rank_1_i, sub_tensor_rank_1_i)|
        assert_eq!(tensor_rank_1_i, sub_tensor_rank_1_i)
    );
}

#[test]
fn zero()
{
    TensorRank1::<8>::zero().iter()
    .for_each(|tensor_rank_1_i|
        assert_eq!(tensor_rank_1_i, &0.0)
    );
}
