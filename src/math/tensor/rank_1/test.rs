use super::TensorRank1;

#[test]
fn copy()
{
    let a = TensorRank1([1.0, 2.0, 3.0]);
    let b = a;
    a.iter().zip(b.iter()).for_each(|(a_i, b_i)|
        assert_eq!(a_i, b_i)
    );
}

#[test]
fn from_iter()
{
    let into_iterator = (0..8).map(|x| x as f64).into_iter();
    let tensor_rank_1 = TensorRank1::<8>::from_iter(into_iterator.clone());
    tensor_rank_1.iter().zip(into_iterator).for_each(|(tensor_rank_1_i, value_i)|
        assert_eq!(tensor_rank_1_i, &value_i)
    );
}

#[test]
fn todo()
{
    todo!("Make comprehensive tests, and traits for <D> and <3> instead of base impls in order to use for mechanics.");
}
