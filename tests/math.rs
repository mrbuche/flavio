use flavio::math::
{
    TensorRank0,
    TensorRank1,
    TensorRank2,
    TensorRank4
};

// "Each file in the tests directory is a separate crate."

// Can you test that only expected items are visible?

#[test]
fn tensor_rank_0()
{
    let _: TensorRank0;
}

#[test]
fn tensor_rank_1()
{
    let _: TensorRank1<3>;
}

#[test]
fn tensor_rank_2()
{
    let _: TensorRank2<3>;
}

#[test]
fn tensor_rank_4()
{
    let _: TensorRank4<3>;
}
