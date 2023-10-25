use flavio::math::
{
    TensorRank0,
    TensorRank1,
    TensorRank1List,
    TensorRank2,
    TensorRank2List,
    TensorRank3,
    TensorRank4
};

// "Each file in the tests directory is a separate crate."

// Can you test that only expected items are visible?

// Can you import and test the traits too?

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
fn tensor_rank_1_list()
{
    let _: TensorRank1List<3, 8>;
}

#[test]
fn tensor_rank_2()
{
    let _: TensorRank2<3>;
}

#[test]
fn tensor_rank_2_list()
{
    let _: TensorRank2List<3, 8>;
}

#[test]
fn tensor_rank_3()
{
    let _: TensorRank3<3>;
}

#[test]
fn tensor_rank_4()
{
    let _: TensorRank4<3>;
}
