#[cfg(feature = "math")]
mod public
{
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
    #[test]
    fn tensor_rank_0()
    {
        let _: TensorRank0;
    }
    #[test]
    fn tensor_rank_1()
    {
        let _: TensorRank1<3, 1>;
    }
    #[test]
    fn tensor_rank_1_list()
    {
        let _: TensorRank1List<3, 1, 8>;
    }
    #[test]
    fn tensor_rank_2()
    {
        let _: TensorRank2<3, 1, 1>;
    }
    #[test]
    fn tensor_rank_2_list()
    {
        let _: TensorRank2List<3, 1, 1, 8>;
    }
    #[test]
    fn tensor_rank_3()
    {
        let _: TensorRank3<3, 1, 1, 1>;
    }
    #[test]
    fn tensor_rank_4()
    {
        let _: TensorRank4<3, 1, 1, 1, 1>;
    }
}