#[cfg(feature = "math")]
mod public
{
    use flavio::math::
    {
        TensorRank0,
        TensorRank1,
        TensorRank1List,
        TensorRank1List2D,
        TensorRank2,
        TensorRank2List,
        TensorRank2List2D,
        TensorRank3,
        TensorRank3List,
        TensorRank3List2D,
        TensorRank3List3D,
        TensorRank4,
        TensorRank4List
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
    fn tensor_rank_1_list_2d()
    {
        let _: TensorRank1List2D<3, 1, 8, 8>;
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
    fn tensor_rank_2_list_2d()
    {
        let _: TensorRank2List2D<3, 1, 1, 8, 8>;
    }
    #[test]
    fn tensor_rank_3()
    {
        let _: TensorRank3<3, 1, 1, 1>;
    }
    #[test]
    fn tensor_rank_3_list()
    {
        let _: TensorRank3List<3, 1, 1, 1, 8>;
    }
    #[test]
    fn tensor_rank_3_list_2d()
    {
        let _: TensorRank3List2D<3, 1, 1, 1, 8, 8>;
    }
    #[test]
    fn tensor_rank_3_list_3d()
    {
        let _: TensorRank3List3D<3, 1, 1, 1, 8, 8, 8>;
    }
    #[test]
    fn tensor_rank_4()
    {
        let _: TensorRank4<3, 1, 1, 1, 1>;
    }
    #[test]
    fn tensor_rank_4_list()
    {
        let _: TensorRank4List<3, 1, 1, 1, 1, 8>;
    }
}
