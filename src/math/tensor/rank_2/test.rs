use super::{
    super::test::{assert_eq, assert_eq_within_tols, TensorError, TestError},
    Tensor, TensorRank0, TensorRank1, TensorRank1List, TensorRank2, TensorRank2List2D, TensorRank4,
};
use crate::{ABS_TOL, REL_TOL};
use std::cmp::Ordering;

fn get_array_dim_2() -> [[TensorRank0; 2]; 2] {
    [[1.0, 2.0], [3.0, 4.0]]
}

fn get_array_dim_3() -> [[TensorRank0; 3]; 3] {
    [[1.0, 4.0, 6.0], [7.0, 2.0, 5.0], [9.0, 8.0, 3.0]]
}

fn get_array_dim_4() -> [[TensorRank0; 4]; 4] {
    [
        [1.0, 4.0, 6.0, 6.0],
        [1.0, 5.0, 1.0, 0.0],
        [1.0, 3.0, 5.0, 0.0],
        [1.0, 4.0, 6.0, 0.0],
    ]
}

fn get_array_dim_5() -> [[TensorRank0; 5]; 5] {
    [
        [0.71090581, 0.46870553, 0.61471311, 0.32058057, 0.24220277],
        [0.81896269, 0.94126734, 0.3150277 , 0.83762531, 0.41263593],
        [0.53827188, 0.25011113, 0.5942029 , 0.10833879, 0.4730989 ],
        [0.9740403 , 0.9954912 , 0.49567574, 0.36914927, 0.73993237],
        [0.24507233, 0.64323682, 0.79238416, 0.97693838, 0.99396137]
    ]
}

fn get_array_dim_9() -> [[TensorRank0; 9]; 9] {
    [
        [2.0, 2.0, 4.0, 0.0, 0.0, 1.0, 1.0, 3.0, 3.0],
        [0.0, 3.0, 1.0, 0.0, 0.0, 1.0, 4.0, 2.0, 1.0],
        [3.0, 0.0, 1.0, 2.0, 0.0, 3.0, 4.0, 4.0, 2.0],
        [4.0, 4.0, 0.0, 2.0, 1.0, 1.0, 0.0, 0.0, 4.0],
        [0.0, 1.0, 0.0, 1.0, 1.0, 3.0, 0.0, 1.0, 1.0],
        [4.0, 2.0, 3.0, 4.0, 2.0, 4.0, 3.0, 0.0, 4.0],
        [1.0, 3.0, 2.0, 0.0, 0.0, 0.0, 2.0, 4.0, 2.0],
        [2.0, 2.0, 2.0, 4.0, 1.0, 2.0, 4.0, 2.0, 2.0],
        [1.0, 2.0, 3.0, 4.0, 0.0, 1.0, 4.0, 2.0, 1.0],
    ]
}

fn get_tensor_rank_1_a() -> TensorRank1<4, 1> {
    TensorRank1::new([1.0, 2.0, 3.0, 4.0])
}

fn get_tensor_rank_1_b() -> TensorRank1<4, 1> {
    TensorRank1::new([5.0, 7.0, 6.0, 8.0])
}

fn get_tensor_rank_2_dim_2() -> TensorRank2<2, 1, 1> {
    TensorRank2::new(get_array_dim_2())
}

fn get_tensor_rank_2_dim_3() -> TensorRank2<3, 1, 1> {
    TensorRank2::new(get_array_dim_3())
}

fn get_tensor_rank_2_dim_4() -> TensorRank2<4, 1, 1> {
    TensorRank2::new(get_array_dim_4())
}

fn get_tensor_rank_2_dim_5() -> TensorRank2<5, 1, 1> {
    TensorRank2::new(get_array_dim_5())
}

fn get_tensor_rank_2_dim_9() -> TensorRank2<9, 1, 1> {
    TensorRank2::new(get_array_dim_9())
}

fn get_other_tensor_rank_2_dim_2() -> TensorRank2<2, 1, 1> {
    TensorRank2::new([[5.0, 6.0], [7.0, 8.0]])
}

fn get_other_tensor_rank_2_dim_3() -> TensorRank2<3, 1, 1> {
    TensorRank2::new([[3.0, 2.0, 3.0], [6.0, 5.0, 2.0], [4.0, 5.0, 0.0]])
}

fn get_other_tensor_rank_2_dim_4() -> TensorRank2<4, 1, 1> {
    TensorRank2::new([
        [3.0, 2.0, 3.0, 5.0],
        [6.0, 5.0, 2.0, 4.0],
        [4.0, 5.0, 0.0, 4.0],
        [4.0, 4.0, 1.0, 6.0],
    ])
}

fn get_diagonal_tensor_rank_2_dim_4() -> TensorRank2<4, 1, 1> {
    TensorRank2::new([
        [3.0, 0.0, 0.0, 0.0],
        [0.0, 5.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 6.0],
    ])
}

fn get_other_tensor_rank_2_dim_9() -> TensorRank2<9, 1, 1> {
    TensorRank2::new([
        [0.0, 4.0, 2.0, 0.0, 1.0, 4.0, 2.0, 4.0, 1.0],
        [1.0, 2.0, 2.0, 1.0, 0.0, 3.0, 0.0, 2.0, 0.0],
        [3.0, 0.0, 2.0, 3.0, 3.0, 0.0, 0.0, 0.0, 2.0],
        [2.0, 3.0, 0.0, 0.0, 1.0, 3.0, 3.0, 4.0, 2.0],
        [0.0, 4.0, 1.0, 3.0, 1.0, 1.0, 1.0, 2.0, 1.0],
        [1.0, 3.0, 0.0, 3.0, 3.0, 2.0, 1.0, 3.0, 4.0],
        [0.0, 0.0, 0.0, 1.0, 0.0, 3.0, 1.0, 3.0, 4.0],
        [2.0, 0.0, 4.0, 3.0, 1.0, 2.0, 0.0, 3.0, 4.0],
        [4.0, 2.0, 0.0, 0.0, 4.0, 0.0, 4.0, 2.0, 2.0],
    ])
}

fn get_other_tensor_rank_2_mul_tensor_rank_1_dim_4() -> TensorRank1<4, 1> {
    TensorRank1::new([51.0, 14.0, 22.0, 27.0])
}

fn get_other_tensor_rank_2_add_tensor_rank_2_dim_4() -> TensorRank2<4, 1, 1> {
    TensorRank2::new([
        [4.0, 6.0, 9.0, 11.0],
        [7.0, 10.0, 3.0, 4.0],
        [5.0, 8.0, 5.0, 4.0],
        [5.0, 8.0, 7.0, 6.0],
    ])
}

fn get_other_tensor_rank_2_sub_tensor_rank_2_dim_4() -> TensorRank2<4, 1, 1> {
    TensorRank2::new([
        [-2.0, 2.0, 3.0, 1.0],
        [-5.0, 0.0, -1.0, -4.0],
        [-3.0, -2.0, 5.0, -4.0],
        [-3.0, 0.0, 5.0, -6.0],
    ])
}

fn get_other_tensor_rank_2_mul_tensor_rank_2_dim_4() -> TensorRank2<4, 1, 1> {
    TensorRank2::new([
        [75.0, 76.0, 17.0, 81.0],
        [37.0, 32.0, 13.0, 29.0],
        [41.0, 42.0, 9.0, 37.0],
        [51.0, 52.0, 11.0, 45.0],
    ])
}

fn get_tensor_rank_1_list() -> TensorRank1List<3, 1, 8> {
    TensorRank1List::new([
        [5.0, 0.0, 0.0],
        [5.0, 5.0, 6.0],
        [3.0, 1.0, 4.0],
        [3.0, 4.0, 2.0],
        [1.0, 0.0, 3.0],
        [1.0, 3.0, 1.0],
        [1.0, 6.0, 0.0],
        [1.0, 1.0, 1.0],
    ])
}

fn get_tensor_rank_2_list_2d() -> TensorRank2List2D<3, 1, 1, 2, 2> {
    TensorRank2List2D::new([
        [
            [[1.0, 4.0, 6.0], [7.0, 2.0, 5.0], [9.0, 8.0, 3.0]],
            [[3.0, 2.0, 3.0], [6.0, 5.0, 2.0], [4.0, 5.0, 0.0]],
        ],
        [
            [[5.0, 2.0, 9.0], [2.0, 4.0, 5.0], [1.0, 3.0, 8.0]],
            [[4.0, 3.0, 2.0], [2.0, 5.0, 4.0], [1.0, 7.0, 1.0]],
        ],
    ])
}

fn get_tensor_rank_2_mul_tensor_rank_2_list_2d() -> TensorRank2List2D<3, 1, 1, 2, 2> {
    TensorRank2List2D::new([
        [
            [[83.0, 60.0, 44.0], [66.0, 72.0, 67.0], [92.0, 76.0, 103.0]],
            [[51.0, 52.0, 11.0], [53.0, 49.0, 25.0], [87.0, 73.0, 43.0]],
        ],
        [
            [[19.0, 36.0, 77.0], [44.0, 37.0, 113.0], [64.0, 59.0, 145.0]],
            [[18.0, 65.0, 24.0], [37.0, 66.0, 27.0], [55.0, 88.0, 53.0]],
        ],
    ])
}

fn get_tensor_rank_4_dim_3() -> TensorRank4<3, 1, 1, 2, 3> {
    TensorRank4::new([
        [
            [[7.0, 3.0, 7.0], [3.0, 2.0, 7.0], [9.0, 8.0, 4.0]],
            [[1.0, 10.0, 7.0], [0.0, 3.0, 3.0], [4.0, 8.0, 8.0]],
            [[0.0, 1.0, 7.0], [1.0, 2.0, 9.0], [3.0, 5.0, 4.0]],
        ],
        [
            [[2.0, 1.0, 8.0], [6.0, 2.0, 6.0], [4.0, 6.0, 2.0]],
            [[7.0, 7.0, 8.0], [8.0, 4.0, 4.0], [10.0, 9.0, 9.0]],
            [[3.0, 3.0, 3.0], [1.0, 4.0, 3.0], [10.0, 9.0, 5.0]],
        ],
        [
            [[9.0, 5.0, 1.0], [7.0, 9.0, 9.0], [5.0, 9.0, 10.0]],
            [[5.0, 9.0, 0.0], [4.0, 5.0, 7.0], [5.0, 4.0, 7.0]],
            [[1.0, 2.0, 7.0], [8.0, 2.0, 6.0], [2.0, 7.0, 5.0]],
        ],
    ])
}

fn get_tensor_rank_4_dim_4() -> TensorRank4<4, 1, 1, 2, 3> {
    TensorRank4::new([[[[ 0.36504245,  0.25958712, -0.49094648,  0.23852932],
        [ 0.34197679,  0.42071977,  0.33126981, -0.1893317 ],
        [ 0.20671662,  0.08686201, -0.00903859, -0.04902273],
        [-0.03723894,  0.47279997, -0.04718978, -0.01714976]],

       [[ 0.48294834, -0.20821728, -0.22394587,  0.2821748 ],
        [-0.45891843, -0.38735857,  0.19808051, -0.21364514],
        [-0.03610032, -0.42469122, -0.0582858 ,  0.43604045],
        [-0.14018656,  0.18051686,  0.26486812,  0.42422411]],

       [[ 0.33261466,  0.08027591, -0.36318408, -0.22953091],
        [-0.36292561,  0.00613202,  0.12206715, -0.23223834],
        [ 0.19674955, -0.19307904, -0.09875444,  0.2027698 ],
        [ 0.37122124, -0.3540528 ,  0.43354083, -0.11340162]],

       [[ 0.24938961, -0.03151568, -0.07170848,  0.01027823],
        [ 0.39285943, -0.27067197, -0.39616398,  0.25810698],
        [ 0.10250197, -0.47695138, -0.42770545,  0.18593659],
        [ 0.18049076, -0.23975184, -0.22326055,  0.340205  ]]],


      [[[-0.45636765, -0.1518267 ,  0.22664406, -0.29907472],
        [ 0.05522072, -0.29506855, -0.22186704, -0.06410391],
        [ 0.12644432,  0.25416653,  0.44748232, -0.08556985],
        [-0.41869128,  0.03754867,  0.17352907,  0.33009593]],

       [[ 0.0585868 ,  0.40857121, -0.3631779 ,  0.13189339],
        [ 0.025262  , -0.40104261,  0.30517377, -0.21726244],
        [ 0.15097828, -0.22304707,  0.33906195,  0.06631058],
        [-0.14161196, -0.19071875, -0.11441335,  0.18769477]],

       [[ 0.3785788 ,  0.48882586, -0.28561109,  0.39631061],
        [ 0.12883195,  0.40472449,  0.38549137, -0.32514343],
        [ 0.15621689,  0.04922377, -0.34870569, -0.26869119],
        [ 0.10058429,  0.07847799, -0.33330574, -0.12381485]],

       [[ 0.21935125, -0.48869409, -0.43990478,  0.13679917],
        [-0.2491771 ,  0.20400182, -0.33036077, -0.40581175],
        [-0.36796176,  0.25286535,  0.36997294,  0.23885918],
        [-0.20720237, -0.38815007, -0.28260394, -0.00389615]]],


      [[[-0.46518456, -0.23029132, -0.24008545,  0.05420127],
        [ 0.28251842,  0.45665523, -0.40047365,  0.18006298],
        [-0.10406522,  0.30923804,  0.35670474, -0.08450248],
        [ 0.01228638, -0.23805422, -0.20934271, -0.07603387]],

       [[ 0.3274761 ,  0.48652302,  0.1590766 ,  0.10642943],
        [-0.31065525, -0.05133154, -0.43861131, -0.07373563],
        [-0.35599307, -0.20935542,  0.29142214,  0.03419599],
        [ 0.13985499, -0.38136637, -0.00683118,  0.47251816]],

       [[-0.16099946, -0.31994012, -0.47474626,  0.0565323 ],
        [ 0.43676038,  0.03845753, -0.26803846, -0.07011935],
        [ 0.40163798,  0.22278568, -0.19230069, -0.2618809 ],
        [ 0.41605806, -0.487748  , -0.26802637, -0.37403794]],

       [[-0.37472943, -0.01402907, -0.45868484,  0.20504375],
        [-0.07664344, -0.24600269, -0.29901455, -0.30379262],
        [-0.15276443, -0.26685095,  0.26490452,  0.3547848 ],
        [-0.35543205,  0.1092247 ,  0.31776305,  0.43800556]]],


      [[[ 0.48331134, -0.39470774, -0.24443201, -0.31482083],
        [-0.01702343,  0.42210837,  0.25275603, -0.4231595 ],
        [-0.30498086,  0.11660921, -0.12087058,  0.09060479],
        [ 0.36313984,  0.23668824, -0.24189811,  0.24523653]],

       [[ 0.21694912,  0.21071948, -0.35453846,  0.00204652],
        [ 0.27753646,  0.26592748, -0.03682113, -0.36669156],
        [-0.15768757,  0.4930917 , -0.30483748,  0.16793988],
        [ 0.28913016, -0.19149777,  0.13052259, -0.19102942]],

       [[ 0.24206344, -0.09209475, -0.43666025,  0.22532292],
        [-0.0671709 , -0.32964919,  0.00830651,  0.11496252],
        [-0.25396396,  0.28746468, -0.02420304,  0.13051961],
        [ 0.35704736,  0.38527746,  0.42540446,  0.21107368]],

       [[ 0.16301161,  0.17698186, -0.3456682 , -0.35423248],
        [ 0.16426966,  0.14919343,  0.4246134 , -0.43421738],
        [ 0.16996703,  0.05641823, -0.38721139, -0.29874544],
        [-0.07514367,  0.33395239,  0.32808475, -0.3366101 ]]]])
}

fn get_tensor_rank_2_div_tensor_rank_4_dim_3() -> TensorRank2<3, 2, 3> {
    TensorRank2::new([
        [-0.8591023283605275, 0.5463144610682097, 0.48148464803521684],
        [0.14461826142457423, 2.8819091589827597, 0.3555608669979796],
        [
            0.29609312727618836,
            -0.4778620587076813,
            -1.3810401169942013,
        ],
    ])
}

fn get_tensor_rank_2_div_tensor_rank_4_dim_4() -> TensorRank2<4, 2, 3> {
    TensorRank2::new([
        [ 4.881756346219273, -4.109013352742245, -5.185049749853222, 2.6476641850025713],
        [ 7.1683318540707255, 0.5386125990539328, 9.852359208151153, 6.6815534237210965],
        [-0.3432022275100972, 3.4629202930908183, -1.8381771381459546, -6.697343615969693],
        [-1.0719872798261108, -15.717993837871234, 10.863632415009409, 13.987579848929212]])
}

#[test]
fn add_tensor_rank_2_to_self() -> Result<(), TestError> {
    assert_eq(
        &(get_tensor_rank_2_dim_4() + get_other_tensor_rank_2_dim_4()),
        &get_other_tensor_rank_2_add_tensor_rank_2_dim_4(),
    )
}

#[test]
fn add_tensor_rank_2_ref_to_self() -> Result<(), TestError> {
    assert_eq(
        &(get_tensor_rank_2_dim_4() + &get_other_tensor_rank_2_dim_4()),
        &get_other_tensor_rank_2_add_tensor_rank_2_dim_4(),
    )
}

#[test]
fn add_tensor_rank_2_to_self_ref() -> Result<(), TestError> {
    assert_eq(
        &(&get_tensor_rank_2_dim_4() + get_other_tensor_rank_2_dim_4()),
        &get_other_tensor_rank_2_add_tensor_rank_2_dim_4(),
    )
}

#[test]
fn add_assign_tensor_rank_2() -> Result<(), TestError> {
    let mut tensor_rank_2 = get_tensor_rank_2_dim_4();
    tensor_rank_2 += get_other_tensor_rank_2_dim_4();
    assert_eq(
        &tensor_rank_2,
        &get_other_tensor_rank_2_add_tensor_rank_2_dim_4(),
    )
}

#[test]
fn add_assign_tensor_rank_2_ref() -> Result<(), TestError> {
    let mut tensor_rank_2 = get_tensor_rank_2_dim_4();
    tensor_rank_2 += &get_other_tensor_rank_2_dim_4();
    assert_eq(
        &tensor_rank_2,
        &get_other_tensor_rank_2_add_tensor_rank_2_dim_4(),
    )
}

#[test]
fn as_array_dim_2() {
    assert_eq!(get_tensor_rank_2_dim_2().as_array(), get_array_dim_2())
}

#[test]
fn as_array_dim_3() {
    assert_eq!(get_tensor_rank_2_dim_3().as_array(), get_array_dim_3())
}

#[test]
fn as_array_dim_4() {
    assert_eq!(get_tensor_rank_2_dim_4().as_array(), get_array_dim_4())
}

#[test]
fn as_array_dim_9() {
    assert_eq!(get_tensor_rank_2_dim_9().as_array(), get_array_dim_9())
}

#[test]
fn div_tensor_rank_4_to_self_dim_3() -> Result<(), TestError> {
    assert_eq(
        &(get_tensor_rank_2_dim_3() / get_tensor_rank_4_dim_3()),
        &get_tensor_rank_2_div_tensor_rank_4_dim_3(),
    )
}

#[test]
fn div_tensor_rank_4_to_self_dim_4() -> Result<(), TestError> {
    assert_eq(
        &(get_tensor_rank_2_dim_4() / get_tensor_rank_4_dim_4()),
        &get_tensor_rank_2_div_tensor_rank_4_dim_4(),
    )
}

#[test]
fn div_tensor_rank_0_to_self() -> Result<(), TestError> {
    (get_tensor_rank_2_dim_4() / 3.3)
        .iter()
        .zip(get_array_dim_4().iter())
        .try_for_each(|(tensor_rank_2_i, array_i)| {
            tensor_rank_2_i.iter().zip(array_i.iter()).try_for_each(
                |(tensor_rank_2_ij, array_ij)| assert_eq(tensor_rank_2_ij, &(array_ij / 3.3)),
            )
        })?;
    Ok(())
}

#[test]
fn div_tensor_rank_0_to_self_ref() -> Result<(), TestError> {
    (&get_tensor_rank_2_dim_4() / 3.3)
        .iter()
        .zip(get_array_dim_4().iter())
        .try_for_each(|(tensor_rank_2_i, array_i)| {
            tensor_rank_2_i.iter().zip(array_i.iter()).try_for_each(
                |(tensor_rank_2_ij, array_ij)| assert_eq(tensor_rank_2_ij, &(array_ij / 3.3)),
            )
        })?;
    Ok(())
}

#[test]
#[allow(clippy::op_ref)]
fn div_tensor_rank_0_ref_to_self() -> Result<(), TestError> {
    (get_tensor_rank_2_dim_4() / &3.3)
        .iter()
        .zip(get_array_dim_4().iter())
        .try_for_each(|(tensor_rank_2_i, array_i)| {
            tensor_rank_2_i.iter().zip(array_i.iter()).try_for_each(
                |(tensor_rank_2_ij, array_ij)| assert_eq(tensor_rank_2_ij, &(array_ij / 3.3)),
            )
        })?;
    Ok(())
}

#[test]
#[allow(clippy::op_ref)]
fn div_tensor_rank_0_ref_to_self_ref() -> Result<(), TestError> {
    (&get_tensor_rank_2_dim_4() / &3.3)
        .iter()
        .zip(get_array_dim_4().iter())
        .try_for_each(|(tensor_rank_2_i, array_i)| {
            tensor_rank_2_i.iter().zip(array_i.iter()).try_for_each(
                |(tensor_rank_2_ij, array_ij)| assert_eq(tensor_rank_2_ij, &(array_ij / 3.3)),
            )
        })?;
    Ok(())
}

#[test]
fn div_assign_tensor_rank_0() -> Result<(), TestError> {
    let mut tensor_rank_2 = get_tensor_rank_2_dim_4();
    tensor_rank_2 /= 3.3;
    tensor_rank_2
        .iter()
        .zip(get_array_dim_4().iter())
        .try_for_each(|(tensor_rank_2_i, array_i)| {
            tensor_rank_2_i.iter().zip(array_i.iter()).try_for_each(
                |(tensor_rank_2_ij, array_ij)| assert_eq(tensor_rank_2_ij, &(array_ij / 3.3)),
            )
        })?;
    Ok(())
}

#[test]
fn div_assign_tensor_rank_0_ref() -> Result<(), TestError> {
    let mut tensor_rank_2 = get_tensor_rank_2_dim_4();
    tensor_rank_2 /= &3.3;
    tensor_rank_2
        .iter()
        .zip(get_array_dim_4().iter())
        .try_for_each(|(tensor_rank_2_i, array_i)| {
            tensor_rank_2_i.iter().zip(array_i.iter()).try_for_each(
                |(tensor_rank_2_ij, array_ij)| assert_eq(tensor_rank_2_ij, &(array_ij / 3.3)),
            )
        })?;
    Ok(())
}

#[test]
fn determinant_dim_2() -> Result<(), TestError> {
    assert_eq(&get_tensor_rank_2_dim_2().determinant(), &-2.0)
}

#[test]
fn determinant_dim_3() -> Result<(), TestError> {
    assert_eq(&get_tensor_rank_2_dim_3().determinant(), &290.0)
}

#[test]
fn determinant_dim_4() -> Result<(), TestError> {
    assert_eq_within_tols(&get_tensor_rank_2_dim_4().determinant(), &36.0)
}

#[test]
fn determinant_dim_9() -> Result<(), TestError> {
    assert_eq_within_tols(&get_tensor_rank_2_dim_9().determinant(), &2398.0)
}

#[test]
fn deviatoric_dim_2() -> Result<(), TestError> {
    let tensor_rank_2 = get_tensor_rank_2_dim_2();
    let trace = tensor_rank_2.trace();
    let deviatoric_tensor_rank_2 = tensor_rank_2.deviatoric();
    assert_eq(&deviatoric_tensor_rank_2.trace(), &0.0)?;
    assert_eq(
        &deviatoric_tensor_rank_2,
        &(tensor_rank_2 - TensorRank2::identity() * (trace / 2.0)),
    )
}

#[test]
fn deviatoric_dim_3() -> Result<(), TestError> {
    let tensor_rank_2 = get_tensor_rank_2_dim_3();
    let trace = tensor_rank_2.trace();
    let deviatoric_tensor_rank_2 = tensor_rank_2.deviatoric();
    assert_eq(&deviatoric_tensor_rank_2.trace(), &0.0)?;
    assert_eq(
        &deviatoric_tensor_rank_2,
        &(tensor_rank_2 - TensorRank2::identity() * (trace / 3.0)),
    )
}

#[test]
fn deviatoric_dim_4() -> Result<(), TestError> {
    let tensor_rank_2 = get_tensor_rank_2_dim_4();
    let trace = tensor_rank_2.trace();
    let deviatoric_tensor_rank_2 = tensor_rank_2.deviatoric();
    assert_eq(&deviatoric_tensor_rank_2.trace(), &0.0)?;
    assert_eq(
        &deviatoric_tensor_rank_2,
        &(tensor_rank_2 - TensorRank2::identity() * (trace / 4.0)),
    )
}

#[test]
fn deviatoric_dim_9() -> Result<(), TestError> {
    let tensor_rank_2 = get_tensor_rank_2_dim_9();
    let trace = tensor_rank_2.trace();
    let deviatoric_tensor_rank_2 = tensor_rank_2.deviatoric();
    assert_eq(&deviatoric_tensor_rank_2.trace(), &0.0)?;
    assert_eq(
        &deviatoric_tensor_rank_2,
        &(tensor_rank_2 - TensorRank2::identity() * (trace / 9.0)),
    )
}

#[test]
fn deviatoric_and_trace_dim_2() -> Result<(), TestError> {
    let tensor_rank_2 = get_tensor_rank_2_dim_2();
    let (deviatoric, trace) = tensor_rank_2.deviatoric_and_trace();
    assert_eq(&tensor_rank_2.trace(), &trace)?;
    assert_eq(&tensor_rank_2.deviatoric(), &deviatoric)
}

#[test]
fn deviatoric_and_trace_dim_3() -> Result<(), TestError> {
    let tensor_rank_2 = get_tensor_rank_2_dim_3();
    let (deviatoric, trace) = tensor_rank_2.deviatoric_and_trace();
    assert_eq(&tensor_rank_2.trace(), &trace)?;
    assert_eq(&tensor_rank_2.deviatoric(), &deviatoric)
}

#[test]
fn deviatoric_and_trace_dim_4() -> Result<(), TestError> {
    let tensor_rank_2 = get_tensor_rank_2_dim_4();
    let (deviatoric, trace) = tensor_rank_2.deviatoric_and_trace();
    assert_eq(&tensor_rank_2.trace(), &trace)?;
    assert_eq(&tensor_rank_2.deviatoric(), &deviatoric)
}

#[test]
fn deviatoric_and_trace_dim_9() -> Result<(), TestError> {
    let tensor_rank_2 = get_tensor_rank_2_dim_9();
    let (deviatoric, trace) = tensor_rank_2.deviatoric_and_trace();
    assert_eq(&tensor_rank_2.trace(), &trace)?;
    assert_eq(&tensor_rank_2.deviatoric(), &deviatoric)
}

#[test]
fn dyad() {
    let tensor_rank_1_a = get_tensor_rank_1_a();
    let tensor_rank_1_b = get_tensor_rank_1_b();
    let tensor_rank_2 = TensorRank2::dyad(&tensor_rank_1_a, &tensor_rank_1_b);
    tensor_rank_2.iter().zip(tensor_rank_1_a.iter()).for_each(
        |(tensor_rank_2_i, tensor_rank_1_a_i)| {
            tensor_rank_2_i.iter().zip(tensor_rank_1_b.iter()).for_each(
                |(tensor_rank_2_ij, tensor_rank_1_b_j)| {
                    assert_eq!(tensor_rank_2_ij, &(tensor_rank_1_a_i * tensor_rank_1_b_j))
                },
            )
        },
    );
}

#[test]
fn error() {
    let a = get_tensor_rank_1_a();
    let b = get_tensor_rank_1_b();
    assert_eq!(a.error(&a, &ABS_TOL, &REL_TOL), None);
    assert_eq!(a.error(&b, &ABS_TOL, &REL_TOL), Some(4));
}

#[test]
fn from_iter() {
    let into_iterator = get_tensor_rank_2_dim_4().0.into_iter();
    let tensor_rank_2 = TensorRank2::<4, 1, 1>::from_iter(get_tensor_rank_2_dim_4().0);
    tensor_rank_2
        .iter()
        .zip(into_iterator)
        .for_each(|(tensor_rank_2_i, value_i)| {
            tensor_rank_2_i
                .iter()
                .zip(value_i.iter())
                .for_each(|(tensor_rank_2_ij, value_ij)| assert_eq!(tensor_rank_2_ij, value_ij))
        });
}

#[test]
fn full_contraction_dim_2() -> Result<(), TestError> {
    assert_eq_within_tols(
        &get_tensor_rank_2_dim_2().full_contraction(&get_other_tensor_rank_2_dim_2()),
        &70.0,
    )
}

#[test]
fn full_contraction_dim_3() -> Result<(), TestError> {
    assert_eq_within_tols(
        &get_tensor_rank_2_dim_3().full_contraction(&get_other_tensor_rank_2_dim_3()),
        &167.0,
    )
}

#[test]
fn full_contraction_dim_4() -> Result<(), TestError> {
    assert_eq_within_tols(
        &get_tensor_rank_2_dim_4().full_contraction(&get_other_tensor_rank_2_dim_4()),
        &137.0,
    )
}

#[test]
fn full_contraction_dim_9() -> Result<(), TestError> {
    assert_eq_within_tols(
        &get_tensor_rank_2_dim_9().full_contraction(&get_other_tensor_rank_2_dim_9()),
        &269.0,
    )
}

#[test]
fn identity() {
    TensorRank2::<9, 1, 1>::identity()
        .iter()
        .enumerate()
        .for_each(|(i, tensor_rank_2_i)| {
            tensor_rank_2_i
                .iter()
                .enumerate()
                .for_each(|(j, tensor_rank_2_ij)| {
                    if i == j {
                        assert_eq!(tensor_rank_2_ij, &1.0)
                    } else {
                        assert_eq!(tensor_rank_2_ij, &0.0)
                    }
                })
        });
}

#[test]
fn inverse_dim_2() -> Result<(), TestError> {
    assert_eq_within_tols(
        &(get_tensor_rank_2_dim_2() * get_tensor_rank_2_dim_2().inverse()),
        &TensorRank2::identity(),
    )
}

#[test]
fn inverse_dim_3() -> Result<(), TestError> {
    assert_eq_within_tols(
        &(get_tensor_rank_2_dim_3() * get_tensor_rank_2_dim_3().inverse()),
        &TensorRank2::identity(),
    )
}

#[test]
fn inverse_dim_4() -> Result<(), TestError> {
    assert_eq_within_tols(
        &(get_tensor_rank_2_dim_4() * get_tensor_rank_2_dim_4().inverse()),
        &TensorRank2::identity(),
    )
}

#[test]
fn inverse_dim_5() -> Result<(), TestError> {
    assert_eq_within_tols(
        &(get_tensor_rank_2_dim_5() * get_tensor_rank_2_dim_5().inverse()),
        &TensorRank2::identity(),
    )
}

#[test]
fn inverse_dim_9() -> Result<(), TestError> {
    assert_eq_within_tols(
        &(get_tensor_rank_2_dim_9() * get_tensor_rank_2_dim_9().inverse()),
        &TensorRank2::identity(),
    )
}

#[test]
fn inverse_and_determinant_dim_2() -> Result<(), TestError> {
    let tensor_rank_2 = get_tensor_rank_2_dim_2();
    let (inverse, determinant) = tensor_rank_2.inverse_and_determinant();
    assert_eq(&determinant, &tensor_rank_2.determinant())?;
    assert_eq(&inverse, &tensor_rank_2.inverse())
}

#[test]
fn inverse_and_determinant_dim_3() -> Result<(), TestError> {
    let tensor_rank_2 = get_tensor_rank_2_dim_3();
    let (inverse, determinant) = tensor_rank_2.inverse_and_determinant();
    assert_eq(&determinant, &tensor_rank_2.determinant())?;
    assert_eq(&inverse, &tensor_rank_2.inverse())
}

#[test]
fn inverse_and_determinant_dim_4() -> Result<(), TestError> {
    let tensor_rank_2 = get_tensor_rank_2_dim_4();
    let (inverse, determinant) = tensor_rank_2.inverse_and_determinant();
    assert_eq(&determinant, &tensor_rank_2.determinant())?;
    assert_eq(&inverse, &tensor_rank_2.inverse())
}

#[test]
fn inverse_transpose_dim_2() -> Result<(), TestError> {
    assert_eq_within_tols(
        &(get_tensor_rank_2_dim_2().transpose() * get_tensor_rank_2_dim_2().inverse_transpose()),
        &TensorRank2::identity(),
    )
}

#[test]
fn inverse_transpose_dim_3() -> Result<(), TestError> {
    assert_eq_within_tols(
        &(get_tensor_rank_2_dim_3().transpose() * get_tensor_rank_2_dim_3().inverse_transpose()),
        &TensorRank2::identity(),
    )
}

#[test]
fn inverse_transpose_dim_4() -> Result<(), TestError> {
    assert_eq_within_tols(
        &(get_tensor_rank_2_dim_4().transpose() * get_tensor_rank_2_dim_4().inverse_transpose()),
        &TensorRank2::identity(),
    )
}

#[test]
fn inverse_transpose_9() -> Result<(), TestError> {
    assert_eq_within_tols(
        &(get_tensor_rank_2_dim_9().transpose() * get_tensor_rank_2_dim_9().inverse_transpose()),
        &TensorRank2::identity(),
    )
}

#[test]
fn inverse_transpose_and_determinant_dim_2() -> Result<(), TestError> {
    let tensor_rank_2 = get_tensor_rank_2_dim_2();
    let (inverse_transpose, determinant) = tensor_rank_2.inverse_transpose_and_determinant();
    assert_eq(&determinant, &tensor_rank_2.determinant())?;
    assert_eq(&inverse_transpose, &tensor_rank_2.inverse_transpose())
}

#[test]
fn inverse_transpose_and_determinant_dim_3() -> Result<(), TestError> {
    let tensor_rank_2 = get_tensor_rank_2_dim_3();
    let (inverse_transpose, determinant) = tensor_rank_2.inverse_transpose_and_determinant();
    assert_eq(&determinant, &tensor_rank_2.determinant())?;
    assert_eq(&inverse_transpose, &tensor_rank_2.inverse_transpose())
}

#[test]
fn inverse_transpose_and_determinant_dim_4() -> Result<(), TestError> {
    let tensor_rank_2 = get_tensor_rank_2_dim_4();
    let (inverse_transpose, determinant) = tensor_rank_2.inverse_transpose_and_determinant();
    assert_eq(&determinant, &tensor_rank_2.determinant())?;
    assert_eq(&inverse_transpose, &tensor_rank_2.inverse_transpose())
}

#[test]
fn is_diagonal() {
    assert!(get_diagonal_tensor_rank_2_dim_4().is_diagonal())
}

#[test]
fn is_not_diagonal() {
    assert!(!get_other_tensor_rank_2_dim_4().is_diagonal())
}

#[test]
fn is_diagonal_identity() {
    assert!(TensorRank2::<3, 0, 0>::identity().is_diagonal())
}

#[test]
fn is_diagonal_zero() {
    assert!(TensorRank2::<4, 1, 1>::zero().is_diagonal())
}

#[test]
fn is_identity_dim_3() {
    assert!(TensorRank2::<3, 0, 0>::identity().is_identity())
}

#[test]
fn is_not_identity_dim_3() {
    assert!(!TensorRank2::<3, 0, 0>::zero().is_identity())
}

#[test]
fn is_identity_dim_4() {
    assert!(TensorRank2::<4, 1, 1>::identity().is_identity())
}

#[test]
fn is_not_identity_dim_4() {
    assert!(!get_diagonal_tensor_rank_2_dim_4().is_identity())
}

#[test]
fn is_zero_dim_3() {
    assert!(TensorRank2::<3, 0, 0>::zero().is_zero())
}

#[test]
fn is_not_zero_dim_3() {
    assert!(!TensorRank2::<3, 0, 0>::identity().is_zero())
}

#[test]
fn is_zero_dim_4() {
    assert!(TensorRank2::<4, 1, 1>::zero().is_zero())
}

#[test]
fn is_not_zero_dim_4() {
    assert!(!get_other_tensor_rank_2_dim_4().is_zero())
}

#[test]
fn iter() {
    get_tensor_rank_2_dim_4()
        .iter()
        .zip(get_array_dim_4().iter())
        .for_each(|(tensor_rank_2_i, array_i)| {
            tensor_rank_2_i
                .iter()
                .zip(array_i.iter())
                .for_each(|(tensor_rank_2_ij, array_ij)| assert_eq!(tensor_rank_2_ij, array_ij))
        });
}

#[test]
fn iter_mut() {
    get_tensor_rank_2_dim_4()
        .iter_mut()
        .zip(get_array_dim_4().iter_mut())
        .for_each(|(tensor_rank_2_i, array_i)| {
            tensor_rank_2_i
                .iter_mut()
                .zip(array_i.iter_mut())
                .for_each(|(tensor_rank_2_ij, array_ij)| assert_eq!(tensor_rank_2_ij, array_ij))
        });
}

#[test]
fn lu_decomposition() {
    let (tensor_l, tensor_u) = get_tensor_rank_2_dim_9().lu_decomposition();
    tensor_l
        .iter()
        .enumerate()
        .zip(tensor_u.iter())
        .for_each(|((i, tensor_l_i), tensor_u_i)| {
            tensor_l_i
                .iter()
                .enumerate()
                .zip(tensor_u_i.iter())
                .for_each(|((j, tensor_l_ij), tensor_u_ij)| match i.cmp(&j) {
                    Ordering::Greater => assert_eq!(tensor_u_ij, &0.0),
                    Ordering::Less => assert_eq!(tensor_l_ij, &0.0),
                    _ => (),
                })
        });
}

#[test]
fn mul_tensor_rank_0_to_self() {
    (get_tensor_rank_2_dim_4() * 3.3)
        .iter()
        .zip(get_array_dim_4().iter())
        .for_each(|(tensor_rank_2_i, array_i)| {
            tensor_rank_2_i
                .iter()
                .zip(array_i.iter())
                .for_each(|(tensor_rank_2_ij, array_ij)| {
                    assert_eq!(tensor_rank_2_ij, &(array_ij * 3.3))
                })
        });
}

#[test]
fn mul_tensor_rank_0_to_self_ref() {
    (&get_tensor_rank_2_dim_4() * 3.3)
        .iter()
        .zip(get_array_dim_4().iter())
        .for_each(|(tensor_rank_2_i, array_i)| {
            tensor_rank_2_i
                .iter()
                .zip(array_i.iter())
                .for_each(|(tensor_rank_2_ij, array_ij)| {
                    assert_eq!(tensor_rank_2_ij, &(array_ij * 3.3))
                })
        });
}

#[test]
#[allow(clippy::op_ref)]
fn mul_tensor_rank_0_ref_to_self() {
    (get_tensor_rank_2_dim_4() * &3.3)
        .iter()
        .zip(get_array_dim_4().iter())
        .for_each(|(tensor_rank_2_i, array_i)| {
            tensor_rank_2_i
                .iter()
                .zip(array_i.iter())
                .for_each(|(tensor_rank_2_ij, array_ij)| {
                    assert_eq!(tensor_rank_2_ij, &(array_ij * 3.3))
                })
        });
}

#[test]
#[allow(clippy::op_ref)]
fn mul_tensor_rank_0_ref_to_self_ref() {
    (&get_tensor_rank_2_dim_4() * &3.3)
        .iter()
        .zip(get_array_dim_4().iter())
        .for_each(|(tensor_rank_2_i, array_i)| {
            tensor_rank_2_i
                .iter()
                .zip(array_i.iter())
                .for_each(|(tensor_rank_2_ij, array_ij)| {
                    assert_eq!(tensor_rank_2_ij, &(array_ij * 3.3))
                })
        });
}

#[test]
fn mul_assign_tensor_rank_0() {
    let mut tensor_rank_2 = get_tensor_rank_2_dim_4();
    tensor_rank_2 *= 3.3;
    tensor_rank_2
        .iter()
        .zip(get_array_dim_4().iter())
        .for_each(|(tensor_rank_2_i, array_i)| {
            tensor_rank_2_i
                .iter()
                .zip(array_i.iter())
                .for_each(|(tensor_rank_2_ij, array_ij)| {
                    assert_eq!(tensor_rank_2_ij, &(array_ij * 3.3))
                })
        });
}

#[test]
fn mul_assign_tensor_rank_0_ref() {
    let mut tensor_rank_2 = get_tensor_rank_2_dim_4();
    tensor_rank_2 *= &3.3;
    tensor_rank_2
        .iter()
        .zip(get_array_dim_4().iter())
        .for_each(|(tensor_rank_2_i, array_i)| {
            tensor_rank_2_i
                .iter()
                .zip(array_i.iter())
                .for_each(|(tensor_rank_2_ij, array_ij)| {
                    assert_eq!(tensor_rank_2_ij, &(array_ij * 3.3))
                })
        });
}

#[test]
fn mul_tensor_rank_1_to_self() {
    (get_tensor_rank_2_dim_4() * get_tensor_rank_1_a())
        .iter()
        .zip(get_other_tensor_rank_2_mul_tensor_rank_1_dim_4().iter())
        .for_each(|(tensor_rank_1_i, res_tensor_rank_1_i)| {
            assert_eq!(tensor_rank_1_i, res_tensor_rank_1_i)
        });
}

#[test]
fn mul_tensor_rank_1_ref_to_self() {
    (get_tensor_rank_2_dim_4() * &get_tensor_rank_1_a())
        .iter()
        .zip(get_other_tensor_rank_2_mul_tensor_rank_1_dim_4().iter())
        .for_each(|(tensor_rank_1_i, res_tensor_rank_1_i)| {
            assert_eq!(tensor_rank_1_i, res_tensor_rank_1_i)
        });
}

#[test]
fn mul_tensor_rank_1_to_self_ref() {
    (&get_tensor_rank_2_dim_4() * get_tensor_rank_1_a())
        .iter()
        .zip(get_other_tensor_rank_2_mul_tensor_rank_1_dim_4().iter())
        .for_each(|(tensor_rank_1_i, res_tensor_rank_1_i)| {
            assert_eq!(tensor_rank_1_i, res_tensor_rank_1_i)
        });
}

#[test]
fn mul_tensor_rank_1_ref_to_self_ref() {
    (&get_tensor_rank_2_dim_4() * &get_tensor_rank_1_a())
        .iter()
        .zip(get_other_tensor_rank_2_mul_tensor_rank_1_dim_4().iter())
        .for_each(|(tensor_rank_1_i, res_tensor_rank_1_i)| {
            assert_eq!(tensor_rank_1_i, res_tensor_rank_1_i)
        });
}

#[test]
fn mul_tensor_rank_2_to_self() {
    (get_tensor_rank_2_dim_4() * get_other_tensor_rank_2_dim_4())
        .iter()
        .zip(get_other_tensor_rank_2_mul_tensor_rank_2_dim_4().iter())
        .for_each(|(tensor_rank_2_i, res_tensor_rank_2_i)| {
            tensor_rank_2_i
                .iter()
                .zip(res_tensor_rank_2_i.iter())
                .for_each(|(tensor_rank_2_ij, res_tensor_rank_2_ij)| {
                    assert_eq!(tensor_rank_2_ij, res_tensor_rank_2_ij)
                })
        });
}

#[test]
fn mul_tensor_rank_2_ref_to_self() {
    (get_tensor_rank_2_dim_4() * &get_other_tensor_rank_2_dim_4())
        .iter()
        .zip(get_other_tensor_rank_2_mul_tensor_rank_2_dim_4().iter())
        .for_each(|(tensor_rank_2_i, res_tensor_rank_2_i)| {
            tensor_rank_2_i
                .iter()
                .zip(res_tensor_rank_2_i.iter())
                .for_each(|(tensor_rank_2_ij, res_tensor_rank_2_ij)| {
                    assert_eq!(tensor_rank_2_ij, res_tensor_rank_2_ij)
                })
        });
}

#[test]
fn mul_tensor_rank_2_to_self_ref() {
    (&get_tensor_rank_2_dim_4() * get_other_tensor_rank_2_dim_4())
        .iter()
        .zip(get_other_tensor_rank_2_mul_tensor_rank_2_dim_4().iter())
        .for_each(|(tensor_rank_2_i, res_tensor_rank_2_i)| {
            tensor_rank_2_i
                .iter()
                .zip(res_tensor_rank_2_i.iter())
                .for_each(|(tensor_rank_2_ij, res_tensor_rank_2_ij)| {
                    assert_eq!(tensor_rank_2_ij, res_tensor_rank_2_ij)
                })
        });
}

#[test]
fn mul_tensor_rank_2_ref_to_self_ref() {
    (&get_tensor_rank_2_dim_4() * &get_other_tensor_rank_2_dim_4())
        .iter()
        .zip(get_other_tensor_rank_2_mul_tensor_rank_2_dim_4().iter())
        .for_each(|(tensor_rank_2_i, res_tensor_rank_2_i)| {
            tensor_rank_2_i
                .iter()
                .zip(res_tensor_rank_2_i.iter())
                .for_each(|(tensor_rank_2_ij, res_tensor_rank_2_ij)| {
                    assert_eq!(tensor_rank_2_ij, res_tensor_rank_2_ij)
                })
        });
}

#[test]
fn mul_tensor_rank_1_list_to_self() {
    (get_tensor_rank_2_dim_3() * get_tensor_rank_1_list())
        .iter()
        .zip(get_tensor_rank_1_list().iter())
        .for_each(|(res_tensor_rank_1, tensor_rank_1)| {
            res_tensor_rank_1
                .iter()
                .zip((get_tensor_rank_2_dim_3() * tensor_rank_1).iter())
                .for_each(|(res_tensor_rank_1_i, tensor_rank_1_i)| {
                    assert_eq!(res_tensor_rank_1_i, tensor_rank_1_i)
                })
        })
}

#[test]
fn mul_tensor_rank_1_list_ref_to_self() {
    (get_tensor_rank_2_dim_3() * &get_tensor_rank_1_list())
        .iter()
        .zip(get_tensor_rank_1_list().iter())
        .for_each(|(res_tensor_rank_1, tensor_rank_1)| {
            res_tensor_rank_1
                .iter()
                .zip((get_tensor_rank_2_dim_3() * tensor_rank_1).iter())
                .for_each(|(res_tensor_rank_1_i, tensor_rank_1_i)| {
                    assert_eq!(res_tensor_rank_1_i, tensor_rank_1_i)
                })
        })
}

#[test]
fn mul_tensor_rank_1_list_to_self_ref() {
    (&get_tensor_rank_2_dim_3() * get_tensor_rank_1_list())
        .iter()
        .zip(get_tensor_rank_1_list().iter())
        .for_each(|(res_tensor_rank_1, tensor_rank_1)| {
            res_tensor_rank_1
                .iter()
                .zip((get_tensor_rank_2_dim_3() * tensor_rank_1).iter())
                .for_each(|(res_tensor_rank_1_i, tensor_rank_1_i)| {
                    assert_eq!(res_tensor_rank_1_i, tensor_rank_1_i)
                })
        })
}

#[test]
fn mul_tensor_rank_1_list_ref_to_self_ref() {
    (&get_tensor_rank_2_dim_3() * &get_tensor_rank_1_list())
        .iter()
        .zip(get_tensor_rank_1_list().iter())
        .for_each(|(res_tensor_rank_1, tensor_rank_1)| {
            res_tensor_rank_1
                .iter()
                .zip((get_tensor_rank_2_dim_3() * tensor_rank_1).iter())
                .for_each(|(res_tensor_rank_1_i, tensor_rank_1_i)| {
                    assert_eq!(res_tensor_rank_1_i, tensor_rank_1_i)
                })
        })
}

#[test]
fn mul_tensor_rank_2_list_2d_to_self() {
    (get_tensor_rank_2_dim_3() * get_tensor_rank_2_list_2d())
        .iter()
        .zip(get_tensor_rank_2_mul_tensor_rank_2_list_2d().iter())
        .for_each(|(tensor_rank_2_list_2d_entry, res_entry)| {
            tensor_rank_2_list_2d_entry
                .iter()
                .zip(res_entry.iter())
                .for_each(|(tensor_rank_2, res)| {
                    tensor_rank_2
                        .iter()
                        .zip(res.iter())
                        .for_each(|(tensor_rank_2_i, res_i)| {
                            tensor_rank_2_i.iter().zip(res_i.iter()).for_each(
                                |(tensor_rank_2_ij, res_ij)| assert_eq!(tensor_rank_2_ij, res_ij),
                            )
                        })
                })
        });
}

#[test]
fn mul_tensor_rank_2_list_2d_to_self_ref() {
    (&get_tensor_rank_2_dim_3() * get_tensor_rank_2_list_2d())
        .iter()
        .zip(get_tensor_rank_2_mul_tensor_rank_2_list_2d().iter())
        .for_each(|(tensor_rank_2_list_2d_entry, res_entry)| {
            tensor_rank_2_list_2d_entry
                .iter()
                .zip(res_entry.iter())
                .for_each(|(tensor_rank_2, res)| {
                    tensor_rank_2
                        .iter()
                        .zip(res.iter())
                        .for_each(|(tensor_rank_2_i, res_i)| {
                            tensor_rank_2_i.iter().zip(res_i.iter()).for_each(
                                |(tensor_rank_2_ij, res_ij)| assert_eq!(tensor_rank_2_ij, res_ij),
                            )
                        })
                })
        });
}

#[test]
fn new() {
    get_tensor_rank_2_dim_4()
        .iter()
        .zip(get_array_dim_4().iter())
        .for_each(|(tensor_rank_2_i, array_i)| {
            tensor_rank_2_i
                .iter()
                .zip(array_i.iter())
                .for_each(|(tensor_rank_2_ij, array_ij)| assert_eq!(tensor_rank_2_ij, array_ij))
        });
}

#[test]
fn norm_dim_2() -> Result<(), TestError> {
    assert_eq(&get_tensor_rank_2_dim_2().norm(), &5.477_225_575_051_661)
}

#[test]
fn norm_dim_3() -> Result<(), TestError> {
    assert_eq(&get_tensor_rank_2_dim_3().norm(), &16.881_943_016_134_134)
}

#[test]
fn norm_dim_4() -> Result<(), TestError> {
    assert_eq(&get_tensor_rank_2_dim_4().norm(), &14.282_856_857_085_7)
}

#[test]
fn norm_dim_9() -> Result<(), TestError> {
    assert_eq(&get_tensor_rank_2_dim_9().norm(), &20.976_176_963_403_03)
}

#[test]
fn size() {
    assert_eq!(
        std::mem::size_of::<TensorRank2::<3, 1, 1>>(),
        std::mem::size_of::<[TensorRank1::<3, 1>; 3]>()
    )
}

#[test]
fn second_invariant() {
    assert_eq!(get_tensor_rank_2_dim_4().second_invariant(), 16.0);
}

#[test]
fn squared_trace_dim_2() -> Result<(), TestError> {
    assert_eq_within_tols(&get_tensor_rank_2_dim_2().squared_trace(), &29.0)
}

#[test]
fn squared_trace_dim_3() -> Result<(), TestError> {
    assert_eq_within_tols(&get_tensor_rank_2_dim_3().squared_trace(), &258.0)
}

#[test]
fn squared_trace_dim_4() -> Result<(), TestError> {
    assert_eq_within_tols(&get_tensor_rank_2_dim_4().squared_trace(), &89.0)
}

#[test]
fn squared_trace_dim_9() -> Result<(), TestError> {
    assert_eq_within_tols(&get_tensor_rank_2_dim_9().squared_trace(), &318.0)
}

#[test]
fn sub_tensor_rank_2_to_self() {
    (get_tensor_rank_2_dim_4() - get_other_tensor_rank_2_dim_4())
        .iter()
        .zip(get_other_tensor_rank_2_sub_tensor_rank_2_dim_4().iter())
        .for_each(|(tensor_rank_2_i, res_tensor_rank_2_i)| {
            tensor_rank_2_i
                .iter()
                .zip(res_tensor_rank_2_i.iter())
                .for_each(|(tensor_rank_2_ij, res_tensor_rank_2_ij)| {
                    assert_eq!(tensor_rank_2_ij, res_tensor_rank_2_ij)
                })
        });
}

#[test]
fn sub_tensor_rank_2_ref_to_self() {
    (get_tensor_rank_2_dim_4() - &get_other_tensor_rank_2_dim_4())
        .iter()
        .zip(get_other_tensor_rank_2_sub_tensor_rank_2_dim_4().iter())
        .for_each(|(tensor_rank_2_i, res_tensor_rank_2_i)| {
            tensor_rank_2_i
                .iter()
                .zip(res_tensor_rank_2_i.iter())
                .for_each(|(tensor_rank_2_ij, res_tensor_rank_2_ij)| {
                    assert_eq!(tensor_rank_2_ij, res_tensor_rank_2_ij)
                })
        });
}

#[test]
fn sub_assign_tensor_rank_2() {
    let mut tensor_rank_2 = get_tensor_rank_2_dim_4();
    tensor_rank_2 -= get_other_tensor_rank_2_dim_4();
    tensor_rank_2
        .iter()
        .zip(get_other_tensor_rank_2_sub_tensor_rank_2_dim_4().iter())
        .for_each(|(tensor_rank_2_i, res_tensor_rank_2_i)| {
            tensor_rank_2_i
                .iter()
                .zip(res_tensor_rank_2_i.iter())
                .for_each(|(tensor_rank_2_ij, res_tensor_rank_2_ij)| {
                    assert_eq!(tensor_rank_2_ij, res_tensor_rank_2_ij)
                })
        });
}

#[test]
fn sub_assign_tensor_rank_2_ref() {
    let mut tensor_rank_2 = get_tensor_rank_2_dim_4();
    tensor_rank_2 -= &get_other_tensor_rank_2_dim_4();
    tensor_rank_2
        .iter()
        .zip(get_other_tensor_rank_2_sub_tensor_rank_2_dim_4().iter())
        .for_each(|(tensor_rank_2_i, res_tensor_rank_2_i)| {
            tensor_rank_2_i
                .iter()
                .zip(res_tensor_rank_2_i.iter())
                .for_each(|(tensor_rank_2_ij, res_tensor_rank_2_ij)| {
                    assert_eq!(tensor_rank_2_ij, res_tensor_rank_2_ij)
                })
        });
}

#[test]
fn trace_dim_2() -> Result<(), TestError> {
    assert_eq(&get_tensor_rank_2_dim_2().trace(), &5.0)
}

#[test]
fn trace_dim_3() -> Result<(), TestError> {
    assert_eq(&get_tensor_rank_2_dim_3().trace(), &6.0)
}

#[test]
fn trace_dim_4() -> Result<(), TestError> {
    assert_eq(&get_tensor_rank_2_dim_4().trace(), &11.0)
}

#[test]
fn trace_dim_9() -> Result<(), TestError> {
    assert_eq(&get_tensor_rank_2_dim_9().trace(), &18.0)
}

#[test]
fn transpose() {
    let tensor_rank_2 = get_tensor_rank_2_dim_4();
    let tensor_rank_2_transpose = tensor_rank_2.transpose();
    tensor_rank_2
        .iter()
        .enumerate()
        .for_each(|(i, tensor_rank_2_i)| {
            tensor_rank_2_i
                .iter()
                .enumerate()
                .for_each(|(j, tensor_rank_2_ij)| {
                    assert_eq!(tensor_rank_2_ij, &tensor_rank_2_transpose[j][i])
                })
        });
    tensor_rank_2_transpose
        .iter()
        .enumerate()
        .for_each(|(i, tensor_rank_2_transpose_i)| {
            tensor_rank_2_transpose_i.iter().enumerate().for_each(
                |(j, tensor_rank_2_transpose_ij)| {
                    assert_eq!(tensor_rank_2_transpose_ij, &tensor_rank_2[j][i])
                },
            )
        });
}

#[test]
fn zero_dim_2() {
    TensorRank2::<2, 1, 1>::zero()
        .iter()
        .for_each(|tensor_rank_2_i| {
            tensor_rank_2_i
                .iter()
                .for_each(|tensor_rank_2_ij| assert_eq!(tensor_rank_2_ij, &0.0))
        });
}

#[test]
fn zero_dim_3() {
    TensorRank2::<3, 1, 1>::zero()
        .iter()
        .for_each(|tensor_rank_2_i| {
            tensor_rank_2_i
                .iter()
                .for_each(|tensor_rank_2_ij| assert_eq!(tensor_rank_2_ij, &0.0))
        });
}

#[test]
fn zero_dim_4() {
    TensorRank2::<4, 1, 1>::zero()
        .iter()
        .for_each(|tensor_rank_2_i| {
            tensor_rank_2_i
                .iter()
                .for_each(|tensor_rank_2_ij| assert_eq!(tensor_rank_2_ij, &0.0))
        });
}

#[test]
fn zero_dim_9() {
    TensorRank2::<9, 1, 1>::zero()
        .iter()
        .for_each(|tensor_rank_2_i| {
            tensor_rank_2_i
                .iter()
                .for_each(|tensor_rank_2_ij| assert_eq!(tensor_rank_2_ij, &0.0))
        });
}
