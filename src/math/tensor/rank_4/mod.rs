#[cfg(test)]
mod test;

use std::ops::
{
    Add,
    AddAssign,
    Div,
    DivAssign,
    Index,
    IndexMut,
    Mul,
    MulAssign,
    Sub,
    SubAssign
};

use super::
{
    rank_0::TensorRank0,
    rank_1::
    {
        TensorRank1,
        TensorRank1Trait
    },
    rank_2::
    {
        TensorRank2,
        TensorRank2Trait
    },
    rank_3::TensorRank3
};

pub struct TensorRank4<const D: usize, const I: usize, const J: usize, const K: usize, const L: usize>
(
    [TensorRank3<D, J, K, L>; D]
);

impl<const D: usize, const I: usize, const J: usize, const K: usize, const L: usize> TensorRank4<D, I, J, K, L>
{
    pub fn iter(&self) -> impl Iterator<Item=&TensorRank3<D, J, K, L>>
    {
        self.0.iter()
    }
    pub fn iter_mut(&mut self) -> impl Iterator<Item=&mut TensorRank3<D, J, K, L>>
    {
        self.0.iter_mut()
    }
    pub fn zero() -> Self
    {
        Self(std::array::from_fn(|_| TensorRank3::zero()))
    }
}

pub trait TensorRank4Trait<const D: usize, T2A, T2B, T2C, T2D>
{
    type OutputAsTensorRank2;
    type OutputContractAllIndicesWithFirstIndicesOf;
    type OutputContractFirstThirdFourthIndicesWithFirstIndicesOf;
    type OutputContractSecondIndexWithFirstIndexOf;
    fn as_tensor_rank_2(&self) -> Self::OutputAsTensorRank2;
    fn contract_all_indices_with_first_indices_of(&self, tensor_rank_2_a: T2A, tensor_rank_2_b: T2B, tensor_rank_2_c: T2C, tensor_rank_2_d: T2D) -> Self::OutputContractAllIndicesWithFirstIndicesOf;
    fn contract_first_third_fourth_indices_with_first_indices_of(&self, tensor_rank_2_a: T2A, tensor_rank_2_c: T2C, tensor_rank_2_d: T2D) -> Self::OutputContractFirstThirdFourthIndicesWithFirstIndicesOf;
    fn contract_second_index_with_first_index_of(&self, tensor_rank_2_b: T2B) -> Self::OutputContractSecondIndexWithFirstIndexOf;
    fn dyad_ij_kl(tensor_rank_2_a: T2A, tensor_rank_2_b: T2B) -> Self;
    fn dyad_ik_jl(tensor_rank_2_a: T2A, tensor_rank_2_b: T2B) -> Self;
    fn dyad_il_jk(tensor_rank_2_a: T2A, tensor_rank_2_b: T2B) -> Self;
    fn dyad_il_kj(tensor_rank_2_a: T2A, tensor_rank_2_b: T2B) -> Self;
    fn new(array: [[[[TensorRank0; D]; D]; D]; D]) -> Self;
    fn inverse(&self) -> Self;
}

impl<const I: usize, const J: usize, const K: usize, const L: usize, const M: usize, const N: usize, const O: usize, const P: usize> TensorRank4Trait<3, &TensorRank2<3, I, M>, &TensorRank2<3, J, N>, &TensorRank2<3, K, O>, &TensorRank2<3, L, P>> for TensorRank4<3, I, J, K, L>
{
    // I don't see how to make this work and retain the configurations here without passing them through somehow
    type OutputAsTensorRank2 = TensorRank2<9, I, J>;
    type OutputContractAllIndicesWithFirstIndicesOf = TensorRank4<3, M, N, O, P>;
    type OutputContractFirstThirdFourthIndicesWithFirstIndicesOf = TensorRank4<3, M, J, O, P>;
    type OutputContractSecondIndexWithFirstIndexOf = TensorRank4<3, I, N, K, L>;
    fn as_tensor_rank_2(&self) -> Self::OutputAsTensorRank2
    {
        // let mut tensor_rank_2 = TensorRank2::zero();
        // self.iter().enumerate().for_each(|(i, self_i)|
        //     self_i.iter().enumerate().for_each(|(j, self_ij)|
        //         self_ij.iter().enumerate().for_each(|(k, self_ijk)|
        //             self_ijk.iter().enumerate().for_each(|(l, self_ijkl)|
        //                 tensor_rank_2[3*i + j][3*k + l] = *self_ijkl
        //             )
        //         )
        //     )
        // );
        // tensor_rank_2
        todo!()
    }
    fn contract_all_indices_with_first_indices_of(&self, tensor_rank_2_a: &TensorRank2<3, I, M>, tensor_rank_2_b: &TensorRank2<3, J, N>, tensor_rank_2_c: &TensorRank2<3, K, O>, tensor_rank_2_d: &TensorRank2<3, L, P>) -> Self::OutputContractAllIndicesWithFirstIndicesOf
    {
        let mut output = TensorRank4::zero();
        self.iter().zip(tensor_rank_2_a.iter()).for_each(|(self_m, tensor_rank_2_a_m)|
            self_m.iter().zip(tensor_rank_2_b.iter()).for_each(|(self_mn, tensor_rank_2_b_n)|
                self_mn.iter().zip(tensor_rank_2_c.iter()).for_each(|(self_mno, tensor_rank_2_c_o)|
                    self_mno.iter().zip(tensor_rank_2_d.iter()).for_each(|(self_mnop, tensor_rank_2_d_p)|
                        output.iter_mut().zip(tensor_rank_2_a_m.iter()).for_each(|(output_i, tensor_rank_2_a_mi)|
                            output_i.iter_mut().zip(tensor_rank_2_b_n.iter()).for_each(|(output_ij, tensor_rank_2_b_nj)|
                                output_ij.iter_mut().zip(tensor_rank_2_c_o.iter()).for_each(|(output_ijk, tensor_rank_2_c_ok)|
                                    output_ijk.iter_mut().zip(tensor_rank_2_d_p.iter()).for_each(|(output_ijkl, tensor_rank_2_dp)|
                                        *output_ijkl += self_mnop*tensor_rank_2_a_mi*tensor_rank_2_b_nj*tensor_rank_2_c_ok*tensor_rank_2_dp
                                    )
                                )
                            )
                        )
                    )
                )
            )
        );
        output
    }
    fn contract_first_third_fourth_indices_with_first_indices_of(&self, tensor_rank_2_a: &TensorRank2<3, I, M>, tensor_rank_2_b: &TensorRank2<3, K, O>, tensor_rank_2_c: &TensorRank2<3, L, P>) -> Self::OutputContractFirstThirdFourthIndicesWithFirstIndicesOf
    {
        let mut output = TensorRank4::zero();
        self.iter().zip(tensor_rank_2_a.iter()).for_each(|(self_q, tensor_rank_2_a_q)|
            output.iter_mut().zip(tensor_rank_2_a_q.iter()).for_each(|(output_i, tensor_rank_2_a_qi)|    
                output_i.iter_mut().zip(self_q.iter()).for_each(|(output_ij, self_qj)|
                    self_qj.iter().zip(tensor_rank_2_b.iter()).for_each(|(self_qjm, tensor_rank_2_b_m)|
                        self_qjm.iter().zip(tensor_rank_2_c.iter()).for_each(|(self_qjmn, tensor_rank_2_c_n)|
                            output_ij.iter_mut().zip(tensor_rank_2_b_m.iter()).for_each(|(output_ijk, tensor_rank_2_b_mk)|
                                output_ijk.iter_mut().zip(tensor_rank_2_c_n.iter()).for_each(|(output_ijkl, tensor_rank_2_c_nl)|
                                    *output_ijkl += self_qjmn*tensor_rank_2_a_qi*tensor_rank_2_b_mk*tensor_rank_2_c_nl
                                )
                            )
                        )
                    )
                )
            )
        );
        output
    }
    fn contract_second_index_with_first_index_of(&self, tensor_rank_2: &TensorRank2<3, J, N>) -> Self::OutputContractSecondIndexWithFirstIndexOf
    {
        let mut output = TensorRank4::zero();
        output.iter_mut().zip(self.iter()).for_each(|(output_i, self_i)|
            self_i.iter().zip(tensor_rank_2.iter()).for_each(|(self_is, tensor_rank_2_s)|
                output_i.iter_mut().zip(tensor_rank_2_s.iter()).for_each(|(output_ij, tensor_rank_2_sj)|
                    *output_ij += self_is * tensor_rank_2_sj
                )
            )
        );
        output
    }
    fn dyad_ij_kl(tensor_rank_2_a: &TensorRank2<3, I, J>, tensor_rank_2_b: &TensorRank2<3, K, L>) -> Self
    {
        tensor_rank_2_a.iter().map(|tensor_rank_2_a_i|
            tensor_rank_2_a_i.iter().map(|tensor_rank_2_a_ij|
                tensor_rank_2_b * tensor_rank_2_a_ij
            ).collect()
        ).collect()
    }
    fn dyad_ik_jl(tensor_rank_2_a: &TensorRank2<3, I, K>, tensor_rank_2_b: &TensorRank2<3, J, L>) -> Self
    {
        tensor_rank_2_a.iter().map(|tensor_rank_2_a_i|
            tensor_rank_2_b.iter().map(|tensor_rank_2_b_j|
                tensor_rank_2_a_i.iter().map(|tensor_rank_2_a_ik|
                    tensor_rank_2_b_j * tensor_rank_2_a_ik
                ).collect()
            ).collect()
        ).collect()
    }
    fn dyad_il_jk(tensor_rank_2_a: &TensorRank2<3, I, L>, tensor_rank_2_b: &TensorRank2<3, J, K>) -> Self
    {
        tensor_rank_2_a.iter().map(|tensor_rank_2_a_i|
            tensor_rank_2_b.iter().map(|tensor_rank_2_b_j|
                tensor_rank_2_b_j.iter().map(|tensor_rank_2_b_jk|
                    tensor_rank_2_a_i * tensor_rank_2_b_jk
                ).collect()
            ).collect()
        ).collect()
    }
    fn dyad_il_kj(tensor_rank_2_a: &TensorRank2<3, I, L>, tensor_rank_2_b: &TensorRank2<3, K, J>) -> Self
    {
        Self::dyad_il_jk(tensor_rank_2_a, &(tensor_rank_2_b.transpose()))
    }
    fn new(array: [[[[TensorRank0; 3]; 3]; 3]; 3]) -> Self
    {
        array.iter().map(|array_i|
            array_i.iter().map(|array_ij|
                array_ij.iter().map(|array_ijk|
                    TensorRank1::new(*array_ijk)
                ).collect()
            ).collect()
        ).collect()
    }
    fn inverse(&self) -> Self
    {
        // self.as_tensor_rank_2().inverse().as_tensor_rank_4()
        todo!()
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize, const L: usize> FromIterator<TensorRank3<D, J, K, L>> for TensorRank4<D, I, J, K, L>
{
    fn from_iter<Ii: IntoIterator<Item=TensorRank3<D, J, K, L>>>(into_iterator: Ii) -> Self
    {
        let mut tensor_rank_4 = Self::zero();
        tensor_rank_4.iter_mut().zip(into_iterator).for_each(|(tensor_rank_4_i, value_i)|
            *tensor_rank_4_i = value_i
        );
        tensor_rank_4
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize, const L: usize> Index<usize> for TensorRank4<D, I, J, K, L>
{
    type Output = TensorRank3<D, J, K, L>;
    fn index(&self, index: usize) -> &Self::Output
    {
        &self.0[index]
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize, const L: usize> IndexMut<usize> for TensorRank4<D, I, J, K, L>
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output
    {
        &mut self.0[index]
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize, const L: usize> Div<TensorRank0> for TensorRank4<D, I, J, K, L>
{
    type Output = Self;
    fn div(mut self, tensor_rank_0: TensorRank0) -> Self::Output
    {
        self /= &tensor_rank_0;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize, const L: usize> Div<&TensorRank0> for TensorRank4<D, I, J, K, L>
{
    type Output = Self;
    fn div(mut self, tensor_rank_0: &TensorRank0) -> Self::Output
    {
        self /= tensor_rank_0;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize, const L: usize> DivAssign<TensorRank0> for TensorRank4<D, I, J, K, L>
{
    fn div_assign(&mut self, tensor_rank_0: TensorRank0)
    {
        self.iter_mut().for_each(|self_i|
            *self_i /= &tensor_rank_0
        );
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize, const L: usize> DivAssign<&TensorRank0> for TensorRank4<D, I, J, K, L>
{
    fn div_assign(&mut self, tensor_rank_0: &TensorRank0)
    {
        self.iter_mut().for_each(|self_i|
            *self_i /= tensor_rank_0
        );
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize, const L: usize> Mul<TensorRank0> for TensorRank4<D, I, J, K, L>
{
    type Output = Self;
    fn mul(mut self, tensor_rank_0: TensorRank0) -> Self::Output
    {
        self *= &tensor_rank_0;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize, const L: usize> Mul<&TensorRank0> for TensorRank4<D, I, J, K, L>
{
    type Output = Self;
    fn mul(mut self, tensor_rank_0: &TensorRank0) -> Self::Output
    {
        self *= tensor_rank_0;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize, const L: usize> MulAssign<TensorRank0> for TensorRank4<D, I, J, K, L>
{
    fn mul_assign(&mut self, tensor_rank_0: TensorRank0)
    {
        self.iter_mut().for_each(|self_i|
            *self_i *= &tensor_rank_0
        );
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize, const L: usize> MulAssign<&TensorRank0> for TensorRank4<D, I, J, K, L>
{
    fn mul_assign(&mut self, tensor_rank_0: &TensorRank0)
    {
        self.iter_mut().for_each(|self_i|
            *self_i *= tensor_rank_0
        );
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize, const L: usize, const M: usize> Mul<TensorRank2<D, L, M>> for TensorRank4<D, I, J, K, L>
// where
//     TensorRank2<D>: TensorRank2Trait<D>
{
    type Output = TensorRank4<D, I, J, K, M>;
    fn mul(self, tensor_rank_2: TensorRank2<D, L, M>) -> Self::Output
    {
        self.iter().map(|self_i|
            self_i.iter().map(|self_ij|
                self_ij * &tensor_rank_2
            ).collect()
        ).collect()
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize, const L: usize, const M: usize> Mul<&TensorRank2<D, L, M>> for TensorRank4<D, I, J, K, L>
// where
//     TensorRank2<D>: TensorRank2Trait<D>
{
    type Output = TensorRank4<D, I, J, K, M>;
    fn mul(self, tensor_rank_2: &TensorRank2<D, L, M>) -> Self::Output
    {
        self.iter().map(|self_i|
            self_i.iter().map(|self_ij|
                self_ij * tensor_rank_2
            ).collect()
        ).collect()
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize, const L: usize> Add for TensorRank4<D, I, J, K, L>
{
    type Output = Self;
    fn add(mut self, tensor_rank_4: Self) -> Self::Output
    {
        self += tensor_rank_4;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize, const L: usize> Add<&Self> for TensorRank4<D, I, J, K, L>
{
    type Output = Self;
    fn add(mut self, tensor_rank_4: &Self) -> Self::Output
    {
        self += tensor_rank_4;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize, const L: usize> Add<TensorRank4<D, I, J, K, L>> for &TensorRank4<D, I, J, K, L>
{
    type Output = TensorRank4<D, I, J, K, L>;
    fn add(self, mut tensor_rank_4: TensorRank4<D, I, J, K, L>) -> Self::Output
    {
        tensor_rank_4 += self;
        tensor_rank_4
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize, const L: usize> AddAssign for TensorRank4<D, I, J, K, L>
{
    fn add_assign(&mut self, tensor_rank_4: Self)
    {
        self.iter_mut().zip(tensor_rank_4.iter()).for_each(|(self_i, tensor_rank_4_i)|
            *self_i += tensor_rank_4_i
        );
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize, const L: usize> AddAssign<&Self> for TensorRank4<D, I, J, K, L>
{
    fn add_assign(&mut self, tensor_rank_4: &Self)
    {
        self.iter_mut().zip(tensor_rank_4.iter()).for_each(|(self_i, tensor_rank_4_i)|
            *self_i += tensor_rank_4_i
        );
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize, const L: usize> Sub for TensorRank4<D, I, J, K, L>
{
    type Output = Self;
    fn sub(mut self, tensor_rank_4: Self) -> Self::Output
    {
        self -= tensor_rank_4;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize, const L: usize> Sub<&Self> for TensorRank4<D, I, J, K, L>
{
    type Output = Self;
    fn sub(mut self, tensor_rank_4: &Self) -> Self::Output
    {
        self -= tensor_rank_4;
        self
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize, const L: usize> SubAssign for TensorRank4<D, I, J, K, L>
{
    fn sub_assign(&mut self, tensor_rank_4: Self)
    {
        self.iter_mut().zip(tensor_rank_4.iter()).for_each(|(self_i, tensor_rank_4_i)|
            *self_i -= tensor_rank_4_i
        );
    }
}

impl<const D: usize, const I: usize, const J: usize, const K: usize, const L: usize> SubAssign<&Self> for TensorRank4<D, I, J, K, L>
{
    fn sub_assign(&mut self, tensor_rank_4: &Self)
    {
        self.iter_mut().zip(tensor_rank_4.iter()).for_each(|(self_i, tensor_rank_4_i)|
            *self_i -= tensor_rank_4_i
        );
    }
}