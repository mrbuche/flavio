use super::{TensorRank0, TensorRank4, TensorRank4List, TensorRank4ListTrait};

fn get_array() -> [[[[[TensorRank0; 3]; 3]; 3]; 3]; 2] {
    [
        [
            [
                [[4.0, 2.0, 4.0], [1.0, 4.0, 3.0], [2.0, 4.0, 4.0]],
                [[2.0, 2.0, 2.0], [3.0, 1.0, 1.0], [1.0, 4.0, 2.0]],
                [[1.0, 2.0, 3.0], [2.0, 2.0, 3.0], [1.0, 1.0, 0.0]],
            ],
            [
                [[2.0, 4.0, 2.0], [1.0, 2.0, 3.0], [3.0, 3.0, 2.0]],
                [[1.0, 1.0, 1.0], [4.0, 2.0, 1.0], [1.0, 4.0, 1.0]],
                [[2.0, 2.0, 4.0], [3.0, 3.0, 1.0], [0.0, 3.0, 3.0]],
            ],
            [
                [[0.0, 1.0, 4.0], [3.0, 3.0, 3.0], [4.0, 4.0, 0.0]],
                [[2.0, 3.0, 1.0], [1.0, 2.0, 0.0], [2.0, 2.0, 4.0]],
                [[3.0, 4.0, 1.0], [2.0, 1.0, 2.0], [4.0, 4.0, 1.0]],
            ],
        ],
        [
            [
                [[2.0, 2.0, 4.0], [0.0, 0.0, 1.0], [1.0, 3.0, 3.0]],
                [[0.0, 3.0, 1.0], [0.0, 0.0, 1.0], [4.0, 2.0, 1.0]],
                [[3.0, 0.0, 1.0], [2.0, 0.0, 3.0], [4.0, 4.0, 2.0]],
            ],
            [
                [[4.0, 4.0, 0.0], [2.0, 1.0, 1.0], [0.0, 0.0, 4.0]],
                [[0.0, 1.0, 0.0], [1.0, 1.0, 3.0], [0.0, 1.0, 1.0]],
                [[4.0, 2.0, 3.0], [4.0, 2.0, 4.0], [3.0, 0.0, 4.0]],
            ],
            [
                [[1.0, 3.0, 2.0], [0.0, 0.0, 0.0], [2.0, 4.0, 2.0]],
                [[2.0, 2.0, 2.0], [4.0, 1.0, 2.0], [4.0, 2.0, 2.0]],
                [[1.0, 2.0, 3.0], [4.0, 0.0, 1.0], [4.0, 2.0, 1.0]],
            ],
        ],
    ]
}

fn get_tensor_rank_4_list() -> TensorRank4List<3, 1, 1, 1, 1, 2> {
    TensorRank4List::new(get_array())
}

#[test]
fn as_array() {
    get_tensor_rank_4_list()
        .as_array()
        .iter()
        .zip(get_array().iter())
        .for_each(|(tensor_rank_4_entry, array_entry)| {
            tensor_rank_4_entry.iter().zip(array_entry.iter()).for_each(
                |(tensor_rank_4_entry_i, array_entry_i)| {
                    tensor_rank_4_entry_i
                        .iter()
                        .zip(array_entry_i.iter())
                        .for_each(|(tensor_rank_4_entry_ij, array_entry_ij)| {
                            tensor_rank_4_entry_ij
                                .iter()
                                .zip(array_entry_ij.iter())
                                .for_each(|(tensor_rank_4_entry_ijk, array_entry_ijk)| {
                                    tensor_rank_4_entry_ijk
                                        .iter()
                                        .zip(array_entry_ijk.iter())
                                        .for_each(|(tensor_rank_4_entry_ijkl, array_entry_ijkl)| {
                                            assert_eq!(tensor_rank_4_entry_ijkl, array_entry_ijkl)
                                        })
                                })
                        })
                },
            )
        });
}

#[test]
fn from_iter() {
    let into_iterator = get_tensor_rank_4_list().0.into_iter();
    let tensor_rank_4_list =
        TensorRank4List::<3, 1, 1, 1, 1, 2>::from_iter(get_tensor_rank_4_list().0.into_iter());
    tensor_rank_4_list
        .iter()
        .zip(into_iterator)
        .for_each(|(tensor_rank_4_entry, array_entry)| {
            tensor_rank_4_entry.iter().zip(array_entry.iter()).for_each(
                |(tensor_rank_4_entry_i, array_entry_i)| {
                    tensor_rank_4_entry_i
                        .iter()
                        .zip(array_entry_i.iter())
                        .for_each(|(tensor_rank_4_entry_ij, array_entry_ij)| {
                            tensor_rank_4_entry_ij
                                .iter()
                                .zip(array_entry_ij.iter())
                                .for_each(|(tensor_rank_4_entry_ijk, array_entry_ijk)| {
                                    tensor_rank_4_entry_ijk
                                        .iter()
                                        .zip(array_entry_ijk.iter())
                                        .for_each(|(tensor_rank_4_entry_ijkl, array_entry_ijkl)| {
                                            assert_eq!(tensor_rank_4_entry_ijkl, array_entry_ijkl)
                                        })
                                })
                        })
                },
            )
        });
}

#[test]
fn iter() {
    get_tensor_rank_4_list()
        .iter()
        .zip(get_array().iter())
        .for_each(|(tensor_rank_4_entry, array_entry)| {
            tensor_rank_4_entry.iter().zip(array_entry.iter()).for_each(
                |(tensor_rank_4_entry_i, array_entry_i)| {
                    tensor_rank_4_entry_i
                        .iter()
                        .zip(array_entry_i.iter())
                        .for_each(|(tensor_rank_4_entry_ij, array_entry_ij)| {
                            tensor_rank_4_entry_ij
                                .iter()
                                .zip(array_entry_ij.iter())
                                .for_each(|(tensor_rank_4_entry_ijk, array_entry_ijk)| {
                                    tensor_rank_4_entry_ijk
                                        .iter()
                                        .zip(array_entry_ijk.iter())
                                        .for_each(|(tensor_rank_4_entry_ijkl, array_entry_ijkl)| {
                                            assert_eq!(tensor_rank_4_entry_ijkl, array_entry_ijkl)
                                        })
                                })
                        })
                },
            )
        });
}

#[test]
fn iter_mut() {
    get_tensor_rank_4_list()
        .iter_mut()
        .zip(get_array().iter_mut())
        .for_each(|(tensor_rank_4_entry, array_entry)| {
            tensor_rank_4_entry
                .iter_mut()
                .zip(array_entry.iter_mut())
                .for_each(|(tensor_rank_4_entry_i, array_entry_i)| {
                    tensor_rank_4_entry_i
                        .iter_mut()
                        .zip(array_entry_i.iter_mut())
                        .for_each(|(tensor_rank_4_entry_ij, array_entry_ij)| {
                            tensor_rank_4_entry_ij
                                .iter_mut()
                                .zip(array_entry_ij.iter_mut())
                                .for_each(|(tensor_rank_4_entry_ijk, array_entry_ijk)| {
                                    tensor_rank_4_entry_ijk
                                        .iter_mut()
                                        .zip(array_entry_ijk.iter_mut())
                                        .for_each(|(tensor_rank_4_entry_ijkl, array_entry_ijkl)| {
                                            assert_eq!(tensor_rank_4_entry_ijkl, array_entry_ijkl)
                                        })
                                })
                        })
                })
        });
}

#[test]
fn new() {
    get_tensor_rank_4_list()
        .iter()
        .zip(get_array().iter())
        .for_each(|(tensor_rank_4_entry, array_entry)| {
            tensor_rank_4_entry.iter().zip(array_entry.iter()).for_each(
                |(tensor_rank_4_entry_i, array_entry_i)| {
                    tensor_rank_4_entry_i
                        .iter()
                        .zip(array_entry_i.iter())
                        .for_each(|(tensor_rank_4_entry_ij, array_entry_ij)| {
                            tensor_rank_4_entry_ij
                                .iter()
                                .zip(array_entry_ij.iter())
                                .for_each(|(tensor_rank_4_entry_ijk, array_entry_ijk)| {
                                    tensor_rank_4_entry_ijk
                                        .iter()
                                        .zip(array_entry_ijk.iter())
                                        .for_each(|(tensor_rank_4_entry_ijkl, array_entry_ijkl)| {
                                            assert_eq!(tensor_rank_4_entry_ijkl, array_entry_ijkl)
                                        })
                                })
                        })
                },
            )
        });
}

#[test]
fn size() {
    assert_eq!(
        std::mem::size_of::<TensorRank4List::<3, 1, 1, 1, 1, 8>>(),
        std::mem::size_of::<[TensorRank4::<3, 1, 1, 1, 1>; 8]>()
    )
}

#[test]
fn zero() {
    TensorRank4List::<3, 1, 1, 1, 1, 8>::zero()
        .iter()
        .for_each(|tensor_rank_4_entry| {
            tensor_rank_4_entry
                .iter()
                .for_each(|tensor_rank_4_entry_i| {
                    tensor_rank_4_entry_i
                        .iter()
                        .for_each(|tensor_rank_4_entry_ij| {
                            tensor_rank_4_entry_ij
                                .iter()
                                .for_each(|tensor_rank_4_entry_ijk| {
                                    tensor_rank_4_entry_ijk.iter().for_each(
                                        |tensor_rank_4_entry_ijkl| {
                                            assert_eq!(tensor_rank_4_entry_ijkl, &0.0)
                                        },
                                    )
                                })
                        })
                })
        });
}
