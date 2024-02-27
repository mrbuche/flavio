//! Finite element library.

mod block;

pub use block::
{
    ElasticBlock,
    ViscoelasticBlock,
    FiniteElementBlock,
    ElasticFiniteElementBlock,
    HyperelasticFiniteElementBlock,
    ViscoelasticFiniteElementBlock,
    HyperviscoelasticFiniteElementBlock,
    element::
    {
        FiniteElement,
        ElasticFiniteElement,
        HyperelasticFiniteElement,
        ViscoelasticFiniteElement,
        HyperviscoelasticFiniteElement,
        linear::
        {
            LinearFiniteElement,
            HyperelasticLinearFiniteElement,
            ViscoelasticLinearFiniteElement,
            HyperviscoelasticLinearFiniteElement,
            tetrahedron::
            {
                Tetrahedron
            }
        }
    }
};

use crate::
{
    constitutive::
    {
        Constitutive,
        Parameters,
        solid::
        {
            Solid,
            elastic::Elastic,
            hyperelastic::Hyperelastic,
            viscoelastic::Viscoelastic,
            elastic_hyperviscous::ElasticHyperviscous,
            hyperviscoelastic::Hyperviscoelastic
        }
    },
    math::
    {
        ContractSecondFourthIndicesWithFirstIndicesOf,
        Convert,
        TensorRank0ListTrait,
        TensorRank1ListTrait,
        TensorRank2Trait,
        TensorRank2List2DTrait
    },
    mechanics::
    {
        CurrentCoordinates,
        DeformationGradient,
        DeformationGradientRate,
        FirstPiolaKirchoffStress,
        FirstPiolaKirchoffTangentStiffness,
        Forces,
        ReferenceCoordinates,
        Scalar,
        Scalars,
        Stiffnesses,
        Vectors
    }
};

type Connectivity<const E: usize, const N: usize> = [[usize; N]; E];
type NodalCoordinates<const D: usize> = CurrentCoordinates<D>;
type NodalForces<const D: usize> = Forces<D>;
type NodalStiffnesses<const D: usize> = Stiffnesses<D>;
type NodalVelocities<const D: usize> = CurrentCoordinates<D>;
type ReferenceNodalCoordinates<const D: usize> = ReferenceCoordinates<D>;
