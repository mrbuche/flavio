mod block;

pub use block::
{
    Block,
    FiniteElementBlock,
    HyperelasticFiniteElementBlock,
    element::
    {
        FiniteElement,
        HyperelasticFiniteElement,
        linear::
        {
            LinearFiniteElement,
            HyperelasticLinearFiniteElement,
            tetrahedron::
            {
                LinearTetrahedron
            }
        }
    }
};

use crate::
{
    constitutive::
    {
        ConstitutiveModel,
        ConstitutiveModelParameters,
        hyperelastic::HyperelasticConstitutiveModel
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
type CurrentNodalCoordinates<const D: usize> = CurrentCoordinates<D>;
type NodalForces<const D: usize> = Forces<D>;
type NodalStiffnesses<const D: usize> = Stiffnesses<D>;
type ReferenceNodalCoordinates<const D: usize> = ReferenceCoordinates<D>;
