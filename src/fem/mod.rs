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
            LinearElement,
            HyperelasticLinearElement,
            ViscoelasticLinearElement,
            HyperviscoelasticLinearElement,
            tetrahedron::
            {
                Tetrahedron
            },
            localization::
            {
                LinearLocalizationElement,
                wedge::Wedge
            },
            surface::
            {
                LinearSurfaceElement,
                triangle::Triangle
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
        TensorRank1Trait,
        TensorRank1List,
        TensorRank1ListTrait,
        TensorRank2,
        TensorRank2List,
        TensorRank2Trait,
        TensorRank2List2DTrait,
        TensorRank3List2D,
        levi_civita
    },
    mechanics::
    {
        Coordinates,
        CurrentCoordinates,
        DeformationGradient,
        DeformationGradients,
        DeformationGradientRate,
        DeformationGradientRates,
        FirstPiolaKirchoffStresses,
        Forces,
        ReferenceCoordinates,
        Scalar,
        Scalars,
        Stiffnesses,
        Vector,
        Vectors,
        Vectors2D
    }
};

type Connectivity<const E: usize, const N: usize> = [[usize; N]; E];
type IntegrationWeights<const G: usize> = Scalars<G>;
type NodalCoordinates<const D: usize> = CurrentCoordinates<D>;
type NodalForces<const D: usize> = Forces<D>;
type NodalStiffnesses<const D: usize> = Stiffnesses<D>;
type NodalVelocities<const D: usize> = CurrentCoordinates<D>;
type ReferenceNodalCoordinates<const D: usize> = ReferenceCoordinates<D>;
