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
        TensorRank1Trait,
        TensorRank1List,
        TensorRank1ListTrait,
        TensorRank1List2D,
        TensorRank1List2DTrait,
        TensorRank2,
        TensorRank2List,
        TensorRank2Trait,
        TensorRank2ListTrait,
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
        FirstPiolaKirchoffTangentStiffnesses,
        FirstPiolaKirchoffRateTangentStiffnesses,
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

type Basis<const I: usize> = Vectors<I, 2>;
type Connectivity<const E: usize, const N: usize> = [[usize; N]; E];
type GradientVectors<const N: usize> = Vectors<0, N>;
type NodalCoordinates<const D: usize> = CurrentCoordinates<D>;
type NodalForces<const D: usize> = Forces<D>;
type NodalStiffnesses<const D: usize> = Stiffnesses<D>;
type NodalVelocities<const D: usize> = CurrentCoordinates<D>;
type Normal<const I: usize> = Vector<I>;
type NormalGradients<const O: usize> = TensorRank2List<3, 1, 1, O>;
type NormalTangents<const O: usize> = TensorRank3List2D<3, 1, 1, 1, O>;
type NormalRate = Vector<1>;
type ParametricGradientOperators<const P: usize> = TensorRank2List<3, 0, 9, P>;
type ProjectedGradientVectors<const G: usize, const N: usize> = Vectors2D<0, N, G>;
type ProjectionMatrix<const Q: usize> = TensorRank2<Q, 9, 9>;
type ReferenceNodalCoordinates<const D: usize> = ReferenceCoordinates<D>;
type ReferenceNormal = Vector<0>;
type ShapeFunctionIntegrals<const P: usize, const Q: usize> = TensorRank1List<Q, 9, P>;
type ShapeFunctionIntegralsProducts<const P: usize, const Q: usize> = TensorRank2List<Q, 9, 9, P>;
type ShapeFunctionsAtIntegrationPoints<const G: usize, const Q: usize> = TensorRank1List<Q, 9, G>;
type StandardGradientOperator<const M: usize, const O: usize> = TensorRank1List<M, 9, O>;
type StandardGradientOperators<const M: usize, const O: usize, const P: usize> = TensorRank1List2D<M, 9, O, P>;
type StandardGradientOperatorsTransposed<const M: usize, const O: usize, const P: usize> = TensorRank1List2D<M, 9, P, O>;