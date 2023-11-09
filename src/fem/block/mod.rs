#[cfg(test)]
mod test;

pub mod element;

use crate::
{
    constitutive::
    {
        ConstitutiveModel,
        ConstitutiveModelParameters
    },
    math::
    {
        ContractSecondFourthIndicesWithFirstIndicesOf,
        TensorRank0ListTrait,
        TensorRank1ListTrait,
        TensorRank2Trait
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
use self::element::FiniteElement;

type CurrentNodalCoordinates<const D: usize> = CurrentCoordinates<D>;
type NodalForces<const D: usize> = Forces<D>;
type NodalStiffnesses<const D: usize> = Stiffnesses<D>;
type ReferenceNodalCoordinates<const D: usize> = ReferenceCoordinates<D>;

pub struct FiniteElementBlock<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
where
    C: ConstitutiveModel<'a>,
    F: FiniteElement<'a, C, G, N>
{
    elements: [F; E],
    phantom_a: std::marker::PhantomData<*const &'a C>
}

pub trait FiniteElementBlockTraits<'a, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
where
    C: ConstitutiveModel<'a>,
    F: FiniteElement<'a, C, G, N>
{
    fn calculate_helmholtz_free_energy(&self) -> Scalar;
    fn calculate_nodal_forces(&self) -> NodalForces<D>;
    fn calculate_nodal_stiffnesses(&self) -> NodalStiffnesses<D>;
    fn get_current_nodal_coordinates(&self) -> &CurrentNodalCoordinates<D>;
    fn get_elements(&self) -> &[F; E];
    fn set_current_nodal_coordinates(&mut self, current_nodal_coordinates: CurrentNodalCoordinates<D>);
}