#[cfg(test)]
mod test;

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
    ABS_TOL,
    REL_TOL,
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

pub trait FiniteElementModel<'a, B, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
where
    B: FiniteElementBlock<'a, C, D, E, F, G, N>,
    C: ConstitutiveModel<'a> + HyperelasticConstitutiveModel,
    F: FiniteElement<'a, C, G, N>
{
    // going to get real tricky to keep static besides 1 block 1 element 1 constitutive
    // because need to hold a list of Blocks which have generic arguments
    //
    // it's the generic arguments of the trait that matters
    // 'a is OK
    // C comes from F which comes from get_elements()
    // D, G, N would also need to go and that seems too hard
    //
    // the alternative would be to implement every possible combination of blocks so each has own generic arguments
    // might be useful for things like hyperelastic mixed with viscoelastic blocks
    //
    // maybe do that here for now? starting with 1 hyperelastic block?
    // and move the solver/BCs here? (not "owned" by blocks)
    //
    // or should solvers stay in block?
    // so that can solve different depending on constitutive model in the block?
    // and have the block solvers take in arguments that would represent BCs or merges to other blocks?
    //
    // not really solving differently, just calling nodal_forces and stuff which already work differently
    // so the solver top-level stuff can stay here?
    //
    // so now the model needs to own the coordinates, since blocks will also share some
    // call the block routines (not elements) since blocks could have different constitutive models (some simpler than others)
    // but remember:
    //     nodes are not owned by blocks
    //     ELEMENTS ARE OWNED BY BLOCKS
    //
    // so then the connectivity (list of nodes for each element) is owned by each block
    // and can stay in Block along with the list of elements (those two things really define a block)
    //
    fn calculate_nodal_forces(&self) -> NodalForces<D>;
    fn get_block(&self) -> &B;
    fn get_current_nodal_coordinates(&self) -> &CurrentNodalCoordinates<D>;
    fn set_current_nodal_coordinates(&mut self, current_nodal_coordinates: CurrentNodalCoordinates<D>);
    fn solve_using_gradient_descent(&mut self);
}

pub trait HyperelasticFiniteElementModel<'a, B, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
where
    B: FiniteElementBlock<'a, C, D, E, F, G, N>,
    C: ConstitutiveModel<'a> + HyperelasticConstitutiveModel,
    F: FiniteElement<'a, C, G, N>,
    Self: FiniteElementModel<'a, B, C, D, E, F, G, N>
{
    fn calculate_helmholtz_free_energy(&self) -> Scalar;
}

pub struct Model<'a, B, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
where
    B: FiniteElementBlock<'a, C, D, E, F, G, N>,
    C: ConstitutiveModel<'a> + HyperelasticConstitutiveModel,
    F: FiniteElement<'a, C, G, N>
{
    block: B,
    phantom_a: std::marker::PhantomData<*const &'a C>,
    phantom_f: std::marker::PhantomData<F>
}

impl<'a, B, C, const D: usize, const E: usize, F, const G: usize, const N: usize>
    FiniteElementModel<'a, B, C, D, E, F, G, N>
    for Model<'a, B, C, D, E, F, G, N>
where
    B: FiniteElementBlock<'a, C, D, E, F, G, N>,
    C: ConstitutiveModel<'a> + HyperelasticConstitutiveModel,
    F: FiniteElement<'a, C, G, N>
{
    fn calculate_nodal_forces(&self) -> NodalForces<D>
    {
        self.get_block().calculate_nodal_forces()
    }
    fn get_block(&self) -> &B
    {
        &self.block
    }
    fn get_current_nodal_coordinates(&self) -> &CurrentNodalCoordinates<D>
    {
        self.get_block().get_current_nodal_coordinates()
    }
    fn set_current_nodal_coordinates(&mut self, current_nodal_coordinates: CurrentNodalCoordinates<D>)
    {
        todo!()
        // self.get_block().current_nodal_coordinates = current_nodal_coordinates;
    }
    fn solve_using_gradient_descent(&mut self)
    {
        todo!()
        // let mut step_size = REL_TOL;
        // let mut coordinates = self.get_current_nodal_coordinates() * 1.0;
        // let mut nodal_forces = self.calculate_nodal_forces();
        // let mut residual = nodal_forces.norm();
        // let mut coordinates_old: CurrentNodalCoordinates<D>;
        // let mut nodal_forces_old: NodalForces<D>;
        // let mut nodal_forces_diff: NodalForces<D>;
        // while residual > ABS_TOL
        // {
        //     coordinates_old = &coordinates * 1.0;
        //     nodal_forces_old = &nodal_forces * 1.0;
        //     self.set_current_nodal_coordinates(
        //         nodal_forces * (-step_size) + coordinates
        //     );
        //     coordinates = self.get_current_nodal_coordinates() * 1.0;
        //     nodal_forces = self.calculate_nodal_forces();
        //     nodal_forces_diff = nodal_forces_old - &nodal_forces;
        //     step_size = (
        //         (coordinates_old - &coordinates).dot(&nodal_forces_diff)
        //     ).abs() / nodal_forces_diff.dot_self();
        //     residual = nodal_forces.norm();
        // }
    }
}