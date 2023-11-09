mod block;

pub use block::
{
    FiniteElementBlock,
    FiniteElementBlockTrait,
    element::
    {
        FiniteElement,
        linear::
        {
            LinearFiniteElement,
            tetrahedron::
            {
                LinearTetrahedron
            }
        }
    }
};
