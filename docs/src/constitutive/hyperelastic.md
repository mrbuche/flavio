# Hyperelastic models

- [Gent model](hyperelastic/gent.md)
- [Neo-Hookean model](hyperelastic/neo-hookean.md)
- [Yeoh model](hyperelastic/yeoh.md)

## Helmholtz free energy density

A hyperelastic constitutive model can be defined by specifying the Helmholtz free energy density as some function of the deformation gradient

\\begin{equation}
    a = a(\\mathbf{F})
\\end{equation}

## Cauchy stress

The Cauchy stress is calculated from the Helmholtz free energy density as

\\begin{equation}
    \\boldsymbol{\\sigma} = \\frac{1}{J}\\frac{\\partial a}{\\partial\\mathbf{F}}\\cdot\\mathbf{F}^T
\\end{equation}

## Tangent stiffness

The tangent stiffness is calculated from the Cauchy stress as

\\begin{equation}
    \\boldsymbol{\\mathcal{T}} = \\frac{\\partial\\boldsymbol{\\sigma}}{\\partial\\mathbf{F}}
\\end{equation}
