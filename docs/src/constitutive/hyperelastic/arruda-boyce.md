# Arruda-Boyce model

**Parameters**

- The bulk modulus \\(\\kappa\\).
- The number density of chains \\(n\\).
- The number of links \\(N_b\\).

**State variables**

- The deformation gradient \\(\\mathbf{F}\\).
- The temperature \\(T\\).

**Notes**

- The shear modulus in the Arruda-Boyce model is given by \\(\mu=nkT\\).
- The Arruda-Boyce model reduces to the [Neo-Hookean model](neo-hookean.md) when \\(N_b\to\infty\\).

## Helmholtz free energy density

\\begin{equation}
    a = ? + \\frac{\\kappa}{2}\\left[\\frac{1}{2}\\left(J^2 - 1\\right) - \\ln J\\right]
\\end{equation}

## Cauchy stress

\\begin{equation}
    \\boldsymbol{\\sigma} = ? + \\frac{\\kappa}{2}\\left(J - \\frac{1}{J}\\right)\\mathbf{1}
\\end{equation}

## Tangent stiffness

\\begin{equation}
    \\mathcal{T}\_{ijkL} = ?
\\end{equation}
