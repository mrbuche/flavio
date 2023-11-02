# Arruda-Boyce model

**Parameters**

- The bulk modulus \\(\\kappa\\).
- The shear modulus \\(\\mu\\).
- The number of links \\(N_b\\).

**State variables**

- The deformation gradient \\(\\mathbf{F}\\).

**Notes**

- The nondimensional end-to-end length per link of the chains is \\(\\gamma=\\sqrt{\\mathrm{tr}(\\mathbf{B}^*)/3N_b}\\).
- The nondimensional force is given by the inverse Langevin function as \\(\\eta(\\gamma)=\\mathcal{L}^{-1}(\\gamma)\\).
- The Arruda-Boyce model reduces to the [Neo-Hookean model](neo-hookean.md) when \\(N_b\to\infty\\).

## Helmholtz free energy density

\\begin{equation}
    a = \\mu N_b\\left[\\gamma\\eta - \\ln\\left(\\frac{\\sinh\\eta}{\\eta}\\right)\\right] + \\frac{\\kappa}{2}\\left[\\frac{1}{2}\\left(J^2 - 1\\right) - \\ln J\\right]
\\end{equation}

## Cauchy stress

\\begin{equation}
    \\boldsymbol{\\sigma} = ? + \\frac{\\kappa}{2}\\left(J - \\frac{1}{J}\\right)\\mathbf{1}
\\end{equation}

## Tangent stiffness

\\begin{equation}
    \\mathcal{T}\_{ijkL} = ?
\\end{equation}
