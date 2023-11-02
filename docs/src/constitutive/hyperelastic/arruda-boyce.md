# Arruda-Boyce model

**Parameters**

- The bulk modulus \\(\\kappa\\).
- The shear modulus \\(\\mu\\).
- The number of links \\(N\_b\\).

**State variables**

- The deformation gradient \\(\\mathbf{F}\\).

**Notes**

- The nondimensional end-to-end length per link of the chains is \\(\\gamma=\\sqrt{\\mathrm{tr}(\\mathbf{B}^*)/3N\_b}\\).
- The nondimensional force is given by the inverse Langevin function as \\(\\eta=\\mathcal{L}^{-1}(\\gamma)\\).
- The initial values are given by \\(\\gamma\_0=\\sqrt{1/3N\_b}\\) and \\(\\eta\_0=\\mathcal{L}^{-1}(\\gamma\_0)\\).
- The Arruda-Boyce model reduces to the [Neo-Hookean model](neo-hookean.md) when \\(N\_b\to\infty\\).

## Helmholtz free energy density

\\begin{equation}
    a = \\mu N\_b\\left[\\gamma\\eta - \\gamma_0\\eta_0 - \\ln\\left(\\frac{\\eta_0\\sinh\\eta}{\\eta\\sinh\\eta_0}\\right) \\right] + \\frac{\\kappa}{2}\\left[\\frac{1}{2}\\left(J^2 - 1\\right) - \\ln J\\right]
\\end{equation}

## Cauchy stress

\\begin{equation}
    \\boldsymbol{\\sigma} = \\frac{\\mu\\eta}{3J\\gamma}\\,{\\mathbf{B}^*}' + \\frac{\\kappa}{2}\\left(J - \\frac{1}{J}\\right)\\mathbf{1}
\\end{equation}

## Tangent stiffness

\\begin{align}
    \\mathcal{T}\_{ijkL} = & \\frac{\\mu\\eta}{3J^{5/3}\\gamma}\\left(\\delta\_{ik}F\_{jL} + \\delta\_{jk}F\_{iL} - \\frac{2}{3}\\,\\delta\_{ij}F\_{kL}- \\frac{5}{3} \\, B\_{ij}'F\_{kL}^{-T} \\right) + \\frac{\\kappa}{2} \\left(J + \\frac{1}{J}\\right)\\delta\_{ij}F\_{kL}^{-T} \\\\ & + \\frac{\\mu\\eta}{9J^{7/3}N\_b\\gamma^2}\\left(\\frac{1}{\\eta\\mathcal{L}'(\\eta)} - \\frac{1}{\\gamma}\\right)B\_{ij}'B\_{km}'F\_{mL}^{-T}
\\end{align}
