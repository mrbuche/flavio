# Neo-Hookean model

**Parameters**

- The bulk modulus \\(\\kappa\\).
- The shear modulus \\(\\mu\\).

**State variables**

- The deformation gradient \\(\\mathbf{F}\\).

## Helmholtz free energy density

\\begin{equation}
    a = \\frac{\\mu}{2}\\left[\\mathrm{tr}(\\mathbf{B}^*) - 3\\right] + \\frac{\\kappa}{2}\\left[\\frac{1}{2}\\left(J^2 - 1\\right) - \\ln J\\right]
\\end{equation}

## Cauchy stress

\\begin{equation}
    \\boldsymbol{\\sigma} = \\frac{\\mu}{J}\\,{\\mathbf{B}^*}' + \\frac{\\kappa}{2}\\left(J - \\frac{1}{J}\\right)\\mathbf{1}
\end{equation}

## Tangent stiffness

\\begin{equation}
    \\mathcal{T}\_{ijkL} = \\frac{\\mu}{J^{5/3}}\\left(\\delta\_{ik}F\_{jL} + \\delta\_{jk}F\_{iL} - \\frac{2}{3}\\,\\delta\_{ij}F\_{kL}- \\frac{5}{3} \\, B\_{ij}'F\_{kL}^{-T} \\right) + \\frac{\\kappa}{2} \\left(J + \\frac{1}{J}\\right)\\delta\_{ij}F\_{kL}^{-T}
\\end{equation}
