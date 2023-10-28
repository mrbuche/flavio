# Mooney-Rivlin model

**Parameters**

- The bulk modulus \\(\\kappa\\).
- The shear modulus \\(\\mu\\).
- The extra modulus \\(\\mu_m\\).

**State variables**

- The deformation gradient \\(\\mathbf{F}\\).

**Notes**

- The Mooney-Rivlin model reduces to the [Neo-Hookean model](neo-hookean.md) when \\(\\mu_m\to 0\\).

## Helmholtz free energy density

\\begin{equation}
    a = \\frac{\\mu - \\mu_m}{2}\\left[\\mathrm{tr}(\\mathbf{B}^* ) - 3\\right] + \\frac{\\mu_m}{2}\\left[I_2(\\mathbf{B}^*) - 3\\right] + \\frac{\\kappa}{2}\\left[\\frac{1}{2}\\left(J^2 - 1\\right) - \\ln J\\right]
\\end{equation}

## Cauchy stress

\\begin{equation}
    \\boldsymbol{\\sigma} = \\frac{\\mu - \\mu_m}{J}\\, {\\mathbf{B}^* }' - \\frac{\\mu_m}{J}\\left(\\mathbf{B}^{* -1}\\right)' + \\frac{\\kappa}{2}\\left(J - \\frac{1}{J}\\right)\\mathbf{1}
\\end{equation}

## Tangent stiffness

\\begin{align}
    \\mathcal{T}\_{ijkL} = & \\frac{\\mu-\\mu_m}{J^{5/3}}\\left(\\delta\_{ik}F\_{jL} + \\delta\_{jk}F\_{iL} - \\frac{2}{3}\\,\\delta\_{ij}F\_{kL}- \\frac{5}{3} \\, B\_{ij}'F\_{kL}^{-T} \\right) + \\frac{\\kappa}{2} \\left(J + \\frac{1}{J}\\right)\\delta\_{ij}F\_{kL}^{-T} \\\\ & - \\frac{\\mu_m}{J}\\left[ \\frac{2}{3}\\,B\_{ij}^{* -1}F\_{kL}^{-T} - B\_{ik}^{* -1}F\_{jL}^{-T} - B\_{ik}^{* -1}F\_{iL}^{-T} + \\frac{2}{3}\\,\\delta\_{ij}\\left(B\_{km}^{* -1}\\right)'F\_{mL}^{-T} - \\left(B\_{ij}^{* -1}\\right)'F\_{kL}^{-T} \\right]
\\end{align}
