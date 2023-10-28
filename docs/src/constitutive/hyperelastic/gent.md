# Gent model

**Parameters**

- The bulk modulus \\(\\kappa\\).
- The shear modulus \\(\\mu\\).
- The extensibility \\(J_m\\).

**State variables**

- The deformation gradient \\(\\mathbf{F}\\).

**Notes**

- The Gent model reduces to the [Neo-Hookean model](neo-hookean.md) when \\(J_m\to\infty\\).

## Helmholtz free energy density

\begin{equation}
    a = -\\frac{\\mu J_m}{2}\\,\\ln\\left[1 - \\frac{\\mathrm{tr}(\\mathbf{B}^* ) - 3}{J_m}\\right] + \\frac{\\kappa}{2}\\left[\\frac{1}{2}\\left(J^2 - 1\\right) - \\ln J\\right]
\end{equation}

## Cauchy stress

\begin{equation}
    \\boldsymbol{\\sigma} = \\frac{J^{-1}\\mu J_m {\\mathbf{B}^* }'}{J_m - \\mathrm{tr}(\\mathbf{B}^* ) + 3} + \\frac{\\kappa}{2}\\left(J - \\frac{1}{J}\\right)\\mathbf{1}
\end{equation}

## Tangent stiffness

\begin{align}
    \\mathcal{T}\_{ijkL} = \\frac{J^{-5/3}\\mu J_m}{J_m - \\mathrm{tr}(\\mathbf{B}^* ) + 3}\\Bigg[ & \\delta\_{ik}F\_{jL} + \\delta\_{jk}F\_{iL} - \\frac{2}{3}\\,\\delta\_{ij}F\_{kL} + \\frac{2{B\_{ij}^* }' F\_{kL}}{J_m - \\mathrm{tr}(\\mathbf{B}^* ) + 3} \\\\ & - \\left(\\frac{5}{3} + \\frac{2}{3}\\frac{\\mathrm{tr}(\\mathbf{B}^* )}{J_m - \\mathrm{tr}(\\mathbf{B}^* ) + 3}\\right) J^{2/3} {B\_{ij}^* }' F\_{kL}^{-T} \\Bigg] + \\frac{\\kappa}{2} \\left(J + \\frac{1}{J}\\right)\\delta\_{ij}F\_{kL}^{-T}
\end{align}
