# Yeoh model

**Parameters**

- The bulk modulus \\(\\kappa\\).
- The shear modulus \\(\\mu=\\mu\_1\\).
- The extra moduli \\(\\mu\_n\\) for \\(n=2\ldots N\\).

**State variables**

- The deformation gradient \\(\\mathbf{F}\\).

**Notes**

- The Yeoh model reduces to the [Neo-Hookean model](neo-hookean.md) when \\(\\mu\_n\\to 0\\) for \\(n=2\ldots N\\).

## Helmholtz free energy density

\\begin{equation}
    a = \\sum\_{n=1}^N \\frac{\\mu_n}{2}\\left[\\mathrm{tr}(\\mathbf{B}^* ) - 3\\right]^n + \\frac{\\kappa}{2}\\left[\\frac{1}{2}\\left(J^2 - 1\\right) - \\ln J\\right]
\\end{equation}

## Cauchy stress

\\begin{equation}
    \\boldsymbol{\\sigma} = \\sum\_{n=1}^N \\frac{n\\mu_n}{J}\\left[\\mathrm{tr}(\\mathbf{B}^* ) - 3\\right]^{n-1}\\,{\\mathbf{B}^*}' + \\frac{\\kappa}{2}\\left(J - \\frac{1}{J}\\right)\\mathbf{1}
\\end{equation}

## Tangent stiffness

\\begin{align}
    \\mathcal{T}\_{ijkL} = & \\sum\_{n=1}^N \\frac{n\\mu_n}{J^{5/3}}\\left[\\mathrm{tr}(\\mathbf{B}^* ) - 3\\right]^{n-1}\\left(\\delta\_{ik}F\_{jL} + \\delta\_{jk}F\_{iL} - \\frac{2}{3}\\,\\delta\_{ij}F\_{kL}- \\frac{5}{3} \\, B\_{ij}'F\_{kL}^{-T} \\right) \\\\ & + \\sum\_{n=2}^N \\frac{2n(n-1)\\mu_n}{J^{7/3}}\\left[\\mathrm{tr}(\\mathbf{B}^* ) - 3\\right]^{n-2}B\_{ij}'B\_{km}'F\_{mL}^{-T} + \\frac{\\kappa}{2} \\left(J + \\frac{1}{J}\\right)\\delta\_{ij}F\_{kL}^{-T}
\\end{align}
