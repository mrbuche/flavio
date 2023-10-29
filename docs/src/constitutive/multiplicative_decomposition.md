# Multiplicative decomposition

Also useful to write down. Sometimes called Kr&ouml;ner-Lee decomposition, usually when referring to elastic-plastic decompositions.

Improve this, make more general (second constitutive model usually does not have stress directly specified, the Mandel stress is what is useful), and no need to go into other residuals if not using them.

\\begin{equation}
    \\mathbf{F} =
    \\mathbf{F}\_1\cdot\\mathbf{F}\_2
\\end{equation}

The Helmholtz free energy density.

\\begin{equation}
    a(\\mathbf{F},\\mathbf{F}\_2) =
    a\_1(\\mathbf{F}\_1) + a\_2(\\mathbf{F}\_2)
\\end{equation}

The Cauchy stress.

\\begin{equation}
    \\boldsymbol{\\sigma}(\\mathbf{F},\\mathbf{F}\_2) =
    J_2^{-1} \\boldsymbol{\\sigma}_1(\\mathbf{F}\_1)
\\end{equation}

The tangent stiffness.

\\begin{equation}
    \\boldsymbol{\\mathcal{T}}(\\mathbf{F},\\mathbf{F}\_2) =
    J_2^{-1} \\boldsymbol{\\mathcal{T}}_1(\\mathbf{F}\_1) \\cdot \\mathbf{F}\_2^{-T}
\\end{equation}

The Helmholtz free energy density is minimized with respect to \\(\\mathbf{F}\_2\\) and simplied to show that Mandel stress in the first consitutive model equals the Kirchoff stress in the second.

\\begin{equation}
    \\mathbf{M}_1 =
    \\mathbf{F}\_1^T\\cdot(J\_1\\boldsymbol{\\sigma}\_1)\\cdot\\mathbf{F}\_1^{-T} =
    \\mathbf{F}\_1^T\\cdot\\mathbf{P}\_1 =
    J\_2 \\boldsymbol{\\sigma}\_2
\\end{equation}

When solving for the \\(\\mathbf{F}\_2\\) that minimizes the Helmholtz free energy under the multiplicative decomposition constraint, the residual can then be defined as the difference between the Mandel stress from the first model and the Kirchoff stress from the second model.

\\begin{equation}
    \\mathbf{R}(\\mathbf{F},\\mathbf{F}\_2) =
    J\_2 \\boldsymbol{\\sigma}\_2 - \\mathbf{F}\_1^T\\cdot\\mathbf{P}\_1
\\end{equation}

Note that this residual is indirectly related to the derivative of the Helmholtz free energy.

\\begin{equation}
    \\frac{\\partial a}{\\partial\\mathbf{F}\_2} =
    \\mathbf{P}\_2 - \\mathbf{F}\_1^T\\cdot\\mathbf{P}\_1\\cdot\\mathbf{F}\_2^{-T}
\\end{equation}

The residual is equivalent to using the direct residual as long as \\(\\mathbf{F}\_2\\) is invertible.

\\begin{equation}
    \\frac{\\partial a}{\\partial\\mathbf{F}\_2}\\cdot\\mathbf{F}\_2^T =
    \\mathbf{R}
\\end{equation}

The tangent associated with this residual can be calculated using another derivative.

\\begin{equation}
    \\frac{\\partial R\_{i'j'}}{\\partial F\_{2,k'L}} =
    J\_2 \\sigma\_{2,i'j'} F^{-T}\_{2,k'L} + J\_2 \\mathcal{T}\_{2,i'j'k'L} + F^{-T}\_{2,i'L} M\_{1,k'j'} + \\mathcal{C}\_{1,sj'mn'} F\_{1,si'} F\_{1,mk'} F^{-T}\_{2,n'L}
\\end{equation}

The other tangent is similarly calculated.

\\begin{equation}
    \\frac{\\partial^2 a}{\\partial F\_{2,i'J}\\partial F\_{2,k'L}} =
    \\mathcal{C}\_{2,i'Jk'L} + \\mathcal{C}\_{2,mn'op'}F\_{1,mi'}F\_{2,n'J}^{-T}F\_{1,ok'}F\_{2,p'L}^{-T} + F\_{2,i'L}^{-T}F\_{1,k'm}^TP\_{1,mn'}F\_{2,n'J}^{-T} + F\_{1,i'm}^TP\_{1,mn'}F\_{2,n'L}^{-T}F\_{2,k'J}^{-T}
\\end{equation}

Note that indices in the intermediate configuartion are denoted with apostrophes.
