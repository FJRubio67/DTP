# DTP: Double Two-Piece Distributions (R package)

[![R package](https://img.shields.io/badge/language-R-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview

`DTP` is an R package implementing the family of **double two-piece (DTP)
distributions** — a flexible class of unimodal distributions that extends the
standard two-piece construction by allowing the tail behaviour to differ
independently on each side of the mode.

While a standard two-piece distribution joins two rescaled half-densities at
the mode (controlling asymmetry through a single skewness parameter), the DTP
family joins two **shape-modified** half-densities, each with its own shape
parameter $\delta_1$ and $\delta_2$. This allows the left and right tails to
have different degrees of heaviness, making the DTP family capable of
modelling both asymmetry and differential tail behaviour simultaneously — a
feature not available in standard two-piece or symmetric unimodal families.

The DTP density takes the form:

$$f_{\mathrm{DTP}}(x;\, \mu, \sigma, \gamma, \delta_1, \delta_2) =
w \cdot f_{\mathrm{L}}(x;\, \mu, \sigma, \delta_1) +
(1-w) \cdot f_{\mathrm{R}}(x;\, \mu, \sigma, \delta_2)$$

where $f_{\mathrm{L}}$ and $f_{\mathrm{R}}$ are rescaled left and right
half-densities derived from a symmetric baseline $f$, and $w$ is a normalising
weight. The baseline distribution is passed as a function argument, making the
package compatible with any symmetric unimodal density available in R.

## Installation

```r
# install.packages("devtools")
devtools::install_github("FJRubio67/DTP")
library(DTP)
```

## Main functions

| Function | Description |
|---|---|
| `ddtp` | Probability density function |
| `pdtp` | Cumulative distribution function |
| `qdtp` | Quantile function |
| `rdtp` | Random number generation |

The baseline distribution is passed via the `f` argument (e.g. `f = dt`,
`f = dnorm`). Two parameterisations of skewness are supported via the `param`
argument: `"tp"` (two-piece scale) and `"eps"` (epsilon-skew).

For full documentation: `?ddtp`

## Quick example

```r
library(DTP)

# Double two-piece t with different left and right tail weights
x <- seq(-6, 6, length.out = 500)
plot(x, ddtp(x, mu = 0, sigma = 1, gamma = 0.3,
             delta1 = 2, delta2 = 8, f = dt, param = "eps"),
     type = "l", ylab = "Density", main = "DTP-t distribution")
```

## Tutorials

- [DTP package](https://rpubs.com/FJRubio/DTP) — illustrative walkthrough with
  examples and parameter interpretation
- [The double two-piece sinh-arcsinh distribution](https://rpubs.com/FJRubio/DTPSAS) —
  DTP applied to the sinh-arcsinh baseline
- [Galton's forecasting competition](https://rpubs.com/FJRubio/Galton) —
  real-data example fitting DTP-t to Galton's (1907) ox weight data

## Related packages

- [twopiece](https://github.com/FJRubio67/twopiece) — the standard two-piece
  family (single shape parameter per side); parent construction
- [TPSAS](https://github.com/FJRubio67/TPSAS) — two-piece sinh-arcsinh
  distributions; related family with joint skewness and tail control

## License

This package is licensed under the [MIT License](LICENSE).
