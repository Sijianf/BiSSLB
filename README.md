
# BiSSLB <img src="https://img.shields.io/badge/R-package-blue.svg" alt="R package badge" height="20"/>

**BiSSLB**: Binary Spike-and-Slab Lasso Biclustering for binary matrices

BiSSLB is an R package for performing sparse biclustering on binary matrices. It uses a Bayesian spike-and-slab prior on matrix factorizations to reveal interpretable biclusters. The core computations are accelerated using `RcppArmadillo`.

---

## 🔧 Installation

To install the development version of the package from GitHub:

```
# Install devtools if not already installed
install.packages("devtools")

# Install BiSSLB from GitHub
devtools::install_github("sijianf/BiSSLB")
```

---

## 📦 Package Overview

The BiSSLB package includes the following key functions:

| Function              | Description                                                         |
|-----------------------|---------------------------------------------------------------------|
| `BiSSLB()`            | Main interface to run biclustering on a binary matrix              |
| `BiSSLB_ladder()`     | Run `BiSSLB` across a grid of spike prior parameters               |
| `main_iterations()`   | Core C++ engine performing iterative updates (not user-facing)     |
| `get_thetas()`        | Estimate sparsity parameters under spike-and-slab prior            |
| `get_logLikelihood()` | Compute the log-likelihood and BIC of current model state       |

---

## 🚀 Quick Start

```r
library(BiSSLB)

# Simulate a binary matrix
set.seed(123)
Y <- matrix(rbinom(100, 1, 0.4), nrow = 10)

# Run BiSSLB with default settings
result <- BiSSLB(Y, K_init = 5)

# Inspect output
names(result)

# View estimated row/column profiles
image(result$A)
image(result$B)
```

---

## 🔁 Parameter Tuning Example

Use `BiSSLB_ladder()` to explore multiple regularization settings:

```r
out_ladder <- BiSSLB_ladder(
  Y,
  tilde_lambda_0s = c(5, 10),
  lambda_0s = c(5, 10),
  tilde_lambda_1 = 1,
  lambda_1 = 1,
  K_init = 10
)

# Access results and BICs
str(out_ladder$results)
out_ladder$BICs
```

---

## 🧠 Method Summary

BiSSLB models a binary matrix Y as:
$P(Y_{ij} = 1) = sigmoid((A B^T)_{ij} + \mu_i)$
The matrices A and B are constrained by spike-and-slab lasso priors to induce sparsity. An optional IBP-style reordering helps emphasize interpretable structures during iteration.

---

## 📘 Documentation

You can view the documentation for all exported functions via:

```r
?BiSSLB
?BiSSLB_ladder
```

---

## 📄 License

This package is licensed under the GPL-3 License. See the `LICENSE` file for details.
