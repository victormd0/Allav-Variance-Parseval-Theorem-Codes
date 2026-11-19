# A Version of Parseval's Theorem for the Allan Variance

This repository contains GNU Octave/MATLAB code developed for the numerical investigation of Parseval-type identities involving the **Allan variance (AVAR)**.

The code accompanies the research paper:

**A Version of Parseval's Theorem for the Allan Variance**  
Victor Mosimann Duarte, François Vernotte, and Enrico Rubiola.

## Overview

The objective of this work is to investigate how the power of a stochastic signal can be recovered directly from its Allan variance, without requiring an explicit reconstruction of its power spectral density (PSD).

The repository focuses on the continuous-time and infinite-discrete theoretical framework developed in the accompanying paper, with numerical experiments performed on finite-length discrete signals.

In particular, the numerical experiments investigate two quantities:

- $a_0$: the power of the increment sequence $y[k]=x[k+1]-x[k]$.
- $b_0$: the power of the original sequence $x[k]$.

The corresponding infinite-discrete Parseval-type identities are

$$
a_0 = \frac{1}{2}\sum_{j=0}^{\infty}\mathrm{AVAR}_y[2^j]
$$

and

$$
b_0 = \lim_{M\to\infty}\frac{1}{3M}
\sum_{m=1}^{M}m^2\mathrm{AVAR}_y[m]
$$


The numerical experiments compare estimates obtained through AVAR inversion with estimates obtained through direct autocorrelation and with theoretical values.

## Repository contents

### Core functions

| File | Description |
|---|---|
| `power_law.m` | Generates synthetic signals with power-law spectral behavior. |
| `AVAR_values_y.m` | Computes Allan variance values from discrete frequency sequences. |
| `AVAR_values_x.m` | Computes Allan variance values from time-error sequences. |
| `a_values_y.m` | Computes coefficients associated with the increment sequence using AVAR-based inversion. |
| `b_values_x.m` | Computes coefficients associated with the original sequence using AVAR-based inversion. |

### Numerical experiments

| File | Description |
|---|---|
| `Test_power_law.m` | Tests the power-law signal generator. |
| `Test_AVAR_y_gen.m` | Tests Allan variance calculations for generated frequency signals. |
| `Test_AVAR_x_gen.m` | Tests Allan variance calculations for generated time-error signals. |
| `Test_a_values_y_gen.m` | Tests AVAR-based estimation of the coefficients $a_n$. |
| `Test_b_values_x_gen.m` | Tests AVAR-based estimation of the coefficients $b_n$. |
| `Test_a_histogram_y_gen.m` | Compares distributions of AVAR-based and autocorrelation-based estimates of $a_0$. |
| `Test_b_histogram_x_gen.m` | Compares distributions of AVAR-based and autocorrelation-based estimates of $b_0$. |
| `Test_a_values_y_vs_N.m` | Investigates the dependence of the $a_n$ estimates on record length. |
| `Test_b_values_x_vs_N.m` | Investigates the dependence of the $b_n$ estimates on record length. |

## Numerical validation

The experiments accompanying the paper consider three power-law spectral models:

| Spectral exponent | Noise type | PSD behavior |
|---|---|---|
| $\alpha=0$ | White noise | $S(f)\propto f^0$ |
| $\alpha=1$ | Blue noise | $S(f)\propto f^1$ |
| $\alpha=2$ | Violet noise | $S(f)\propto f^2$ |

### Distribution of the estimates

For each spectral model, $K=100$ independent realizations of length $N=2^{10}=1024$ are generated.

The resulting estimates are compared through:

- Histograms of the AVAR-based estimates.
- Histograms of the direct autocorrelation estimates.
- Superimposed distributions.
- Comparisons between the two estimates for each realization.
- Sample means, standard deviations, and Pearson correlation coefficients.

### Dependence on record length

The convergence experiments use $K=100$ independent realizations for each record length.

The investigated ranges are:

- **$a_0$:** $N=2^8,2^9,\ldots,2^{20}$.
- **$b_0$:** $N=2^8,2^9,\ldots,2^{15}$.

The smaller maximum record length for $b_0$ is due to the computational cost of the corresponding numerical estimator, rather than a theoretical restriction.

## Requirements

The code is written for **GNU Octave**, using MATLAB-compatible syntax where possible.

The numerical functions rely on standard matrix and vector operations.

Some test scripts include plotting routines and, depending on the version, Windows-specific functionality for displaying simulation progress.

MATLAB compatibility has not been systematically verified.

## Usage

Download or clone the repository:

```bash
git clone https://github.com/victormd0/Allav-Variance-Parseval-Theorem-Codes.git
```

Open GNU Octave and navigate to the repository directory.

The core functions can then be called directly, while the `Test_*.m` files provide numerical experiments and validation routines.

For example, the following commands illustrate the basic workflow:

```matlab
K = 100;
N = 1024;
alpha = 0;
h = 1;

x = power_law(K, N, alpha, h);

avar = AVAR_values_x(x, 1:floor((N-1)/2));

b = b_values_x(x, 0);
```

The exact arguments and computational requirements depend on the function being used. Larger record lengths, particularly in the experiments involving $b_0$, may require substantial execution time.

## Reproducibility

The repository provides the computational routines used to investigate the Parseval-type identities and to compare AVAR inversion with direct autocorrelation.

Because the experiments use randomly generated signals, individual numerical results may vary between executions. Agreement should be assessed statistically rather than by expecting identical numerical outputs.

The theoretical results and their assumptions are discussed in the accompanying paper.

## Citation

If you use this code in academic work, please cite the accompanying paper:

Victor Mosimann Duarte, François Vernotte, and Enrico Rubiola, *A Version of Parseval's Theorem for the Allan Variance*.

Full publication information will be added when available.

## Authors

- Victor Mosimann Duarte
- François Vernotte
- Enrico Rubiola

## Acknowledgment of AI assistance

ChatGPT Plus (OpenAI) was used to assist with code development, documentation, and language editing. The mathematical ideas and theoretical derivations presented in the accompanying paper were developed by the authors, who remain responsible for the scientific content and numerical results.
