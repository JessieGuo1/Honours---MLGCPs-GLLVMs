# Fast Computation Methods for MLGCPS and GLLVMs

Statistics Honours research repository comparing computational methods for fitting **Multivariate Log-Gaussian Cox Processes (MLGCPs)** and **Generalised Linear Latent Variable Models (GLLVMs)**.

This project evaluates the accuracy and computational cost of standard and alternative estimation procedures using simulated count data and ecological point-process data.

## Research overview

The repository contains R and R Markdown workflows for:

- fitting Poisson GLLVMs with one and two latent variables;
- comparing Laplace Approximation and Variational Approximation methods;
- implementing Alternating GLM Estimation;
- implementing penalised alternatives with `glmnet`;
- comparing Joint Maximum-Likelihood Estimation variants;
- fitting Univariate and Multivariate Log-Gaussian Cox Processes;
- assessing fitted models using Conditional Log-Likelihood and KL Divergence;
- comparing computation times, latent-variable representations, and fitted intensity maps.

The code uses the `gllvm` package for standard GLLVM fitting, and the saved `.RDATA` files contain fitted models, simulated datasets, summary statistics, and plotting objects.

## Repository structure

- `GLLVM Fitting.Rmd` — fits GLLVMs, Alternating GLM models, Joint Maximum-Likelihood models, and penalised alternatives to ecological count data.
- `GLLVM Simulations.Rmd` — generates simulated datasets and compares fitting methods over repeated runs.
- `mlgcp - results.Rmd` — creates thesis figures for intensity maps, timing, Log-Likelihood, and biplots.
- `lgcp fitting - gam + scampr.Rmd` — fits Univariate Log-Gaussian Cox Processes using GAM and `scampr` approaches.
- `lgcp fitting - glm + glmnet.Rmd` — fits Univariate Log-Gaussian Cox Processes using GLM and penalised GLM approaches.
- `alt_glm.R`, `alt_glmnet.R`, `mlgcp_alt_glm.R`, and `mlgcp_alt_glmnet.R` — custom Alternating Estimation implementations.
- `jmle.R` — Joint Maximum-Likelihood Estimation implementation.
- `build_basis_functions.R`, `proc.R`, `get_and_format_flora_data.R`, and `obtain_format_data.R` — support functions for basis construction and data processing.
- `gllvm_models/`, `glm_models/`, `glmnet_models/`, `jmle_models/` — directories for saved model fits and simulation summaries.
- `figures/`, `results/`, and `presentation/` — generated plots, result objects, and presentation materials.
- `*.RDATA` — saved datasets and fitted objects used throughout the analyses.

## Methods compared

### Generalised Linear Latent Variable Models

- `gllvm` with Laplace Approximation (LA)
- `gllvm` with Variational Approximation (VA)
- Alternating GLM fitting
- Alternating `glmnet` fitting
- Joint Maximum-Likelihood Estimation
- Penalised Estimation with multiple penalty strengths

### Spatial and Multivariate Models

- Univariate Log-Gaussian Cox Processes
- Multivariate Log-Gaussian Cox Processes
- Basis-function approaches
- Joint and individual latent-variable update schemes - 'basis function updates'
- Joint and one-quadrature formulations

## Requirements

An R installation and RStudio or another R Markdown-compatible environment are recommended. The workflows use packages including:

- `gllvm`
- `glmnet`
- `MASS`
- `dplyr`
- `ggplot2`
- `tmap`
- `spatstat`
- `mgcv`
- `scampr`
- `stringr`

A `glmmTMB` source archive is also included for reproducibility.

Install missing packages from CRAN where applicable, for example:

```r
install.packages(c("gllvm", "glmnet", "MASS", "dplyr", "ggplot2",
                   "tmap", "spatstat", "mgcv", "stringr"))
```

Package availability and model dependencies may vary by operating system and R version.

## Getting started

Clone the repository and open it in RStudio:

```bash
git clone https://github.com/JessieGuo1/Honours---MLGCPs-GLLVMs.git
cd Honours---MLGCPs-GLLVMs
```

Open an R Markdown file, set the repository root as the working directory, and run it with **Knit** or by executing chunks in sequence. Many files load saved `.RDATA` objects before producing plots or summaries, so preserve the repository directory structure.

Suggested starting points:

1. `GLLVM Fitting.Rmd` for fitting and comparison of GLLVM methods.
2. `GLLVM Simulations.Rmd` for Simulation-based Evaluation.
3. `lgcp fitting - gam + scampr.Rmd` or `lgcp fitting - glm + glmnet.Rmd` for LGCP fitting.
4. `mlgcp - results.Rmd` for result and figure generation.

## Reproducibility

The simulation scripts set seeds and store generated datasets, fitted model objects, timing measurements, and summary quantities. Relative paths are used extensively, so scripts should be run from the repository root.

Model fitting can be computationally intensive. Several workflows save intermediate fitted objects and permit previously generated `.RDATA` files to be loaded rather than rerunning every model.

## Author

Jessie Guo, David Warton

## Licence

No licence is currently specified. Please contact the repository owner before reusing the code or results in another project.
