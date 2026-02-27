# SeMPLE mixed-effects

SeMPLE mixed-effects is a simulation‑based inference framework for stochastic nonlinear mixed‑effects models. It extends the SeMPLE methodology (see [henhagg/semple](https://github.com/henhagg/semple)) by learning surrogate models for the likelihood function and the posterior distribution of parameters in hierarchical systems.

This repository contains the code accompanying version v2 of the the paper at https://arxiv.org/abs/2504.11279 (the current version, v1, will be updated soon).

PEPSDI computations can be reproduced via the `gfp.jl` script in the PEPSDI fork [henhagg/PEPSDI](https://github.com/henhagg/PEPSDI).

---

## Repository structure

- `algorithm/` – core R code implementing SeMPLE algorithms and helper functions. Key scripts include:
  - `semple_mem.R` – entry point for mixed‑effects SeMPLE computations.

- `models/` – model definitions and generators for simulated data.
  - Subfolders for each model variant (`ornstein_uhlenbeck*`, `mrna*`, etc.).
  - Each contains an R script with the model definition to simulate data, prior definitions and observed/generated data.
  - `ornstein_uhlenbeck_kalman/src/` – Julia code implementing an exact Kalman filter for the Ornstein–Uhlenbeck model. This provides a benchmark for likelihood computation and is used for comparison studies in the paper.

- `run_scripts/` – high‑level scripts that perform entire analyses, e.g.:
  - `run_ornstein_uhlenbeck*.R`. These call the functions in `algorithm/` for a given model and save output under `results/`. - Running these scripts with the current settings will reproduce all SeMPLE results in the paper.

- `plot_scripts/` – R scripts used to plot figures. Running `plot_results.R` will recreate all figures in the paper.

- `results/` – directory populated by the `run_scripts`. Each algorithm run is saved in a subfolder under `results/<model>/<observation>/<subfolder>`.

- `renv/` – the [renv](https://rstudio.github.io/renv/) environment used to create a reproducible environment. Associated with the `renv.lock` file.


## Getting started

### Prerequisites

1. **R** Ensure a compiler and `Rcpp`/`rstan` toolchains are installed.
2. Clone the repository and set the working directory to the project root.

```sh
git clone https://github.com/henhagg/semple_mem.git
cd semple_mem
```

3. Activate the R environment and restore packages:

```r
# install.packages("renv")
# start R from project root
renv::activate()
renv::restore()
```

This will install all required packages (e.g. `mvtnorm`, `ggplot2`, etc.).

### Running analyses

Each `run_scripts/*.R` file is self‑contained. The output is written to `results/<model>/<observation>/<subfolder>` by default.

Example:

```r
source("run_scripts/run_mrna_fix_tzero.R")
```

For BIC or runtime studies (Fig. 6 and 7 in the paper) run the corresponding scripts:

```r
source("run_scripts/run_bic_computations.R")
source("run_scripts/runtime_comparison.R")
```

### Generating plots

Execute the plotting scripts in `plot_scripts/` to recreate figures. Running

```r
source("plot_scripts/plot_results.R")
```
will recreate all figures in the paper.

### PEPSDI calculations

To reproduce the Julia‑based PEPSDI benchmarks, setup the Julia project [henhagg/PEPSDI](https://github.com/henhagg/PEPSDI) and run

```sh
cd /path/to/PEPSDI
julia gfp.jl
```
The PEPSDI results reported in the paper can be found under `results/PEPSDI`.

## Licensing

This project is released under the MIT License (see `LICENSE`).

---
