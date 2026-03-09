# pmvalsampsize 

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/pmvalsampsize)](https://CRAN.R-project.org/package=pmvalsampsize)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/pmvalsampsize)](https://CRAN.R-project.org/package=pmvalsampsize)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/pmvalsampsize)](https://CRAN.R-project.org/package=pmvalsampsize)
[![License: GPL (>= 3)](https://img.shields.io/badge/License-GPL%20(%3E%3D%203)-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

## Overview

**pmvalsampsize** computes the minimum sample size required for the external validation of an existing multivariable prediction model. It implements the criteria proposed by:

- **Archer et al. (2020)** for models with a **continuous outcome**
- **Riley et al. (2021)** for models with a **binary outcome**

A series of criteria define the sample size needed to ensure precise estimation of key measures of prediction model performance, allowing conclusions to be drawn about whether the model is potentially accurate and useful in a given target population.

> **Note:** Support for survival (time-to-event) outcome models is currently in development.

## Installation

Install the released version from CRAN:

```r
install.packages("pmvalsampsize")
```

Or install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("JoieEnsor/pmvalsampsize")
```

## Usage

The package exports a single main function, `pmvalsampsize()`, which calculates the minimum sample size for external validation of a prediction model.

### Binary outcome example

Calculate the minimum sample size to externally validate a logistic prediction model for a binary outcome (e.g. mechanical heart valve failure). Based on previous evidence, the outcome prevalence is anticipated to be 1.8% and the reported C-statistic was 0.8.

```r
library(pmvalsampsize)

pmvalsampsize(
  type       = "b",
  prevalence = 0.018,
  cstatistic = 0.8,
  lpcstat    = -4.7,
  oeciwidth  = 1,
  seed       = 1234,
  simobs     = 10000
)
```

### Key parameters

| Parameter | Description |
|-----------|-------------|
| `type` | Outcome type: `"b"` for binary, `"c"` for continuous |
| `prevalence` | Anticipated prevalence of the outcome (binary models only) |
| `cstatistic` | Reported C-statistic of the existing model |
| `cstatciwidth` | Target CI width for the C-statistic (default: 0.1) |
| `cslope` | Expected calibration slope (default: 1) |
| `csciwidth` | Target CI width for the calibration slope (default: 0.2) |
| `oe` | Expected O:E ratio (default: 1) |
| `oeciwidth` | Target CI width for the O:E ratio (default: 0.2) |
| `simobs` | Number of observations for simulation (default: 1,000,000) |
| `graph` | Logical; if `TRUE`, produces a plot of the sample size calculations |
| `seed` | Random seed for reproducibility |

For the full list of parameters, see `?pmvalsampsize`.

## Background

External validation is a critical step in the lifecycle of a clinical prediction model. It assesses how well a model developed in one population performs in a new, independent dataset. Selecting an adequate sample size for a validation study is essential to ensure reliable conclusions about model performance.

`pmvalsampsize` formalises this process by providing minimum sample sizes based on the following criteria:

1. **Criterion 1** – Precise estimation of the calibration slope
2. **Criterion 2** – Precise estimation of the O:E ratio (binary/survival) or mean squared error (continuous)
3. **Criterion 3** – Precise estimation of the C-statistic (binary/survival) or R² (continuous)
4. **Criterion 4** *(binary only)* – Precise estimation of net benefit at a specified risk threshold

The minimum sample size recommended is the **maximum** across all applicable criteria.

## References

Archer L, Snell K, Ensor J, Hudda MT, Collins GS, Riley RD. Minimum sample size for external validation of a clinical prediction model with a continuous outcome. *Statistics in Medicine.* 2020. <https://doi.org/10.1002/sim.8766>

Riley RD, Debray TPA, Collins G, Archer L, Ensor J, van Smeden M, Snell KIE. Minimum sample size for external validation of a clinical prediction model with a binary outcome. *Statistics in Medicine.* 2021. <https://doi.org/10.1002/sim.9025>

## Citation

If you use `pmvalsampsize` in your work, please cite:

```bibtex
@Manual{,
  title  = {pmvalsampsize: Sample Size for External Validation of a Prediction Model},
  author = {Joie Ensor},
  year   = {2023},
  note   = {R package version 0.1.0},
  url    = {https://CRAN.R-project.org/package=pmvalsampsize},
}
```

## Related packages

- [**pmsampsize**](https://CRAN.R-project.org/package=pmsampsize) – Sample size for *development* of a prediction model

## License

GPL (>= 3) — see [LICENSE](LICENSE) for details.

## Author

**Joie Ensor** — University of Birmingham ([j.ensor@bham.ac.uk](mailto:j.ensor@bham.ac.uk))

Bugs and feature requests can be submitted via the [GitHub Issues](https://github.com/JoieEnsor/pmvalsampsize/issues) page.
