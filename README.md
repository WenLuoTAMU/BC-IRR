
# bcirr

`bcirr` is an R package for estimating between-case incidence rate
ratios (IRRs) using generalized estimating equations (GEEs). It is
designed specifically for across-participants multiple-baseline
single-case experimental designs (SCEDs) with small sample sizes and
implements small-sample corrections, including the Kauermann-Carroll
correction. The package includes the `OtaDuPaul2002` and `Cihak2010`
datasets to demonstrate count and percentage outcomes, respectively.

## Installation

You can install the development version of `bcirr` directly from GitHub:

``` r
# Install devtools if you have not already
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install bcirr
devtools::install_github("WenLuoTAMU/BC-IRR")
```

## Usage Example 1: Count Outcome

The package includes the dataset `OtaDuPaul2002` to demonstrate analysis
of count/frequency outcomes. The following script loads the data,
formats it, and runs the GEE model.

``` r
library(bcirr)
library(tidyverse)

# 1. Load the dataset
data("OtaDuPaul2002")

# 2. Prepare the data
# Select the relevant columns and rename them for the model
data <- OtaDuPaul2002 %>%
  select(Subject, Session, Y, Phase) %>%
  rename(
    id = Subject,
    time = Session,
    outcome = Y,
    binary = Phase
  )

# 3. Fit the GEE model
# Family is poisson because the outcome is count/frequency data
fit_auto <- bc_irr_gee(
  formula = outcome ~ binary * time,
  id = "id",
  data = data,
  family = poisson,
  outcome_type = "count", # Default
  se_method = "kc" # Kauermann-Carroll correction
)

# 4. View the results
summary(fit_auto)
```

## Usage Example 2: Percentage Outcome

The package includes the dataset `Cihak2010` to demonstrate analysis of
percentage outcomes. The following script loads the data, formats it,
and runs the GEE model.

``` r
library(bcirr)
library(tidyverse)

# 1. Load the dataset
data("Cihak2010")

# 2. Prepare the data
# Select the relevant columns and rename them for the model
data_percentage <- Cihak2010 %>%
  select(Subject, Session, Y, Phase, Total) %>%
  rename(
    id = Subject,
    time = Session,
    outcome = Y,
    binary = Phase,
    total = Total
  )

# 3. Fit the GEE model
# Percentage outcomes require the total number of intervals
fit_percentage <- bc_irr_gee(
  formula = outcome ~ binary * time,
  id = "id",
  data = data_percentage,
  family = poisson,
  outcome_type = "percentage", # Outcome recorded on a 0–100 scale
  total_var = "total",
  se_method = "md" # Mancl & DeRouen correction
)

# 4. View the results
summary(fit_percentage)
```

For proportion outcomes (coded from 0 to 1), use
`outcome_type = "proportion"` instead.

## Citation

If you use `bcirr` in your work, please cite the following preprints:

Li, H., Li, C., Luo, W., & Baek, E. (2025, October 17). *Between-case
incidence rate ratio for single case experimental designs with count
outcomes: A Monte Carlo simulation with conditional and marginal models*
\[Preprint\]. <https://doi.org/10.31234/osf.io/rx96g_v1>

Luo, W., Li, H., Baek, E., & Li, C. (2025, October 24). *Between-case
incidence rate ratio: A design comparable effect size for count outcomes
in single case experimental designs* \[Preprint\].
<https://doi.org/10.31219/osf.io/9gdkn_v2>

## Authors

- Chendong Li, James Madison University (JMU)
- Wen Luo, Texas A&M University (TAMU)
- Eunkyeng Baek, Texas A&M University (TAMU)
- Haoran Li, University of Minnesota (UMN)
- Lu Chen, Texas A&M University (TAMU)

## AI Usage Disclosure

Generative AI tools, including OpenAI’s ChatGPT (versions 5.4-5.6), were
used in an assistive capacity during the development of this R package
and its accompanying documentation. These tools supported selected
tasks, including reviewing code structure, suggesting potential
approaches to debugging and refactoring, assisting with documentation,
and improving the clarity and readability of written content.

All AI-assisted outputs were thoroughly reviewed, revised, and validated
by the authors. The authors independently made all substantive
methodological, statistical, architectural, and software-design
decisions. All code and documentation were examined and tested by the
authors before being incorporated into the package.

The authors take full responsibility for the accuracy, originality,
functionality, reproducibility, licensing compliance, and ethical and
legal integrity of all materials included in the package.
