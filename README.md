# bcirr

`bcirr` is an R package for estimating between-case incidence rate
ratios (IRRs) using generalized estimating equations (GEEs). It is
designed specifically for across-participants multiple-baseline
single-case experimental designs (SCEDs) with small sample sizes and
implements small-sample corrections, including the Kauermann-Carroll
correction. The package includes the `OtaDuPaul2002` dataset for
demonstration.

## Installation

You can install the development version of `bcirr` directly from GitHub:

    # Install devtools if you have not already
    if (!requireNamespace("devtools", quietly = TRUE)) {
      install.packages("devtools")
    }

    # Install bcirr
    devtools::install_github("WenLuoTAMU/BC-IRR")

## Usage Example

The package includes the dataset `OtaDuPaul2002` for demonstration. The
following script loads the data, formats it, and runs the bias-corrected
GEE model.

    library(bcirr)
    library(tidyverse)

    # 1. Load the dataset
    data("OtaDuPaul2002")

    # 2. Prepare the data
    # We select the relevant columns and rename them for the model
    data <- OtaDuPaul2002 %>%
      select(Subject, Session, Y, Phase) %>%
      rename(
        id = Subject,
        time = Session,
        outcome = Y,
        binary = Phase
      )

    # 3. Fit the GEE model with bias correction
    # Family is poisson because the outcome is count/frequency data
    fit_auto <- bc_irr_gee(
      formula = outcome ~ binary * time,
      id = "id",
      data = data,
      family = poisson,
      se_method = "kc" # Kauermann-Carroll correction
    )

    # 4. View the results
    summary(fit_auto)

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
