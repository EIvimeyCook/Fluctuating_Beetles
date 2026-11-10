<div align="center">
 <h1>Fluctuating_Beetles</h1>
</div>

<!-- badges: start -->
[![DOI](https://img.shields.io/badge/DOI-10.1093%2Fjeb%2Fvoad009-blue)](https://doi.org/10.1093/jeb/voad009)
[![Dryad](https://img.shields.io/badge/Dryad-10.5061%2Fdryad.f1vhhmgz7-blue)](https://doi.org/10.5061/dryad.f1vhhmgz7)
[![Zenodo](https://img.shields.io/badge/Zenodo-10.5281%2Fzenodo.10118423-blue)](https://doi.org/10.5281/zenodo.10118423)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](LICENSE.md)
<!-- badges: end -->
Code for the paper:

> Ivimey-Cook, E. R., Piani, C., Hung, W. T., & Berg, E. C. (2024). Genetic
> background and thermal regime influence adaptation to novel environment in the
> seed beetle, *Callosobruchus maculatus*. *Journal of Evolutionary Biology*,
> 37(1), 1–13. <https://doi.org/10.1093/jeb/voad009>

## Overview

How a population adapts to a novel environment depends both on the genetic
variation it starts with and on the environment it experiences while adapting.
Fluctuating thermal regimes are of particular interest under climate change,
because a population adapting to variable temperatures may follow a very
different trajectory from one adapting to a constant mean — even when the mean is
identical.

Using the seed beetle *Callosobruchus maculatus*, this study tests how **genetic
background** and **thermal regime** jointly shape adaptation to a novel
environment. Beetles from two genetic backgrounds evolved for 19 generations at
either a constant 29 °C or a fluctuating daily cycle (mean 33 °C, range 26–40 °C),
and were then assayed in their evolved environment or in the other environment.
The study measures consequences for body mass, development time, lifetime
reproductive success, and age-specific reproduction.

## Authors

Edward R. Ivimey-Cook, Claudio Piani, Wei-Tse Hung, and Elena C. Berg

## Affiliations

School of Biological Sciences, University of East Anglia, UK
The American University of Paris, Paris, France

## Data

Data have been uploaded to Dryad: <https://doi.org/10.5061/dryad.f1vhhmgz7>

Two files are read by the analysis, `bodymass.csv` and `lifehistory.csv`.

## Code

Code for loading libraries, loading data, modelling, and visualising.

**`01-load_packages.R`** — Packages for use in the analysis and graphing of model
outputs. Also defines `modelchecker(check, filename)`, a helper that runs `DHARMa`
zero-inflation, dispersion and uniformity tests on a list of candidate `glmmTMB`
models. It keeps models with zero-inflation and dispersion ratios between 0.80 and
1.10, a non-missing AIC and passing convergence checks. Uniformity results are
recorded but not used as a filter. The kept models are sorted by AIC. It opens
both tables with `View()`, so it needs an interactive session, and writes them to
`<filename>all.models.csv` and `<filename>final.models.csv`. `filename` sets the
prefix for these files, so give one in each `modelchecker()` call in
`03-model.R` (lines 246, 294, 337 and 494).

**`02-load_clean_data.R`** — Loading and cleaning data used in the analysis and
graphing.

**`03-model.R`** — Models for body mass, development time, lifetime reproductive
success, and age-specific reproduction. Includes code for all pairwise
comparisons.

**`04-figures.R`** — Plots estimated means from the fitted models and saves four
PNGs to `figures/`: `bodymass.png` and `development_time.png` (means plus sexual
dimorphism contrasts), `lrs.png` (LRS means) and `fitness.png` (LRS and
age-specific reproduction). Create the `figures/` and `outputs/` folders before
running; `03-model.R` writes model-selection CSVs to `outputs/`.

**Note:** `04-figures.R` styles the figures with `theme_publication()`, a custom
ggplot2 theme. Define it before running the script, for example
`theme_publication <- function(...) ggplot2::theme_bw(...)`; the published
figures used the custom theme.

## Running the analysis

The scripts are written for this project layout:

```
project/
  bodymass.csv        # from Dryad
  lifehistory.csv     # from Dryad
  code/               # move the four .R scripts here
  figures/            # create this empty folder
  outputs/            # create this empty folder
```

Set `project/` as the working directory. All paths are relative to it, and the
scripts source `code/01-load_packages.R` and `code/02-load_clean_data.R`.

Work through `03-model.R` and `04-figures.R` interactively in the same R session,
section by section, as `04-figures.R` plots the models fitted in `03-model.R`.
The body mass and development time models are both named `m1a`, so plot the
body mass model before fitting the development time model.

## R environment

`tidyverse`, `hablar`, `glmmTMB`, `DHARMa`, `emmeans`, `ggstatsplot`, `easystats`,
`MASS`, `ggeffects`, `patchwork`, `magrittr`, `MuMIn`, `broom`, `sjPlot`, `renv`,
`janitor`, `dabestr`, `brms`.

## Related work

- [**Heatwave_Beetles**](https://github.com/EIvimeyCook/Heatwave_Beetles) —
  long-term evolution under heatwave conditions in the same study system

## Citation

> Ivimey-Cook, E. R., Piani, C., Hung, W. T., & Berg, E. C. (2024). Genetic
> background and thermal regime influence adaptation to novel environment in the
> seed beetle, *Callosobruchus maculatus*. *Journal of Evolutionary Biology*,
> 37(1), 1–13. <https://doi.org/10.1093/jeb/voad009>

A machine-readable [`CITATION.cff`](CITATION.cff) is included.

## Contact

Edward R. Ivimey-Cook — <e.ivimeycook@gmail.com> —
[ORCID 0000-0003-4910-0443](https://orcid.org/0000-0003-4910-0443)

## License
Code is licensed under the [MIT License](LICENSE.md).
