<div align="center">
 <h1>Fluctuating_Beetles</h1>
</div>

<!-- badges: start -->
[![Paper](https://img.shields.io/badge/paper-10.1093%2Fjeb%2Fvoad009-blue)](https://doi.org/10.1093/jeb/voad009)
[![Dryad](https://img.shields.io/badge/Dryad-10.5061%2Fdryad.f1vhhmgz7-blue)](https://doi.org/10.5061/dryad.f1vhhmgz7)
[![Zenodo](https://img.shields.io/badge/Zenodo-10.5281%2Fzenodo.10118423-blue)](https://doi.org/10.5281/zenodo.10118423)
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
environment, measuring consequences for body mass, development time, lifetime
reproductive success, and age-specific reproduction.

## Authors

Edward R. Ivimey-Cook, Claudio Piani, Wen Ting Hung, and Elena C. Berg

## Affiliations

School of Biological Sciences, University of East Anglia, UK
The American University of Paris, Paris, France

## Data

Data have been uploaded to Dryad: <https://doi.org/10.5061/dryad.f1vhhmgz7>

Two files are read by the analysis, `bodymass.csv` and `lifehistory.csv`.

## Code

Code for loading libraries, loading data, modelling, and visualising.

**`01-load_packages.R`** — Packages for use in the analysis and graphing of model
outputs. Also defines `modelchecker()`, a helper that screens a list of candidate
`glmmTMB` models on `DHARMa` diagnostics (zero-inflation, dispersion, uniformity),
checks convergence, sorts by AIC, and writes both the full candidate set and the
filtered subset to CSV.

**`02-load_clean_data.R`** — Loading and cleaning data used in the analysis and
graphing.

**`03-model.R`** — Models for body mass, development time, lifetime reproductive
success, and age-specific reproduction. Includes code for all pairwise
comparisons.

**`04-figures.R`** — Creates figures for all the pairwise interactions and model
outputs, written to a `figures/` directory.

Scripts should be run in numerical order.

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
