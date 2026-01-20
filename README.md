# Sampled Red List Index (SRLI).

The README file and the Codes are generated using Claude Opus 4.5 and reviewed by Kwanmok Kim.

## Overview
R implementation of the Sampled Red List Index methodology for assessing biodiversity trends, based on Baillie et al. (2008).

This repository contains R scripts to calculate the Red List Index (RLI) and determine minimum adequate sample sizes for monitoring extinction risk trends, following the methodology described in:

> Baillie, J.E.M., Collen, B., Amin, R., Akcakaya, H.R., Butchart, S.H.M., Brummitt, N., Meagher, T.R., Ram, M., Hilton-Taylor, C., & Mace, G.M. (2008). Toward monitoring global biodiversity. *Conservation Letters*, 1(1), 18-26. https://doi.org/10.1111/j.1755-263X.2008.00009.x

This code was used for the project CBD Headline indicator 4.2 Genetic Diversity Indicator.

This was necessary to calculate the minimum sample size (number of species) that could represent the genetic diversity of South Korea.  

The SBSTTA report suggests following this approach (Baillie et al. 2008), so the following codes were recreated (using Claude OPUS 4.5). Jan. 2026.

## Background

The Red List Index (RLI) measures trends in extinction risk based on IUCN Red List assessments. It enables tracking of overall extinction risk for taxonomic groups over time. The Sampled Red List Index (SRLI) extends this approach by using stratified sampling, making it feasible to assess large, species-rich taxonomic groups without requiring complete censuses.

Codes are meant to be compatible with Korean fonts hence requires additional package installation. 

Therefore, you can skip some of the functions and steps if you have a clean data. 


### Key Concepts

- **RLI** ranges from 0 (all species extinct) to 1 (all species Least Concern)
- **ΔRLI** measures the change in extinction risk between two time points
- **5% Sign Rule**: Sample size is adequate when <5% of bootstrap replicates show the wrong trend direction

## RLI Weights

Following the Baillie et al. (2008) methodology, IUCN categories are assigned the following weights:

| Category | Weight |
|----------|--------|
| Least Concern (LC) | 0 |
| Near Threatened (NT) | 1 |
| Vulnerable (VU) | 2 |
| Endangered (EN) | 3 |
| Critically Endangered (CR) | 4 |
| Regionally Extinct / Extinct in Wild / Extinct (RE/EW/EX) | 5 |

**Note:** Data Deficient (DD) species are excluded from RLI calculations as they carry no information about threat status.

## Repository Contents

```
SRLI_BaillieEtal2008/
├── README.md
├── .gitignore
└── sampleSize_Baillie2008_trendAnalysis.R # Trend analysis and visualization
```

## Key Functions

### `calc_rli()`
Calculates the Red List Index for a given dataset.

### `bootstrap_srli()`
Performs bootstrap analysis to determine confidence intervals and assess sample size adequacy.

### `recommend_n_sign()`
Recommends minimum sample size based on the Baillie 5% sign rule.

### `plot_srli_single_panel()` / `plot_srli_two_panel()`
Visualizes SRLI results with confidence intervals.

## Usage

```r
# Source the functions
source("sampleSize_Baillie2008.R")
source("sampleSize_Baillie2008_trendAnalysis.R")

# Calculate RLI
rli_t1 <- calc_rli(data, time = "t1")
rli_t2 <- calc_rli(data, time = "t2")

# Bootstrap analysis
bt <- bootstrap_srli(t1, t2, R = 50000, mode = "overlap", 
                     n_min = 20, step_n = 10, seed = 10)

# Plot results
plot_srli_single_panel(bt, wrong_thresh = 0.05)

# Get recommended sample size
n_rec <- recommend_n_sign(bt, wrong_thresh = 0.05)
```

## Data Requirements

Your data should include:
- Species identifiers
- IUCN threat categories for each time point (t1, t2)
- Categories can be in Korean format (e.g., "최소 관심(LC)") — the code extracts English codes automatically

## Methodology Notes

### Sample Size Determination

Following Baillie et al. (2008), this implementation:
1. Performs 50,000 bootstrap replicates at each sample size
2. Calculates the proportion showing the wrong trend direction
3. Identifies the minimum sample size where error rate falls below 5%
4. Recommends a minimum of 900 assessed species (or 1,500 total to allow for up to 40% Data Deficient)

### Overlap Mode for ΔRLI

When calculating change between time periods, "overlap" mode ensures that only species assessed at both time points are included, providing a true measure of status change rather than sampling artifacts.

## Application to South Korean Biodiversity

This code was developed to analyze the Korean National Red List (국가생물적색목록) to:
- Assess trends in extinction risk for Korean species
- Determine adequate sample sizes for future monitoring
- Support biodiversity indicator development aligned with international standards

## License

MIT License - see LICENSE file for details.

## Citation

If you use this code, please cite both this repository and the original methodology:

```
Baillie, J.E.M., Collen, B., Amin, R., Akcakaya, H.R., Butchart, S.H.M., 
Brummitt, N., Meagher, T.R., Ram, M., Hilton-Taylor, C., & Mace, G.M. (2008). 
Toward monitoring global biodiversity. Conservation Letters, 1(1), 18-26.
https://doi.org/10.1111/j.1755-263X.2008.00009.x
```

## References

- Baillie, J.E.M., et al. (2008). Toward monitoring global biodiversity. *Conservation Letters*, 1(1), 18-26.
- Butchart, S.H.M., et al. (2004). Measuring global trends in the status of biodiversity: Red List Indices for birds. *PLoS Biology*, 2, 2294-2304.
- Butchart, S.H.M., et al. (2007). Improvements to the Red List Index. *PLoS ONE*, 2, e140.
- IUCN (2001). IUCN Red List Categories and Criteria: Version 3.1. IUCN, Gland, Switzerland.

## Contact

Kwanmok Kim
kwanmok.kim@gmail.com
Korea Environment Institute (KEI)
