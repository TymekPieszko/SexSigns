# Detecting and quantifying rare sex in natural populations (Pieszko et al., 2025)

Repository accompanying the manuscript (add DOI): an ARG-based analytical framework to explore reproductive mode signaures. The directories contain:

- `sim_pipeline/`
  - Snakemake pipeline used to run simulations using models of asexuality (`scripts/model_GC.slim` and `scripts/model_CO.slim`).
- `sexsign_utils/`
  - Functions for calculating classical statistics, genealogical statistics and summarising tree composition. Note that `calc_tcomp_vTracked` and `calc_tcomp_vRanks` (in `calc.py`) calculate tree composition in an equivalent manner.
- `stats/`
  - Scripts to calculate statistics and summarise tree composition across the space of simulated scenarios. Note that the parameter `bias` in the codebase is defined as 1 - *β*, where *β* denotes the bias parameter described in the manuscript.
- `plots/`
  - Scripts to generate all plots.
- `validation/`
  - Validation of simulation models with analytical results; see SI Appendix, section A.
- `abc/`
  - Analysis of the power of different metrics to predict the rate of sex using approximate Bayesian computation.
- `yeast/`
  - Application of the new anlytical framework to the haplotype-resolved genomic dataset for *Saccharomyces cerevisiae* (Loegler et al. (2025); [DOI](https://doi.org/10.1038/s41586-025-09637-0)).   