# Detecting and quantifying rare sex in natural populations (Pieszko et al., 2026)

Repository accompanying the manuscript (currently on bioRxiv: [DOI](https://doi.org/10.1101/2025.06.03.657731)). The directories contain:

- `sim_pipeline/`
  - Snakemake pipeline used to run simulations under models of asexuality (`scripts/model_GC.slim` and `scripts/model_CO.slim`).
- `sexsigns_utils/`
  - Functions for calculating classical and genealogical statistics and for summarising tree composition. Note that `calc_tcomp_vTracked` and `calc_tcomp_vRanks` (in `calc.py`) calculate tree composition in an equivalent manner.
- `stats/`
  - Scripts used to calculate statistics and summarise tree composition across the space of simulated scenarios. Note that the parameter `bias` in the codebase is defined as 1 - *β*, where *β* denotes the bias parameter described in the manuscript.
- `plots/`
  - Scripts used to generate all plots.
- `validation/`
  - Validation of simulation models against analytical results; see SI Appendix, section A.
- `abc/`
  - Evaluation of different metrics for predicting the rate of sex using approximate Bayesian computation.
- `yeast/`
  - Application of the analytical framework to a haplotype-resolved genomic dataset for *Saccharomyces cerevisiae* (Loegler et al., 2025; [DOI](https://doi.org/10.1038/s41586-025-09637-0)).
- `rotaria/`
  - Application of the analytical framework to two *Rotaria* bdelloid rotifer species: variant calling (`calls/`), phased blocks, co-phased segments and SNP polarisation (`data/`), ARG inference and genealogical statistics (`arg_based/`), classical statistics (`classical/`) and ABC (`abc/`). Simulations for the ABC were run with `sim_pipeline/` (mutation rates tuned with `sim_pipeline/tune_mut_to_rotaria.py`) and summarised with `stats/calc_stats_for_abc.py`.
