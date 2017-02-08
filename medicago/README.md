In this directory:

- `data/` - Where the raw data goes.
- `medicago_data_setup.html` - some initial summaries of the data
- `medicago_data_setup.Rmd` - the code that produced those summaries
- `run_on_medicago.R` - computes local PC coordinates and MDS coordinates.
- `reports/` - Summary reports showing the results of running the algorithm with different parameter combinations
    (e.g., window choice; number of PCs) on the data.

- `summarize_run.Rmd` - templated report used to produce these.
- `compare-parameter-runs.R` - postprocessing to compare different parameter combinations.

First-round analyses (without using the R package):

- plots
- Medicago_distance_all_chr.R
- Medicago_MDS.R
- Medicago_PCA_win104.R
- Medicago_recode_and_cov.R

