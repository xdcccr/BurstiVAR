# Supplement D: Two-burst simulation check

This folder reproduces the B = 2 columns of appendix Table F1 (`tab:supplement_2burst`) and the convergence results reported in Supplement D of the BurstiVAR paper.

The B = 3 comparison columns of Table F1 are not produced here. They come from the Study 1 cell with T_b = 3; see the `study1_burstivar` folder.

## Design

B = 2 measurement bursts, T_b = 3 time points per burst, T = 6 total time points, with N in {100, 500} persons and 100 Monte Carlo replications per condition. The generating model is the Study 1 model restricted to its first two bursts: burst-specific intercept means (0, 1) for burst 1 and (1, 1.5) for burst 2, intercept SDs (1, 1) and (1.5, 2), AR means 0.3 and 0.2, cross-regression means -0.15 (y1 to y2) and 0 (y2 to y1), random-effect SD 0.1 for all four VAR parameters, innovation variances 1.0 with innovation covariance 0.3. The true values are listed in `TrueParameters_mlGVARNoCorNoME_BurstIntercept_2Burst_3T.csv`.

## Files

- `mlGVARNoCorNoME_BurstIntercept_2Burst_3T_DataGeneration.R` — generates the simulated datasets.
- `mlGVARNoCorNoME_BurstIntercept_2Burst.txt` — JAGS model for the two-burst BurstiVAR.
- `mlGVARNoCorNoME_BurstIntercept_2Burst_3T_Modelfit.R` — fits the JAGS model to every replication in parallel.
- `posteriorSummaryStats.R` — helper sourced by the model-fit script; summarizes posteriors (mean, PSD, quantile and HDI intervals, ESS, R-hat).
- `Generate_Summ_BurstIntercept_2Burst_3T.R` — aggregates the per-replication result tables into the Monte Carlo summary files.
- `TrueParameters_mlGVARNoCorNoME_BurstIntercept_2Burst_3T.csv` — true parameter values used for evaluation.
- `results/` — the archived output files behind the paper:
  - `..._MCfile_nT6_nP100.csv`, `..._MCfile_nT6_nP500.csv` — per-replication true value, posterior mean, posterior SD, 95% CI limits, and coverage/significance flags for each parameter.
  - `..._MCfileSumm_nT6_nP100.csv`, `..._MCfileSumm_nT6_nP500.csv` — aggregated recovery statistics per parameter (mean estimate, RMSE, relative bias, SD of estimates, mean posterior SD, RDSE, coverage, power or Type I error, missing proportion). These two files feed the B = 2 columns of Table F1.

## Pipeline

Run the scripts from this directory, in this order:

1. `mlGVARNoCorNoME_BurstIntercept_2Burst_3T_DataGeneration.R` — writes 100 simulated datasets per N to `./data/` and the true-parameter CSV.
2. `mlGVARNoCorNoME_BurstIntercept_2Burst_3T_Modelfit.R` — fits the JAGS model to each dataset and writes a posterior summary table plus the coda samples (`.RData`) per replication to `./result/`. The script skips replications whose result table already exists, so it can be rerun to fill in failures.
3. `Generate_Summ_BurstIntercept_2Burst_3T.R` — reads the result tables in `./result/` and writes the `MCfile` and `MCfileSumm` CSVs there. These should match the copies archived in `results/`.

Convergence was assessed from the RHAT and ESS columns of the per-replication posterior summary tables produced in step 2; the convergence prose in Supplement D summarizes those columns.

## Compute environment

R 4.3.2, JAGS 4.3.1, rjags. MCMC settings per replication: 2 chains, 5000 adaptation iterations, 5000 burn-in iterations, 20000 sampling iterations, with chain seeds (r, r + 500) for replication r.

## Runtime and caveats

- The published B = 2 run was executed fully on a local machine with `num_cores = 2` (2026-02-03), unlike the larger Study 1 grid. With 2 cores, expect the 200 fits (2 conditions x 100 replications) to take several hours in total; raise `num_cores` in the model-fit script to shorten this.
- The model-fit script saves the full coda samples for every replication, which takes substantial disk space. Delete the `.RData` files after aggregation if space is a concern.
- Data generation is seeded (`set.seed(r)` per replication), so regenerated datasets are reproducible. MCMC results may still differ in late decimal places across JAGS builds and platforms.
- Scripts carry a provenance header stating the original archive path. Apart from that header and the path-configuration line (absolute paths replaced so the scripts run from this directory), the scripts are identical to the versions that produced the published results.
- Some comments and console messages in the scripts are in Chinese; they do not affect execution.
