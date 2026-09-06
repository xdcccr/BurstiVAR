# Study 1: BurstiVAR parameter recovery

This folder reproduces the Study 1 simulation of the BurstiVAR paper and the
four tables built from it:

- Table 2 (`tab:study1_trend`): burst-related (trend) parameters — burst-intercept
  means and between-person SDs.
- Table 3 (`tab:study1_results`): dynamic process parameters — fixed effects,
  random-effect SDs, innovation covariance.
- Appendix Table B1 (`tab:study1_trend_supp`): coverage, power/Type-I error,
  empirical SD, and RDSE for the trend parameters.
- Appendix Table B2 (`tab:study1_supp_metrics`): the same supplementary metrics
  for the dynamic process parameters.

## Design

Bivariate multilevel VAR with B = 3 measurement bursts and freely estimated
burst intercepts. Conditions cross burst length Tb in {3, 5, 20} (total
T = 9, 15, 60) with sample size N in {100, 500}. The cross-regression
mu_phi12 (y2 -> y1) is 0 in the data-generating model. 100 Monte Carlo
replications per cell.

## Contents

- `mlGVARNoCorNoME_BurstIntercept_3Burst.txt` — the JAGS model fit in every
  condition. The Tb = 5 and Tb = 20 fits referenced this same model under the
  filename `mlGVARNoCorNoME_BurstIntercept.txt`; the two files are
  byte-identical (only the name differs).
- `posteriorSummaryStats.R` — helper sourced by the model-fit scripts
  (posterior summary table from the coda samples).
- `tb3/` — pipeline for the Tb = 3 cells (data generation, model fit,
  summary aggregation).
- `tb5_tb20/` — pipeline for the Tb = 5 and Tb = 20 cells.
- `results/` — the per-condition Monte Carlo summary CSVs the table scripts
  read. Currently holds the four Tb = 5 / Tb = 20 summaries, verbatim copies
  of the archived files that produced the published tables.
- `tables/` — Python scripts that emit the paper's LaTeX table blocks from the
  summary CSVs.

## Pipeline order

Each pipeline folder runs datagen -> fit -> summary. Run every script with the
working directory set to its own folder (the staged copies replace the original
absolute paths with the working directory; each file's header comment records
this edit).

### Tb = 3 (`tb3/`)

1. Copy the model and helper in first (the fit script looks for them in its
   working directory):
   `mlGVARNoCorNoME_BurstIntercept_3Burst.txt` and `posteriorSummaryStats.R`
   from this folder into `tb3/`.
2. `mlGVARNoCorNoME_BurstIntercept_Flexible_DataGeneration.R` — simulates the
   nT = 9 data into `tb3/data/` (creates the folder itself), seeds
   `set.seed(r)` for r = 1..100.
3. `mlGVARNoCorNoME_BurstIntercept_Flexible_Modelfit.R` — fits the JAGS model
   to each replication in parallel, writing per-replication posterior summary
   tables and coda samples to `tb3/result/`.
4. `Generate_Summ_of_MCfile_withPower_fixed.R` — aggregates the 100
   per-replication tables of each cell into
   `mlGVARNoCorNoME_BurstIntercept_3Burst_3T_MCfileSumm_nT9_nP{100,500}.csv`
   (written inside `tb3/result/`).

### Tb = 5 and Tb = 20 (`tb5_tb20/`)

1. Copy the model and helper in first: `posteriorSummaryStats.R`, and
   `mlGVARNoCorNoME_BurstIntercept_3Burst.txt` renamed to
   `mlGVARNoCorNoME_BurstIntercept.txt` (the name the fit script expects; as
   noted above the content is byte-identical — it is the same 3-burst model).
2. Create the `tb5_tb20/data/` folder, then run
   `mlGVARNoCorNoME_BurstIntercept_CREqual0_DataGeneration.R` — simulates the
   nT = 15 and nT = 60 data, seeds `set.seed(r)` for r = 1..100.
3. `mlGVARNoCorNoME_BurstIntercept_CREqual0_Modelfit_Remaining.R` — the
   "_Remaining" in the name means it scans `result/` and fits whatever
   replications are missing. On a fresh folder that is all 400 (4 cells x 100
   replications), which is exactly how the published run used it.
4. `Generate_Summ_of_MCfile_withPower_fixed.R` — aggregates each cell into
   `mlGVARNoCorNoME_BurstIntercept_CREqual0_MCfileSumm_nT{15,60}_nP{100,500}.csv`
   (written inside `tb5_tb20/result/`). The four archived copies of these
   files are in `results/`.

### Tables (`tables/`)

Both scripts read the summary CSVs from `../results/`.

- `gen_table_latex.py` — writes the four full LaTeX table blocks:
  `table3_trend.tex` (`tab:study1_trend`, manuscript Table 2),
  `table2_restructured.tex` (`tab:study1_results`, manuscript Table 3),
  `trendsupp.tex` (`tab:study1_trend_supp`, Table B1), and
  `suppmetrics_restructured.tex` (`tab:study1_supp_metrics`, Table B2).
  The output filenames are historical; the `\label` in each block is the
  authoritative mapping.
- `extract_table_metrics.py` — prints a plain-text per-cell metrics report for
  the dynamic process parameters, used to check table values. Run it from the
  `tables/` directory.

## Compute environment

R 4.3.2, JAGS 4.3.1, rjags. R packages: rjags, coda, mvtnorm, Matrix, dplyr,
parallel, tidyverse, data.table. MCMC settings in the fit scripts: 2 chains,
5000 adaptation iterations, 5000 burn-in, 20000 sampling iterations, with
chain seeds (r, r + 500) for replication r.

## Runtime expectations

Data generation and summary aggregation take minutes. Model fitting is the
expensive step: each replication is a full JAGS run, and cost grows with N and
total T. With the parallel settings in the scripts (12 cores for Tb = 3,
4 cores for Tb = 5/20 — the lower count avoids memory exhaustion on the long
series), the Tb = 3 cells finish in hours while the full Tb = 5/20 grid takes
on the order of days. Each replication also saves its coda samples as .RData,
so `result/` grows large; budget disk space accordingly.

## Caveats

- Tb = 3 correction: the original 2026-05-09 Tb = 3 run generated the y2
  burst-intercept SDs as (1, 1.5, 2) instead of the designed (1, 2, 3). The
  `tb3/` scripts are the corrected 2026-09-04 versions (the fix is the single
  y2 sigma line, marked CORRECTED in the data-generation script) and are the
  versions whose results appear in the revised manuscript. The corrected
  Tb = 3 summary CSVs (`..._3Burst_3T_MCfileSumm_nT9_nP{100,500}.csv`) in
  `results/` come from the 2026-09-04 corrected re-run (200/200 replications;
  under the corrected design all burst-intercept SDs recover with |rBias|
  below 3% and coverage between .89 and .98 at Tb = 3). In short, the
  under-recovery originally reported for the two largest y2 SDs in this cell
  (relative bias of about -.25 and -.34 with near-zero coverage) was an
  artifact of summarizing the estimates against the intended rather than the
  actually generated values, not a property of the model. Tables regenerated
  from `results/` therefore reflect the corrected Tb = 3 values, which
  supersede the originally reported ones for the two affected y2 SD rows.
- All data in this study are simulated. No participant data are included or
  needed.
- Staged scripts are byte-identical to the versions that produced the
  results, except for a provenance header and the replacement of absolute
  working-directory paths; each file's header states exactly what changed.
