# Study 2: BurstiVAR versus GoBurstiVAR (Gompertz burst trend)

This folder reproduces the Study 2 results of the BurstiVAR paper: the nested
comparison in which data are generated from GoBurstiVAR (a Gompertz function
of burst number for the burst intercepts, plus bivariate VAR(1) dynamics for
the within-burst deviations) and every simulated dataset is fit by both
BurstiVAR (free burst intercepts) and GoBurstiVAR (the true, correctly
specified model).

## Tables reproduced

Identify tables by their LaTeX `\label`; the generator's output filenames
carry an earlier in-manuscript numbering and do not match the final table
numbers.

| Manuscript table | Label | Generated file |
|---|---|---|
| Table 4 | `tab:study2_gompertz` | `table6_gompertz.tex` |
| Table 5 | `tab:study2_trend` | `table5_trend.tex` |
| Table 6 | `tab:study2_var` | `table4_var.tex` |
| Table C1 | `tab:study2_supp_gompertz` | `supp2_gompertz.tex` |
| Table C2 | `tab:study2_supp_trend` | `supp2_trend.tex` |
| Table C3 | `tab:study2_supp_var` | `supp2_var.tex` |

## Design

- B = 5 bursts in all conditions; 100 replications per cell.
- Cells, named by the total number of time points nT and the number of
  persons nP: `nT25_nP100` (T_b = 5, N = 100), `nT25_nP500` (T_b = 5,
  N = 500), `nT100_nP100` (T_b = 20, N = 100). Time points per burst are
  T_b = nT / 5, so nT = 25 gives T_b = 5 and nT = 100 gives T_b = 20.
- Generating VAR fixed effects: mu_phi11 = 0.30, mu_phi22 = 0.20,
  mu_phi21 = -0.15, mu_phi12 = -0.10; all VAR random-effect SDs 0.10;
  innovation variances 1.0 with covariance 0.3. Gompertz parameters as in
  the data-generation script (asymptote, displacement, growth rate for each
  of the two variables).

## Folder layout

- `code/` - data generation, the two JAGS model files, the two model-fit
  drivers, the posterior summary helper, and the three aggregation scripts.
- `results/MCfiles/` - archived per-cell aggregates from the original runs,
  verbatim. `*_MCfile_*.csv` holds the per-replication point estimates,
  posterior SDs, interval limits, and coverage/power flags for every
  monitored parameter; `*_MCfileSumm_*.csv` holds the across-replication
  summaries (rBias, RMSE, Coverage, RDSE, Power, and so on) that the tables
  are built from. Files beginning `mlGVARNoCorNoME_BurstIntercept_5Burst_on_`
  are the BurstiVAR fits, files beginning
  `mlGVARNoCorNoME_GompertzBurst_TrueModel_` the GoBurstiVAR fits.
- `results/TrendRecovery_B5/` - archived aggregates, verbatim, for the
  trend-recovery analysis: the BurstiVAR free burst intercepts evaluated
  against their Gompertz-implied true values (see Caveats).
- `tables/` - `gen_study2_tables.py` writes the six table blocks from
  `../results/`; `verify_study2_tables.py` re-checks every table cell
  against the same CSVs and exits nonzero on any mismatch.

All CSVs here are simulation summaries. No participant data of any kind are
included or needed.

## Pipeline order

All R scripts are run from `code/` (their paths are relative to that
directory). Steps 1-3 follow the author's per-cell workflow: each script was
run once per cell, editing only the nT/nP header lines between runs. The
staged scripts preserve the header values from the last run, so set nT/nP
(or nT_list/nP_list) to the cell you want before running.

1. `mlGVARNoCorNoME_GompertzBurst_DataGeneration.R` - generates
   `Data_*.csv` into `code/data/` for one (nT, nP) cell. Replication r uses
   `set.seed(r)`.
2. `mlGVARNoCorNoME_BurstIntercept_Modelfit.R` (BurstiVAR) and
   `mlGVARNoCorNoME_GompertzBurst_Modelfit.R` (GoBurstiVAR) - fit each
   replication with rjags, writing a per-replication result table (CSV) and
   the coda samples (RData) into `code/result/`. Both scripts skip
   replications whose result table already exists, so they can be re-run to
   fill in failures.
3. `Generate_Summ_GompertzBurst_BurstIntercept.R` and
   `Generate_Summ_GompertzBurst_TrueModel.R` - aggregate `code/result/` into
   the MCfile and MCfileSumm CSVs (written into `code/result/`; the archived
   copies live in `results/MCfiles/`).
4. `Generate_Summ_BurstIntercept_5Burst_TrendRecovery.R` - loops over all
   three cells itself; reads the BurstiVAR result tables from `code/result/`
   and writes the trend-recovery MCfile/MCfileSumm CSVs to
   `code/TrendRecovery_B5/` (archived copies in `results/TrendRecovery_B5/`).
5. `tables/gen_study2_tables.py` then `tables/verify_study2_tables.py` - run
   from `tables/`. These read the archived CSVs in `../results/`, so they
   work without re-running any simulation. At staging time the verification
   passed (764 cells, 0 mismatches) and the six regenerated `.tex` blocks
   were byte-identical to the published versions.

Regenerated summaries land in `code/result/` and `code/TrendRecovery_B5/`,
never in `results/`, so a re-run cannot overwrite the archived numbers.

## Compute environment and runtime

- R 4.3.2, JAGS 4.3.1, rjags (plus coda, mvtnorm, parallel; the aggregation
  scripts also use dplyr/tidyverse/data.table).
- MCMC per replication: 2 chains, 5000 adaptation + 5000 burn-in + 20000
  sampling iterations; chain RNG seeds (r, r + 500) where r is the
  replication number.
- Data generation and the aggregation/table scripts run in minutes or less.
  Model fitting is the heavy step: 100 replications per cell per model on a
  handful of parallel workers takes on the order of hours to a day per cell,
  longer for nT = 100 or nP = 500.

## Caveats

1. GoBurstiVAR fitting is memory-heavy. On a 32 GB machine use at most 8
   parallel workers; a 16-worker run crashed silently (out of memory). The
   staged scripts are set to `num_cores <- 6`.
2. A fourth cell (T_b = 20, N = 500; i.e. nT100_nP500) was generated but
   never fitted or reported. The preserved loop headers in the
   data-generation and GoBurstiVAR model-fit scripts reflect that
   last-edited state.
3. The T_b = 5 versus T_b = 20 distinction comes entirely from
   T_per_burst = nT / 5 inside the scripts; there is no separate T_b
   setting.
4. Trend recovery (Table 5 / `tab:study2_trend` and C2): under the Gompertz
   data-generating process, BurstiVAR's free burst intercepts have
   well-defined estimands - the population mean and between-person SD of
   each burst's intercept. `Generate_Summ_BurstIntercept_5Burst_TrendRecovery.R`
   derives these true values by Monte Carlo integration with 2,000,000 draws
   (seed 20260716) over the person-level Gompertz parameter distribution
   used in data generation.
5. In `verify_study2_tables.py`, the regression checks against the old
   manuscript tables are disabled (`RUN_OLD_REGRESSION = False`); they
   would require the manuscript source, which is not distributed.

Each staged script begins with a provenance header stating its original
location in the project archive and exactly what was changed for this
repository (path configuration only). Everything else is identical to the
versions that produced the published results.
