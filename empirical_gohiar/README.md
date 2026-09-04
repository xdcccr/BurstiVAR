# empirical_gohiar: BurstiVAR applied to the Go-HIAR intervention study

This folder contains the analysis code for the empirical application of BurstiVAR
reported in the paper. It reproduces:

- Table 7 (`tab:gohiar_sensitivity`): population VAR coefficients, innovation
  (co)variances, and MCMC diagnostics across the full design, three within-burst
  sparsity subsamples, the daily-sampling variant, and the fit without burst
  structure.
- Table D1 (`tab:gohiar_covfull`): full posterior summary of the fixed effects
  in the group-covariate model.
- Table D2 (`tab:gohiar_covfull_re`): random-effect SDs and innovation
  (co)variances of the group-covariate model.
- Table E1 (`tab:gohiar_trend_sens`): the trend side of the same auxiliary fits
  as Table 7 (pooled baseline means, phase offsets, between-person SDs of the
  burst intercepts).
- The phase-means figure (`fig:gohiar_phase_means`): the group-covariate fit
  script exports `GoHIAR_BurstiVAR_4Burst_GroupCovFull_GroupMeans.csv`, the
  reconstructed group-by-burst posterior means that the figure plots. The
  plotting script itself is not part of this folder.

## Data availability (code only)

The intervention data are private. This folder ships no participant data.

The scripts read four pre-processed person-level R data files:
`datblock_Meaning.Rdata`, `datblock_relationship.Rdata`,
`datblock_accomplishment.Rdata`, and `datCov.Rdata`. These files were provided
by the original study team and were prepared upstream in a companion project;
no builder script ships here. Requests for the data should be directed to the
original study team. To run the scripts, place the four files in a
subdirectory named `data_private/` inside this folder.

The analytic sample (N = 108; 54 Control, 54 PPI+Meditation) arrives
pre-restricted inside those inputs. The scripts perform no additional
participant selection beyond what is visible in the code (the group-covariate
script keeps groups 1 and 3, which together are the same 108 persons).

## Design facts encoded in the scripts

- B = 4 study phases (bursts): P1 Pre-intervention, P2 Intervention,
  P3 Post-intervention I, P4 Post-intervention II.
- Four 6-hour blocks per day; burst boundaries at cumulative slots
  56 / 116 / 168 / 224 (so T = 224 slots in the full design).
- Bivariate completeness: any occasion with only one of the two variables
  (Meaning of life, Relationship quality) observed is set fully missing;
  JAGS handles the missing values.
- MCMC seeds 2026 and 2027 (one per chain), fixed in every script.

## Which script backs which table column

All six fits are pooled (no group covariate) except the last row. The columns
refer to Tables 7 and E1, which report the same six fits.

| Column / result | Fit script | JAGS model |
|---|---|---|
| Full | `GoHIAR_BurstiVAR_4Burst_BaselineDelta_Modelfit.R` | `GoHIAR_BurstiVAR_4Burst_BaselineDelta.txt` |
| Tb = 20 | `GoHIAR_BurstiVAR_4Burst_BaselineDelta_Subsampled_Tb20_Modelfit.R` | same |
| Tb = 5 and Tb = 3 | `GoHIAR_BurstiVAR_4Burst_BaselineDelta_Subsampled_Modelfit.R` (fits both conditions in one run) | same |
| Daily | `GoHIAR_BurstiVAR_4Burst_BaselineDelta_Subsampled_Daily5_Modelfit.R` (5 occasions per burst, 24-hour lag) | same |
| No burst | `GoHIAR_mlVAR_NoBurst_Modelfit.R` | `GoHIAR_mlVAR_NoBurst.txt` |
| Tables D1, D2, and the phase-means figure | `GoHIAR_BurstiVAR_4Burst_GroupCovFull_Modelfit.R` | `GoHIAR_BurstiVAR_4Burst_GroupCovFull.txt` |

`posteriorSummaryStats.R` (helper by Zita Oravecz, with credits to Kruschke's
DBDA2E utilities) provides `summarizePost()`, which every fit script uses to
build its posterior summary table.

`GoHIAR_BurstiVAR_4Burst_GroupCovFull_descriptives.csv` is the aggregate
descriptives table written by the group-covariate script: per burst and group
(Pooled / Control / Intervention), the mean, median, min, max, and SD across
persons of the number of bivariate-complete observations. It contains
group-level summaries only, no person-level rows, and is included here so the
observed-data description can be checked without access to the raw data.

## Pipeline

There is no data-generation step in this folder; the inputs are the private
files described above.

1. Place the four `.Rdata` files in `data_private/`.
2. Run any fit script from this directory (each is standalone), e.g.
   `Rscript GoHIAR_BurstiVAR_4Burst_BaselineDelta_Modelfit.R`.
   Each script builds its Y array from the inputs, runs JAGS, and writes a
   results folder next to the scripts (`BurstiVAR_results_BaselineDelta/`,
   `BurstiVAR_results_BaselineDelta_Subsampled/`, `BurstiVAR_results_NoBurst/`,
   or `BurstiVAR_results_GroupCovFull/`) containing the posterior summary
   table (`*_resulttable.csv`), the raw MCMC samples (`*_codaSamples.RData`),
   the Y array, and fit-specific CSV exports (SD-difference tables, the
   descriptives table, group means, and group dynamics for the covariate fit).
3. The paper tables are read off the `*_resulttable.csv` files (posterior
   means and 95% equal-tailed credible intervals) plus the printed R-hat and
   effective-sample-size diagnostics in the console log.

The within-burst subsampling in the Tb20 / Tb5+Tb3 / Daily5 scripts is
deterministic (a person-specific sliding-window rule with fixed tie-breaking,
no random number draws), so the subsampled fits are exactly reproducible from
the same inputs.

## Compute environment and runtime

R 4.3.2, JAGS 4.3.1, rjags. Every fit uses 2 chains with 5000 adaptation,
5000 burn-in, and 20000 sampling iterations (thin = 1), chain seeds 2026 and
2027.

On a modern desktop, the full-data no-burst fit took about 40 minutes; the
full-data BurstiVAR fits (224 slots, N = 108, with missing-value imputation)
are of the same order or somewhat longer. The subsampled fits (12 to 80 time
points per person) run in minutes to tens of minutes each. The
Tb5+Tb3 script runs two fits back to back.

## Caveats

- Exact numerical reproduction requires the private input files; with them,
  the fixed seeds make every fit reproducible to the digit reported.
- The scripts contain comments in Chinese alongside English; the code and all
  parameter names are English.
- The commented-out "Roar Cluster" path blocks inside some scripts document
  the original computing environment and can be ignored.
