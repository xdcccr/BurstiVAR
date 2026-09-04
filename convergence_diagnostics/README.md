# convergence_diagnostics

This folder reproduces the MCMC convergence diagnostics reported in the paper:
Appendix Table A1 (`tab:convergence_summary`) and the diagnostics prose in
Supplement D (the two-burst minimum design).

## Contents

Scripts (each carries a provenance header describing its exact modifications;
all are otherwise identical to the versions that produced the published results):

- `Extract_ESS_RHAT_Summary.R` — Round 1 extraction. Reads the per-replication
  `*_resulttable_*.csv` files written by the simulation fits, collects ESS and
  R-hat for every parameter, and writes per-condition raw and summary CSVs plus
  `GRAND_SUMMARY_ESS_RHAT.csv`.
- `Extract_ESS_RHAT_Round2.R` — Round 2 extraction for the five conditions
  Round 1 missed. Writes `GRAND_SUMMARY_ESS_RHAT_Round2.csv`.
- `CORRECTED_audit_Study1_Tb3_B3_from_MLGVAR2605.py` — the correct audit of the
  Study 1 Tb = 3 cell (B = 3 bursts, T = 9). Prints its overview to stdout; it
  writes no files.
- `_verify_2burst_for_suppD.py` — recomputes the Supplement D two-burst numbers
  from the raw CSVs in this folder, as a check on the quoted values.

Data files:

- `GRAND_SUMMARY_ESS_RHAT.csv` — Round 1 condition-level overview (verbatim
  archive copy).
- `GRAND_SUMMARY_ESS_RHAT_Round2.csv` — Round 2 condition-level overview
  (verbatim archive copy; see the legend below about its `Study1_Tb3` rows).
- `SuppD_TwoBurst_N100_ESS_RHAT_raw.csv`, `SuppD_TwoBurst_N100_ESS_RHAT_summary.csv`,
  `SuppD_TwoBurst_N500_ESS_RHAT_raw.csv`, `SuppD_TwoBurst_N500_ESS_RHAT_summary.csv`
  — per-replication ESS/R-hat records and per-parameter summaries for the
  two-burst design (verbatim copies, renamed; see the legend below).

## Legend: the Tb = 3 / two-burst mislabeling

Read this before using the `Study1_Tb3` labels anywhere in these files.

1. The four `SuppD_TwoBurst_*` CSVs, and the two `Study1_Tb3_*` rows inside
   `GRAND_SUMMARY_ESS_RHAT_Round2.csv`, were computed from the B = 2 two-burst
   minimum design (file prefix `..._2Burst_3T`, nT = 6) reported in
   Supplement D — not from the B = 3 Study 1 Tb = 3 cell. The original archive
   filenames said `Study1_Tb3` and were mislabeled. The repo copies of the four
   CSVs are renamed accordingly; the GRAND CSV is kept verbatim for provenance,
   with this caveat standing in place of a rename. The `Condition` column
   inside the renamed raw CSVs likewise still reads `Study1_Tb3_*`, because the
   file contents are untouched.
2. The true B = 3, Tb = 3 audit is
   `CORRECTED_audit_Study1_Tb3_B3_from_MLGVAR2605.py`. It reads the
   per-replication resulttables from the Tb = 3 re-run (MLGVAR2605) and prints
   its condition overviews to stdout.
3. In Round 1, the `Study1_Tb3` block of `Extract_ESS_RHAT_Summary.R` pointed
   at a folder that did not exist, so it produced no output. That is why
   `GRAND_SUMMARY_ESS_RHAT.csv` has no Tb = 3 rows. Only its Study 1 Tb = 5 and
   Tb = 20 rows and its Study 2 GoBurstiVAR rows are valid.
4. Round 2 contributed the Study 2 BurstiVAR rows (and the mislabeled
   `Study1_Tb3` rows described in point 1).
5. The Tb = 3 diagnostics in the manuscript will be recomputed from the
   corrected re-run at revision.

So, to assemble Table A1: Study 1 Tb = 5 and Tb = 20 and Study 2 GoBurstiVAR
come from `GRAND_SUMMARY_ESS_RHAT.csv`; Study 2 BurstiVAR comes from
`GRAND_SUMMARY_ESS_RHAT_Round2.csv`; Study 1 Tb = 3 comes from the stdout of
the corrected Python audit. The Supplement D two-burst prose numbers come from
the `SuppD_TwoBurst_*` files (checked by `_verify_2burst_for_suppD.py`).

## Pipeline order

These scripts sit at the end of the pipeline. Upstream, in the study folders of
this repository, data generation and model fitting produce one
`*_resulttable_nT<T>_nP<N>_r<r>.csv` per replication; each resulttable already
contains an ESS column and an RHAT column per parameter. This folder only
aggregates them:

1. Run the simulations (other repo folders) to produce the per-replication
   resulttables.
2. `Extract_ESS_RHAT_Summary.R`, then `Extract_ESS_RHAT_Round2.R` — run from
   this directory, with the `./MLGVAR2601/...` result folders (as named in the
   scripts' path-configuration blocks) pointing at your local simulation
   output.
3. `CORRECTED_audit_Study1_Tb3_B3_from_MLGVAR2605.py` — run from this
   directory, with `./MLGVAR2605/result` holding the Tb = 3 (B = 3, T = 9)
   resulttables.
4. `_verify_2burst_for_suppD.py` — run from this directory as staged; it reads
   the four `SuppD_TwoBurst_*_raw.csv` files shipped here.

## Runtime

The extraction and audit scripts only read small CSVs; each finishes in a few
minutes at most. Step 4 runs in seconds on the shipped files. The expensive
part is the upstream MCMC fitting, documented in the study folders.

## Compute environment

R 4.3.2, JAGS 4.3.1, rjags. All fits used 2 chains with 5000 adaptation,
5000 burn-in, and 20000 sampling iterations; replication r used seeds
(r, r + 500). The two Python scripts need only the standard library
(Python 3.8+).
