# BurstiVAR

BurstiVAR is a multilevel vector autoregressive (VAR) model for measurement burst designs, in which participants complete short bursts of intensive longitudinal measurement separated by long gaps. The model estimates within-burst multilevel VAR dynamics (autoregression and cross-regression) while free person- and burst-specific intercepts absorb across-burst trends in a form-free way, so no parametric growth curve has to be assumed between bursts. GoBurstiVAR is a nested parametric variant that replaces the free burst intercepts with a Gompertz growth function, trading robustness for efficiency when the trend form is known. This repository contains the simulation code, the aggregate simulation results, the empirical analysis code, the figure scripts, and two runnable tutorials behind the paper.

## Citation

Xiong, X., Pesigan, I. J. A., Oravecz, Z., Lanza, S. T., & Chow, S.-M. (under review). *BurstiVAR for capturing multilevel vector autoregressive processes with trends in measurement burst designs.*

## Repository map

| Folder | Reproduces |
| --- | --- |
| `study1_burstivar/` | Study 1 parameter-recovery simulation: manuscript Tables 2 (`tab:study1_trend`) and 3 (`tab:study1_results`) plus the supplemental metrics table, across the Tb = 3 and Tb = 5/20 conditions |
| `study2_goburstivar/` | Study 2 nested GoBurstiVAR efficiency comparison: Table 4 (`tab:study2_var`), the trend and Gompertz recovery tables (Tables 5 and 6), and their supplement counterparts |
| `supplementD_twoburst/` | The two-burst (B = 2) supplemental condition: the B = 2 columns of Table F1 (`tab:supplement_2burst`) and the Supplement D convergence results |
| `convergence_diagnostics/` | Appendix Table A1 (`tab:convergence_summary`), the ESS/R-hat extraction across studies, and the Supplement D convergence audit |
| `empirical_gohiar/` | Empirical GoHiAR intervention analysis (code only, see Data availability): Table 7 (`tab:gohiar_sensitivity`), Tables D1/D2 (`tab:gohiar_covfull`, `tab:gohiar_covfull_re`), Table E1 (`tab:gohiar_trend_sens`), and the inputs to Figure 2 |
| `figures/` | Plotting scripts for Figure 1 (`fig:figure1`) and Figure 2 (`fig:gohiar_phase_means`) |
| `tutorials/burstivar/` | A small, runnable BurstiVAR example (simulate, fit, summarize) with validated output |
| `tutorials/goburstivar/` | A runnable GoBurstiVAR vs. BurstiVAR comparison on one simulated Gompertz-burst dataset |

Each folder has its own README with the exact pipeline, the table-label mapping, and folder-specific caveats.

## Environment

All results were produced with R 4.3.2, JAGS 4.3.1, and the `rjags` package. Unless a folder README says otherwise, MCMC used 2 chains with 5000 adaptation, 5000 burn-in, and 20000 sampling iterations, with chain seeds (r, r + 500) for replication r.

Runtime notes:

- **Studies 1 and 2.** Full simulation cells (200 replications each) take hours to days on a multicore machine. The fit scripts parallelize across replications. Keep the worker count at 8 or below on a 32 GB machine for the GoBurstiVAR fits, since more workers can exhaust memory.
- **Empirical analysis.** A single full-data fit takes roughly 40 minutes; the subsampled sensitivity fits are faster.
- **Tutorials.** The BurstiVAR tutorial runs in about a minute. The GoBurstiVAR tutorial runs in about five minutes. Both use reduced MCMC settings that their READMEs state explicitly.

## Reproduction quickstart

Every simulation pipeline follows the same order, run from inside the relevant folder:

1. **Data generation.** Run the `*_DataGeneration.R` script to simulate replications (seeded, so output is reproducible).
2. **Model fitting.** Run the `*_Modelfit.R` script to fit the JAGS model to each replication.
3. **Summary.** Run the `Generate_Summ_*.R` script to collapse posterior summaries into the `MCfile`/`MCfileSumm` CSVs.
4. **Tables.** Run the Python table scripts (in `study1_burstivar/tables/` and `study2_goburstivar/tables/`) to regenerate the LaTeX tables from those CSVs.

The scripts were patched to run from their own directories with relative paths; each carries a provenance header stating its original location and the exact lines changed. The staged `results/` CSVs are the aggregate outputs that produced the published tables, so steps 1–3 can be skipped if you only want to rebuild the tables.

## Data availability

All simulation data can be regenerated exactly from the seeded data-generation scripts, and the aggregate summary CSVs used for the published tables are included. The empirical intervention data are not publicly available; requests should be directed to the original study team. The `empirical_gohiar/` folder therefore ships code only, with the private input files replaced by documented `data_private/` placeholders. No participant-level data are included anywhere in this repository.

## Tutorials

Start with `tutorials/burstivar/` for a minimal end-to-end BurstiVAR analysis, then `tutorials/goburstivar/` to see how the parametric Gompertz variant compares against the free-intercept model on the same data. Both tutorials include their real, validated console output.
