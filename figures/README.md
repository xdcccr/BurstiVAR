# figures/

This folder contains the scripts that produce the manuscript's two active figures. The manuscript (`main.tex`) includes the figure files by exactly these names: `figure1_horizontal_final_short.pdf` and `figure3_gohiar_groupmeans.pdf`.

## Figure-to-script map

| Manuscript figure | LaTeX label | Figure file | Producing script |
|---|---|---|---|
| Figure 1 | `fig:figure1` | `figure1_horizontal_final_short.pdf` | `figure1_v2_dotted_3person.py` |
| Figure 2 | `fig:gohiar_phase_means` | `figure3_gohiar_groupmeans.pdf` | `figure3_gohiar_groupmeans_plot.py` |

Note on naming: the second script and its output PDF carry "figure3" in their filenames for historical reasons. In the compiled manuscript the file is included as Figure 2. The filenames are kept as they are so they match the include line in `main.tex` and the original project files.

## Figure 1 (trend-dynamics confound illustration)

`figure1_v2_dotted_3person.py` is fully self-contained. It simulates its own data with fixed seeds. Panel A shows two 15-occasion series generated with opposing parameter configurations (slope beta = 0.42 with AR phi = 0.2, versus beta = 0.40 with phi = 0.5). Panel B shows three simulated persons, each with person-specific burst intercepts over 5 bursts of 10 occasions and a distinct across-burst trajectory. No input files are needed.

Run from this folder:

```
python figure1_v2_dotted_3person.py
```

It writes `figure1_horizontal_final_short.pdf` (plus a PNG preview, `fig1_v2_preview.png`) to the current working directory. Runtime is a few seconds.

## Figure 2 (phase-level posterior means by group)

`figure3_gohiar_groupmeans_plot.py` plots the reconstructed group-by-phase posterior means from the empirical intervention application (the CovFull model, with group covariates on the phase means and the VAR coefficients).

Input: `BurstiVAR_results_GroupCovFull/GoHIAR_BurstiVAR_4Burst_GroupCovFull_GroupMeans.csv`, resolved relative to the script's own location. This CSV is an aggregate posterior summary only: 16 rows (4 bursts x 2 variables x 2 groups) with columns Burst, Phase, Variable, Group, Mean, LL95, UL95. It contains no person-level data. It is written by `GoHIAR_BurstiVAR_4Burst_GroupCovFull_Modelfit.R` in the empirical model-fitting pipeline (R 4.3.2, JAGS 4.3.1, rjags; 2 chains with 5000 adaptation, 5000 burn-in, and 20000 sampling iterations). The CSV is not staged inside this folder. To rerun the script, obtain the CSV from the empirical-analysis part of this repository (or regenerate it by running the fitting pipeline, which requires the raw study data) and place it at the path above.

Run from this folder, with the CSV in place:

```
python figure3_gohiar_groupmeans_plot.py
```

It writes `figure3_gohiar_groupmeans.pdf` (plus a PNG preview, `fig3_new.png`) to the current working directory. Runtime is a few seconds.

## Environment

Python 3.12 with numpy and matplotlib. Figure 2's script uses only the standard library and matplotlib. Both scripts use the Agg backend, so no display is needed.

## Modifications relative to the originals

Each script carries a provenance header stating its original project path. The only edits are the header itself and the replacement of absolute output paths with relative filenames, so that the scripts run from this folder. All other content is identical to the versions that produced the published figures.
