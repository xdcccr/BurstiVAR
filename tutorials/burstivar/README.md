# BurstiVAR tutorial

This folder contains the tutorial promised in the paper's Transparency and
Openness statement. It walks an applied researcher through simulating a small
measurement-burst dataset, fitting the BurstiVAR model with JAGS, and reading
the posterior summary. Start with `tutorial.md`; the same code is collected in
`tutorial.R`, which runs top to bottom from this directory.

## Files

- `tutorial.Rmd` — the tutorial source (R Markdown); knitting it runs all
  code and embeds the output.
- `tutorial.md` — knitted from `tutorial.Rmd`
  (`knitr::knit("tutorial.Rmd", output = "tutorial.md")`), so every number
  shown comes from an actual run. Read this one on GitHub.
- `tutorial.html` — the same knitted tutorial as a standalone, self-contained
  HTML page (GitHub displays HTML files as source; download it and open in a
  browser).
- `tutorial.R` — the same code as a runnable script, for readers who prefer
  to run and modify it directly.
- `mlGVARNoCorNoME_BurstIntercept_3Burst.txt` — the JAGS model file. This is
  the actual file that produced the published Study 1 three-burst results
  (only a provenance header was added).
- `posteriorSummaryStats.R` — the posterior summary utility used throughout
  the paper (only a provenance header was added).

Running `tutorial.R` writes two files into this directory:
`tutorial_simulated_data.csv` (the simulated dataset) and
`tutorial_result_summary.csv` (the posterior summary table).

## Relation to the paper

The tutorial does not reproduce any paper table. It demonstrates, at tutorial
scale, the model and pipeline behind the Study 1 results (Tables
`\label{tab:study1_results}`, `\label{tab:study1_trend}`, and
`\label{tab:study1_supp_metrics}` in the manuscript): the same generating
values, the same JAGS model file, and the same summary utility, but with a
single replication of N = 50 persons, B = 3 bursts of Tb = 5 occasions, and
reduced MCMC settings. The folders for the simulation studies elsewhere in
this repository reproduce the tables themselves.

## Pipeline order

1. Data generation (Section 2 of `tutorial.R`): draw person-specific
   parameters, then simulate the burst time series.
2. Model fit (Sections 3 and 4): fit the JAGS model with rjags.
3. Summary (Section 5): `summarizePost()` from `posteriorSummaryStats.R`,
   printed next to the generating values and written to CSV.

## Runtime and settings

The script takes about a minute on a typical desktop; the MCMC itself took
about 10 seconds in our validation run (R 4.3.2, JAGS 4.3.1, rjags).

The tutorial deliberately uses REDUCED MCMC settings: 2 chains, 1000
adaptation, 1000 burn-in, 2000 sampling iterations. The paper's results used
2 chains with 5000 adaptation, 5000 burn-in, and 20000 sampling iterations,
with chain seeds (r, r + 500) for replication r. Use the paper's settings for
any real analysis; the reduced settings are for a fast first contact with the
model.

## Caveats

- The tutorial fits one simulated dataset with N = 50. Its estimates track
  the realized sample of 50 persons, so an interval can miss the population
  generating value without indicating model bias. `tutorial.md` discusses two
  such rows in the validated output.
- The reduced MCMC settings yield modest effective sample sizes for the
  dynamic parameters. They are fine for the tutorial, not for reporting.
- The design advice in the tutorial's final section summarizes the paper's
  Study 1 findings: Tb = 3 is a floor with caveats, Tb of 5 or more is
  recommended, and a larger N raises precision but cannot substitute for
  burst length.
