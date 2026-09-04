# GoBurstiVAR tutorial

A hands-on tutorial for GoBurstiVAR, the parametric variant of BurstiVAR in
which the free burst-level means are replaced by a person-specific Gompertz
growth curve over burst number (asymptote, displacement, and growth rate per
variable). Work through the BurstiVAR tutorial in `../burstivar/` first; this
one builds on it.

This folder does not reproduce a paper table. It demonstrates, at small
scale, the model comparison behind the paper's Study 2 (the table labeled
`tab:study2_var`): data are generated from GoBurstiVAR and the same dataset
is fit with both GoBurstiVAR and BurstiVAR. The full 100-replication Study 2
pipeline that produces the published numbers lives in
`../../study2_goburstivar/`.

## Files

- `tutorial.md` — the tutorial text, with real output from a validated run
  of `tutorial.R` pasted in.
- `tutorial.R` — the companion script; runs everything end to end.
- `mlGVARNoCorNoME_GompertzBurst.txt` — GoBurstiVAR JAGS model (the
  published Study 2 version; provenance header at top).
- `mlGVARNoCorNoME_BurstIntercept_5Burst.txt` — BurstiVAR free-intercept
  JAGS model for B = 5 (the published Study 2 version; provenance header at
  top).
- `posteriorSummaryStats.R` — posterior summary utilities (`summarizePost`),
  as used throughout the paper's simulations.

## Pipeline order

`tutorial.R` runs the whole pipeline in one pass, from this directory:

1. **Data generation** — one Gompertz-burst dataset with B = 5 bursts,
   T_b = 5 time points per burst, N = 40 persons (generating logic reused
   from the Study 2 data-generation script).
2. **Fit 1** — the true model, GoBurstiVAR
   (`mlGVARNoCorNoME_GompertzBurst.txt`).
3. **Fit 2** — the free-intercept model, BurstiVAR
   (`mlGVARNoCorNoME_BurstIntercept_5Burst.txt`), on the same dataset.
4. **Summary** — side-by-side comparison of the VAR estimates and posterior
   SDs, Gompertz parameter recovery, burst-mean recovery, and convergence
   diagnostics.

Run it as

```
Rscript tutorial.R
```

with this folder as the working directory. The script writes its outputs
(`tutorial_data_gompertzburst.csv`, `tutorial_fit_goburstivar.csv`,
`tutorial_fit_burstivar.csv`, `tutorial_var_comparison.csv`) into the same
folder; they are generated at run time and are not stored in the repository.

## Runtime expectations

With the reduced tutorial MCMC (2 chains, 1000 adaptation + 1000 burn-in +
2000 sampling iterations), the validated run on a Windows desktop took about
4.3 minutes for the GoBurstiVAR fit and about 15 seconds for the BurstiVAR
fit, roughly 5 minutes end to end. The nonlinear Gompertz
likelihood makes GoBurstiVAR clearly slower and more memory-hungry than
BurstiVAR. At paper scale (N = 500, full MCMC below) a single GoBurstiVAR
replication takes several minutes to tens of minutes, and on a 32 GB machine
at most 8 parallel workers should be used for the Gompertz model (16 workers
exhausted RAM).

## Compute environment

Published results were produced with R 4.3.2, JAGS 4.3.1, and rjags. Paper
MCMC settings: 2 chains, 5000 adaptation + 5000 burn-in + 20000 sampling
iterations, chain seeds (r, r + 500) for replication r. The tutorial reduces
this to 1000 / 1000 / 2000 with seeds (1, 501); use the paper settings for
any real analysis.

## Caveats

- The tutorial's N = 40 and reduced MCMC are for demonstration. Single-run
  output is noisy; the paper's efficiency findings rest on 100 replications
  per cell.
- With the reduced adaptation phase JAGS may warn that adaptation is
  incomplete. For the demonstration this is acceptable.
- With only 5 bursts, the between-person SDs of the Gompertz parameters are
  weakly identified; expect wide posteriors on those even in the true model.
