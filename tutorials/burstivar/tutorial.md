# BurstiVAR tutorial

This tutorial shows how to simulate a small measurement-burst dataset, fit the
BurstiVAR model to it with JAGS, and read the output. All code below is
collected in `tutorial.R` in this folder, which runs top to bottom from this
directory. You need R (tested with 4.3.2), JAGS (tested with 4.3.1), and the R
packages `rjags`, `coda`, and `mvtnorm`. The whole script runs in about a
minute on a typical desktop (the MCMC itself took about 10 seconds on ours).

Two other files in this folder are copies of the actual files used to produce
the published Study 1 results: `mlGVARNoCorNoME_BurstIntercept_3Burst.txt`
(the JAGS model) and `posteriorSummaryStats.R` (the posterior summary
utility). Only a provenance header was added to each.

## 1. What BurstiVAR is

Measurement-burst designs collect short bursts of intensive measurements
(for example, a week of daily surveys), separated by long gaps (for example,
months). Two things happen on these two timescales. Within a burst, the
variables fluctuate around a stable level and influence each other from one
occasion to the next. Across bursts, the levels themselves can drift as people
develop.

BurstiVAR is a multilevel bivariate vector autoregressive (VAR) model built
for exactly this structure. It has two parts.

First, each person gets a separate intercept for each burst, on each
variable. These burst intercepts are piecewise constant: flat within a burst,
free to jump between bursts. Because every burst gets its own person-specific
level, any across-burst trend is absorbed without assuming a functional form.
The trend does not have to be linear, or Gompertz, or anything else you would
have to specify in advance.

Second, within a burst, the occasion-to-occasion dynamics operate on the
deviations from the person's current burst intercept. Each variable has a
person-specific autoregressive (AR) effect, its carryover from one occasion to
the next, and a person-specific cross-lagged (CR) effect, the influence of the
other variable at the previous occasion. The AR and CR coefficients are held
constant across bursts. The VAR process restarts at the beginning of each
burst, so nothing dynamic carries across the gap; only the person's
parameters do.

At the between-person level, all person-specific parameters (the burst
intercepts and the AR and CR coefficients) are drawn from normal
distributions. The model therefore estimates, for each parameter, a
group-level mean (`Level2Mean`) and a group-level SD (`Level2Sigma`), plus
the innovation covariance matrix (`sigma_innovation`).

## 2. Simulate a small burst dataset

We simulate data from the model itself, using the same generating logic as
the paper's Study 1 pipeline. The design: B = 3 bursts, Tb = 5 occasions per
burst, N = 50 persons, two variables (`y1`, `y2`). The generating values are
the paper's Study 1 values. Burst-intercept means step up across bursts, from
(0, 1) in burst 1 to (1, 1.5) in burst 2 and (2, 2) in burst 3, and their
between-person SDs grow from 1 to 1.5 to 2. The AR effects are .3 and .2, the
cross-lagged effect of `y1` on `y2` is -.15, the reverse effect is 0, and the
innovations have variance 1 with covariance .3.

```r
library(mvtnorm)
library(rjags)

nBurst       <- 3    # number of bursts (B)
nT_per_burst <- 5    # occasions per burst (Tb)
nT           <- nBurst * nT_per_burst   # total occasions per person (15)
nP           <- 50   # persons (N)
ydim         <- 2    # bivariate outcome

# Burst-intercept means and SDs (paper's Study 1 design)
Intercept_mu_list    <- list()
Intercept_sigma_list <- list()
for (b in 1:nBurst) {
  Intercept_mu_list[[b]]    <- c((b - 1) * 1.0, (b - 1) * 0.5 + 1.0)
  Intercept_sigma_list[[b]] <- c(1.0 + (b - 1) * 0.5, 1.0 + (b - 1) * 0.5)
}

# Within-burst dynamics
AR_mean1  <-  0.3    # autoregression of y1
AR_mean2  <-  0.2    # autoregression of y2
CR12_mean <- -0.15   # cross-lagged effect y1 -> y2
CR21_mean <-  0      # cross-lagged effect y2 -> y1

sigma_AR <- c(0.1, 0.1)   # between-person SDs of the AR coefficients
sigma_CR <- c(0.1, 0.1)   # between-person SDs of the CR coefficients

sigma_innovation  <- 1.0  # innovation SD (variance 1)
error_correlation <- 0.3  # innovation correlation (covariance 0.3)
```

The person-specific parameters live in one vector per person. The order
matters, because the JAGS model uses the same order: the six burst intercepts
first, then the four dynamic coefficients.

```r
# Parameter order (matches the JAGS model file):
#   [1:2] burst-1 intercepts, [3:4] burst-2, [5:6] burst-3,
#   [7] AR y1, [8] AR y2, [9] CR y1->y2, [10] CR y2->y1
n_intercept_params <- nBurst * ydim
n_total_params     <- n_intercept_params + 4

param_means <- c(unlist(Intercept_mu_list),
                 AR_mean1, AR_mean2, CR12_mean, CR21_mean)
param_sds   <- c(unlist(Intercept_sigma_list),
                 sigma_AR[1], sigma_AR[2], sigma_CR[1], sigma_CR[2])

# Person parameters are drawn independently, as in Study 1.
Sigma_full <- diag(param_sds) %*% diag(n_total_params) %*% diag(param_sds)

Sigma_innovation <- matrix(
  c(sigma_innovation^2,
    error_correlation * sigma_innovation^2,
    error_correlation * sigma_innovation^2,
    sigma_innovation^2),
  2, 2)
```

Now generate the data. For each person we draw their parameters, then walk
through the bursts. At the first occasion of a burst the series starts at the
burst intercept plus a fresh innovation. After that, each occasion's deviation
from the intercept follows the VAR(1) process.

```r
set.seed(1)

burst_boundaries <- seq(0, nT, by = nT_per_burst)   # 0, 5, 10, 15

Y       <- array(NA, c(nP, nT, ydim))
Y_trend <- array(NA, c(nP, nT, ydim))
Y_dev   <- array(NA, c(nP, nT, ydim))

person_params <- matrix(NA, nP, n_total_params)
for (i in 1:nP) {
  person_params[i, ] <- rmvnorm(1, param_means, Sigma_full)
}
AR1  <- person_params[, n_intercept_params + 1]
AR2  <- person_params[, n_intercept_params + 2]
CR12 <- person_params[, n_intercept_params + 3]
CR21 <- person_params[, n_intercept_params + 4]

for (i in 1:nP) {
  for (b in 1:nBurst) {
    T_start <- burst_boundaries[b] + 1
    T_end   <- burst_boundaries[b + 1]
    intercept_idx <- ((b - 1) * ydim + 1):(b * ydim)

    # First occasion of the burst: intercept + a fresh innovation.
    Y_trend[i, T_start, 1:ydim] <- person_params[i, intercept_idx]
    Y_dev[i, T_start, ] <- rmvnorm(1, mean = c(0, 0), sigma = Sigma_innovation)
    Y[i, T_start, ]     <- Y_dev[i, T_start, ] + Y_trend[i, T_start, ]

    # Remaining occasions: VAR(1) on the deviations.
    for (t in (T_start + 1):T_end) {
      Y_trend[i, t, 1:ydim] <- person_params[i, intercept_idx]
      mean_dev <- c(
        AR1[i] * Y_dev[i, t - 1, 1] + CR21[i] * Y_dev[i, t - 1, 2],  # y1
        AR2[i] * Y_dev[i, t - 1, 2] + CR12[i] * Y_dev[i, t - 1, 1]   # y2
      )
      Y_dev[i, t, ] <- rmvnorm(1, mean = mean_dev, sigma = Sigma_innovation)
      Y[i, t, ]     <- Y_trend[i, t, ] + Y_dev[i, t, ]
    }
  }
}

dat_long <- data.frame(
  id   = rep(1:nP, each = nT),
  time = rep(0:(nT - 1), times = nP),
  y1   = as.vector(t(Y[, , 1])),
  y2   = as.vector(t(Y[, , 2]))
)
write.csv(dat_long, "tutorial_simulated_data.csv", row.names = FALSE)
```

## 3. Fit the model with rjags

The model file `mlGVARNoCorNoME_BurstIntercept_3Burst.txt` in this folder is
the file that produced the published Study 1 three-burst results. It is worth
opening: Part 1 is the within-person model, written out burst by burst, Part 2
maps the ten person parameters onto the burst intercepts and the AR and CR
coefficients, and Part 3 gives the priors.

The model expects the data as an array `Y` (persons by occasions by 2),
the number of persons `P`, the cumulative burst boundaries `T1`, `T2`, `T3`
(here 5, 10, 15), and the Wishart prior inputs `W0` and `df_res` for the
innovation precision matrix.

```r
jags_data <- list(
  Y      = Y,
  P      = nP,
  T1     = as.integer(1 * nT_per_burst),
  T2     = as.integer(2 * nT_per_burst),
  T3     = as.integer(3 * nT_per_burst),
  W0     = diag(2),
  df_res = 3
)

parameters <- c("Level2Mean", "Level2Sigma", "sigma_innovation")

# Chain-specific RNG seeds, following the paper's (r, r + 500) convention.
inits <- list(
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 1),
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 501)
)

jagsModel <- jags.model("mlGVARNoCorNoME_BurstIntercept_3Burst.txt",
                        data = jags_data, inits = inits,
                        n.chains = 2, n.adapt = 1000)
update(jagsModel, n.iter = 1000)                       # burn-in
codaSamples <- coda.samples(jagsModel,
                            variable.names = parameters,
                            n.iter = 2000)
```

A note on the MCMC settings. The paper ran 2 chains with 5000 adaptation
iterations, 5000 burn-in iterations, and 20000 sampling iterations. Here we
use 1000 adaptation, 1000 burn-in, and 2000 sampling iterations, roughly an
eighth of the paper's run, so the fit finishes in seconds and you can
experiment freely. That is enough to see the model recover the generating
values on this small, well-behaved dataset. It is not enough for reporting
results. With short chains the effective sample sizes are smaller and the
interval endpoints less stable, so for a real analysis use the paper's
settings, which still only take minutes at this sample size, and confirm
convergence (R-hat close to 1, adequate ESS) before interpreting anything.

## 4. Summarize the posterior and read the output

`posteriorSummaryStats.R` provides `summarizePost()`, the same utility used
for the paper's results. It returns one row per monitored parameter with the
posterior mean, posterior SD, 95% credible interval (both equal-tailed
percentile and HDI versions), effective sample size, and R-hat.

```r
source("posteriorSummaryStats.R")
resulttable <- summarizePost(codaSamples)
write.csv(resulttable, "tutorial_result_summary.csv")
```

The rows map onto the model like this:

| Rows | Meaning | Generating value |
|---|---|---|
| `Level2Mean[1]`, `[2]` | burst-1 intercept means (y1, y2) | 0, 1 |
| `Level2Mean[3]`, `[4]` | burst-2 intercept means (y1, y2) | 1, 1.5 |
| `Level2Mean[5]`, `[6]` | burst-3 intercept means (y1, y2) | 2, 2 |
| `Level2Mean[7]`, `[8]` | AR means (y1, y2) | .3, .2 |
| `Level2Mean[9]` | CR mean, y1 -> y2 | -.15 |
| `Level2Mean[10]` | CR mean, y2 -> y1 | 0 |
| `Level2Sigma[1]`-`[6]` | between-person SDs of the burst intercepts | 1, 1, 1.5, 1.5, 2, 2 |
| `Level2Sigma[7]`-`[10]` | between-person SDs of AR and CR | .1 each |
| `sigma_innovation[.,.]` | innovation covariance matrix | var 1, cov .3 |

Here is the actual output of `tutorial.R` on this dataset (posterior mean,
posterior SD, 95% equal-tailed credible interval, ESS, and R-hat, next to the
generating value in the `true` column):

```
                       true   mean   PSD PCI 2.50% PCI 97.50%  ESS  RHAT
Level2Mean[1]          0.00  0.101 0.188    -0.271      0.470 2500 1.000
Level2Mean[2]          1.00  1.221 0.168     0.888      1.544 2641 1.000
Level2Mean[3]          1.00  0.936 0.197     0.550      1.328 2488 1.000
Level2Mean[4]          1.50  1.062 0.222     0.623      1.489 3176 1.000
Level2Mean[5]          2.00  2.073 0.358     1.356      2.773 3685 1.000
Level2Mean[6]          2.00  2.358 0.336     1.694      3.008 3646 1.002
Level2Mean[7]          0.30  0.261 0.073     0.123      0.407  100 1.000
Level2Mean[8]          0.20  0.174 0.061     0.047      0.290  151 1.000
Level2Mean[9]         -0.15 -0.197 0.058    -0.317     -0.088  129 1.002
Level2Mean[10]         0.00 -0.045 0.060    -0.161      0.073   94 1.002
Level2Sigma[1]         1.00  1.144 0.156     0.864      1.471 1358 1.000
Level2Sigma[2]         1.00  1.065 0.134     0.822      1.342 1582 1.000
Level2Sigma[3]         1.50  1.249 0.158     0.975      1.587 1719 1.000
Level2Sigma[4]         1.50  1.505 0.169     1.214      1.867 2589 1.000
Level2Sigma[5]         2.00  2.433 0.266     1.968      3.031 3052 1.000
Level2Sigma[6]         2.00  2.368 0.251     1.931      2.905 3203 1.000
Level2Sigma[7]         0.10  0.124 0.050     0.052      0.240  150 1.003
Level2Sigma[8]         0.10  0.123 0.050     0.052      0.243  136 1.004
Level2Sigma[9]         0.10  0.111 0.043     0.051      0.209  162 1.001
Level2Sigma[10]        0.10  0.095 0.034     0.047      0.182  176 1.000
sigma_innovation[1,1]  1.00  1.162 0.078     1.020      1.326  362 1.000
sigma_innovation[2,1]  0.30  0.306 0.052     0.207      0.409  487 1.001
sigma_innovation[1,2]  0.30  0.306 0.052     0.207      0.409  487 1.001
sigma_innovation[2,2]  1.00  1.001 0.064     0.884      1.136  917 1.000
```

How to read this. Take `Level2Mean[3]`, the burst-2 intercept mean of `y1`.
Its generating value is 1, and the posterior mean should land close to that,
with a 95% credible interval covering it. The interval means what it says in
plain Bayesian terms: given the model and the data, the parameter lies in
that range with 95% probability. The same reading applies to the AR rows
(`Level2Mean[7]`, `[8]`) and the CR rows (`[9]`, `[10]`). For `Level2Mean[10]`
the generating value is 0, so a well-behaved fit gives an interval that
straddles zero. The `Level2Sigma` rows tell you how much people differ: the
intercept SDs grow across bursts by design, and the AR and CR SDs are small
(.1), which makes them the hardest parameters to pin down at this sample
size.

The recovery here is good. The dynamic means come out close to the truth
(AR: .261 and .174 against .3 and .2; CR: -.197 against -.15), the interval
for `Level2Mean[9]` excludes zero as it should, and the interval for
`Level2Mean[10]` straddles zero as it should. The innovation covariance is
recovered almost exactly (.306 against .3). All R-hats are at or near 1.00
even at these reduced settings, though the ESS values for the dynamic
parameters (roughly 100 to 180) show why the full settings are needed for
reportable results.

Two rows deserve a closer look. The interval for `Level2Mean[4]` (posterior
mean 1.06) misses the population value of 1.5, and `sigma_innovation[1,1]`
(1.16) sits a bit above 1. Neither is model bias. With N = 50, the 50 persons
actually drawn need not match the population: for this seed, the realized
sample mean of the burst-2 intercept of `y2` is 1.06, and the model estimates
exactly that. The same holds elsewhere, for instance the realized sample SD
behind `Level2Sigma[5]` is 2.41 and the estimate is 2.43. A single small
replication tracks its own sample. The paper's Study 1 results average over
100 replications with the full MCMC settings and show good coverage of the
generating values.

## 5. Practical guidance from the paper's simulations

The paper's Study 1 varied burst length and sample size, and its conclusions
matter when you plan your own burst design.

Burst length has a lower bound. Tb = 3 is the floor, and it comes with
caveats: with only three occasions per burst, more replications fail to
converge and the dynamic parameters are estimated with visible bias and wide
intervals. Tb = 5 or more is the recommendation. At five occasions per burst
the dynamics are recovered well, and longer bursts mostly buy extra
precision.

More persons help, but they do not substitute for burst length. Raising N
tightens the credible intervals for every parameter, and the burst-intercept
means in particular benefit. Yet the information about each person's AR and
CR coefficients comes from the within-burst transitions, and no number of
additional persons creates transitions that were never observed. If you must
trade off, protect the number of occasions per burst first.

Finally, remember what the burst intercepts do and do not give you. They
absorb any across-burst trend without a functional form, which is the point
of the model. But they describe the trend only at the bursts you measured.
If you need to interpolate between bursts or extrapolate beyond them, you
need a growth model on top, and that is a modeling decision the data alone
will not make for you.
