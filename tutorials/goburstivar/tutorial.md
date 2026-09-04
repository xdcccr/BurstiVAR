# GoBurstiVAR: a hands-on tutorial

This tutorial shows how to fit **GoBurstiVAR**, the parametric variant of
BurstiVAR in which the free burst-level means are replaced by a person-specific
Gompertz growth curve over burst number. It assumes you have already worked
through the **BurstiVAR tutorial** in `../burstivar/`, which introduces the
measurement-burst setup, the free-intercept JAGS model, and the fitting
workflow. Here we build directly on that material: we simulate a small dataset
whose burst means truly follow a Gompertz curve, fit it **twice** (once with
GoBurstiVAR, once with BurstiVAR), and compare what the restriction buys and
what it risks.

Everything below is run by the companion script `tutorial.R` in this folder.
The two JAGS model files and `posteriorSummaryStats.R` are the same files used
for the published Study 2 results (see the provenance headers). All output
shown in this document is real output from running `tutorial.R` with the
settings described here.

## 1. What GoBurstiVAR is

BurstiVAR treats the burst-level means as **free parameters**: with B bursts
and two variables there are 2B fixed effects, plus a person-specific random
intercept for each burst and variable. Nothing constrains how the means move
from burst to burst. That is the model's strength: it imposes no functional
form on developmental change across bursts.

GoBurstiVAR replaces those free means with a **Gompertz growth curve over
burst number**. For person *i*, variable *v*, and burst *b* = 1, ..., B, the
burst-level mean is

```
trend_iv(b) = theta1_iv * exp( -theta2_iv * exp( -b * theta3_iv ) )
```

with three person-specific parameters per variable:

- `theta1` — the **asymptote**, the level growth approaches;
- `theta2` — the **displacement**, which controls where on the curve the
  process starts;
- `theta3` — the **growth rate**, how fast the process approaches the
  asymptote.

Each person draws their own triple from a normal population distribution, so
the model estimates a population mean and a between-person SD for all six
parameters (three per variable). Everything else — the within-burst VAR(1) on
the deviations, the deviation reset at each burst start, the innovation
covariance matrix — is identical to BurstiVAR.

When the random-effect structures are taken as independent, GoBurstiVAR is
**nested** under BurstiVAR's intercept structure: the Gompertz curve
constrains the 2B free burst means (10 here, with B = 5) to lie on a
3-parameter curve per variable (6 parameters). The paper's Study 2 asks what
that restriction is worth. The answer is an asymmetry:

- **When the Gompertz form holds**, both models are correct and the
  restriction buys efficiency. In the paper, the clearest gain was power on
  the small cross-lagged effect `mu_phi12 = -0.10` at T_b = 5, N = 100:
  **.98 for GoBurstiVAR versus .84 for BurstiVAR**. The gap closed at N = 500
  or T_b = 20.
- **When the Gompertz form fails**, GoBurstiVAR is misspecified, while
  BurstiVAR — which assumes no form at all — remains correct.

So BurstiVAR is the safe default, and GoBurstiVAR is the sharper tool for the
case where you have good reason to believe in smooth, monotone,
asymptoting growth across bursts. Section 5 returns to this choice.

## 2. Simulating a Gompertz-burst dataset

We generate one dataset with **B = 5 bursts, T_b = 5 time points per burst
(25 in total), and N = 40 persons**. Forty persons is deliberately small so
the whole tutorial runs in minutes; the paper's Study 2 used N = 100 and 500
with 100 replications per cell. The generating logic is reused from the
Study 2 data-generation script
(`paper_simulation_archive/Study2_Gompertz/code/mlGVARNoCorNoME_GompertzBurst_DataGeneration.R`).

Generating values:

| Parameter | y1 | y2 |
|---|---|---|
| Gompertz mean (asymptote, displacement, growth rate) | 10, 3, 1.0 | 8, 2.5, 0.8 |
| Gompertz between-person SDs | 2, 0.5, 0.2 | 1.5, 0.4, 0.15 |
| Autoregression mean (`mu_phi11`, `mu_phi22`) | 0.30 | 0.20 |
| Cross-lag means | `mu_phi21` (y1 -> y2) = -0.15 | `mu_phi12` (y2 -> y1) = -0.10 |
| SD of every VAR coefficient | 0.10 | 0.10 |
| Innovation SD / correlation | 1.0 / 0.3 | 1.0 / 0.3 |

These are the paper's Study 2 values, including the small cross-lag
`mu_phi12 = -0.10` on which the power difference showed up. Generation works
exactly as in the BurstiVAR tutorial except for where the burst means come
from: for each person we first evaluate their Gompertz curve at b = 1, ..., 5
to get the burst means, then run the person's VAR(1) on the deviations within
each burst, resetting the deviation process at every burst start.

Running Step 1 of `tutorial.R` prints the burst-level means implied by the
mean Gompertz parameters — the growth pattern buried in the data:

```
Simulated dataset: 40 persons x 25 time points ( 5 bursts x 5 points ).
Population burst-level means implied by the mean Gompertz parameters:
    [,1]  [,2]  [,3]  [,4]  [,5]
y1 3.317 6.663 8.613 9.465 9.800
y2 2.602 4.829 6.377 7.225 7.642
```

(Columns 1 through 5 are bursts 1 through 5.)

Both variables rise steeply over the first two or three bursts and then level
off toward their asymptotes, y1 toward 10 and y2 toward 8. This is the
signature Gompertz shape the parametric model will try to exploit.

## 3. Fitting the same dataset twice

Step 2 and Step 3 of `tutorial.R` fit the same `Y` array with two JAGS models:

1. **GoBurstiVAR** (`mlGVARNoCorNoME_GompertzBurst.txt`) — the true model
   here. It estimates 10 level-2 means: the six Gompertz parameters, then the
   four VAR coefficients.
2. **BurstiVAR** (`mlGVARNoCorNoME_BurstIntercept_5Burst.txt`) — the
   free-intercept model from the BurstiVAR tutorial, with 14 level-2 means:
   ten burst intercepts, then the four VAR coefficients.

Both take the same data list (`Y`, `P`, the Wishart prior inputs, and the
burst boundaries `T1`–`T5`), and we monitor `Level2Mean`, `Level2Sigma`, and
`sigma_innovation` in both. Because the two models order their parameters
differently, keep these index maps at hand when reading the output:

| GoBurstiVAR | meaning | true value |
|---|---|---|
| `Level2Mean[1:3]` | y1 asymptote, displacement, growth rate | 10, 3, 1.0 |
| `Level2Mean[4:6]` | y2 asymptote, displacement, growth rate | 8, 2.5, 0.8 |
| `Level2Mean[7]`, `[8]` | `mu_phi11`, `mu_phi22` | 0.30, 0.20 |
| `Level2Mean[9]`, `[10]` | `mu_phi21` (y1 -> y2), `mu_phi12` (y2 -> y1) | -0.15, -0.10 |

| BurstiVAR | meaning | true value |
|---|---|---|
| `Level2Mean[1:10]` | burst intercepts (b1 y1, b1 y2, ..., b5 y1, b5 y2) | Gompertz-implied |
| `Level2Mean[11]`, `[12]` | `mu_phi11`, `mu_phi22` | 0.30, 0.20 |
| `Level2Mean[13]`, `[14]` | `mu_phi21` (y1 -> y2), `mu_phi12` (y2 -> y1) | -0.15, -0.10 |

**MCMC settings.** For a tutorial-sized run we use 2 chains with 1000
adaptation, 1000 burn-in, and 2000 sampling iterations (chain seeds 1 and
501). The paper used 2 chains with 5000 adaptation, 5000 burn-in, and 20000
sampling iterations, seeds (r, r + 500) for replication r, under R 4.3.2 /
JAGS 4.3.1 / rjags. Use the paper's settings for anything beyond a
demonstration. With the reduced adaptation JAGS may warn that adaptation is
incomplete; for this demonstration that is acceptable.

**Runtime and memory warning.** The Gompertz model evaluates a nonlinear
curve inside the likelihood at every node, so it is noticeably **slower and
more memory-hungry** than the free-intercept model. Timings from the run
behind this document (Windows desktop, R 4.3.2, JAGS 4.3.1):

```
GoBurstiVAR fit took 4.3 minutes.
BurstiVAR fit took 0.2 minutes.
```

The free-intercept model finished in about 15 seconds; the Gompertz model
needed over 4 minutes on the same data with the same MCMC settings, roughly a
factor of 20.

At paper scale (N = 500, full MCMC) a single GoBurstiVAR replication takes on
the order of several minutes to tens of minutes, and the memory footprint adds
up across parallel workers: on a 32 GB machine, running 16 parallel
GoBurstiVAR fits exhausted RAM and silently killed workers. **Use at most 8
parallel workers on 32 GB** for the Gompertz model. The free-intercept
BurstiVAR model is much lighter and tolerates more workers.

## 4. Comparing the two fits

### 4a. VAR parameters side by side

The substantive parameters — the dynamics — are the fairest comparison,
because both models estimate exactly the same four VAR coefficients and their
between-person SDs. Step 4a prints them side by side, with the ratio of
posterior SDs (BurstiVAR / GoBurstiVAR) as a simple efficiency measure:

```
            parameter  true GoBurstiVAR_mean GoBurstiVAR_PSD BurstiVAR_mean BurstiVAR_PSD PSD_ratio_Bi_over_Go
1    mu_phi11 (AR y1)  0.30           0.2639          0.0454         0.3056        0.0638                 1.41
2    mu_phi22 (AR y2)  0.20           0.1795          0.0452         0.1821        0.0608                 1.35
3 mu_phi21 (y1 -> y2) -0.15          -0.2088          0.0429        -0.2608        0.0515                 1.20
4 mu_phi12 (y2 -> y1) -0.10          -0.1177          0.0415        -0.1801        0.0458                 1.10
5         sigma_phi11  0.10           0.0986          0.0337         0.0958        0.0357                 1.06
6         sigma_phi22  0.10           0.1366          0.0495         0.1986        0.0604                 1.22
7         sigma_phi21  0.10           0.1272          0.0449         0.1611        0.0504                 1.12
8         sigma_phi12  0.10           0.0977          0.0307         0.0882        0.0296                 0.96
```

Two things to notice. First, the posterior means of the VAR fixed effects are
close between the two models and close to the truth: modeling the burst means
correctly, by either route, protects the dynamics. This mirrors the paper's
finding that accuracy and calibration of the VAR fixed effects were
substantively identical across the two specifications. Second, the posterior
SDs show where the models differ. In this run BurstiVAR's posterior SDs on
the four VAR fixed effects are 10 to 41 percent larger than GoBurstiVAR's
(the `PSD_ratio_Bi_over_Go` column): when its form is right, the parametric
model wrings more information out of the same 25 observations per person. In
the paper's 100-replication design that efficiency gain concentrated on the
small cross-lag `mu_phi12`, where tighter posteriors turned into a 14-point
power advantage (.98 vs .84 at T_b = 5, N = 100). A single small dataset
shows a noisy version of that picture, which is exactly why the paper needed
replications: treat the table above as an illustration of the comparison, not
as an estimate of the efficiency gain.

### 4b. Gompertz parameter recovery

GoBurstiVAR's own reward is that it returns interpretable growth parameters.
Step 4b compares the estimated population means and SDs of the six Gompertz
parameters against the generating values:

```
                 parameter true_mean est_mean    PSD true_SD est_SD
1    theta1 y1 (asymptote)      10.0   10.409 0.4020    2.00  2.362
2 theta2 y1 (displacement)       3.0    3.256 0.2280    0.50  0.701
3  theta3 y1 (growth rate)       1.0    1.020 0.0550    0.20  0.147
4    theta1 y2 (asymptote)       8.0    7.502 0.2755    1.50  1.479
5 theta2 y2 (displacement)       2.5    2.759 0.1974    0.40  0.540
6  theta3 y2 (growth rate)       0.8    0.912 0.0622    0.15  0.179
```

With only 5 bursts and 40 persons the population means of the Gompertz
parameters are recovered reasonably well, while their between-person SDs are
estimated with a lot of uncertainty — five points on a curve is not much
information to separate asymptote, displacement, and rate variability. Expect
wide posteriors on the Gompertz SDs whenever B is small.

### 4c. Burst-level means: two routes to the same curve

BurstiVAR estimates the burst means directly; GoBurstiVAR implies them
through the curve. Step 4c compares both against the truth. Note a subtlety:
because the Gompertz function is nonlinear in the person-specific parameters,
the true population burst mean is E[g(b; theta_i)] over persons, which is not
exactly the curve evaluated at the mean parameters. The script obtains the
true values by Monte Carlo, and the GoBurstiVAR column is a plug-in of the
posterior-mean parameters into the curve:

```
   burst variable true_MC BurstiVAR_free GoBurstiVAR_curve
1      1       y1    3.38           3.37              3.22
2      2       y1    6.55           6.69              6.81
3      3       y1    8.42           8.88              8.93
4      4       y1    9.29           9.79              9.85
5      5       y1    9.68          10.02             10.20
6      1       y2    2.65           2.54              2.48
7      2       y2    4.80           4.81              4.81
8      3       y2    6.28           6.08              6.27
9      4       y2    7.11           6.82              6.98
10     5       y2    7.54           7.16              7.29
```

Both models track the growth pattern. This is the sense in which BurstiVAR is
"form-free but not blind": its free intercepts recover whatever trajectory
the burst means follow, Gompertz included, without being told the form.

### 4d. Innovation covariance and convergence

```
             parameter true GoBurstiVAR_mean BurstiVAR_mean
1    innovation var y1  1.0            1.158          1.184
2 innovation cov y1,y2  0.3            0.333          0.286
3    innovation var y2  1.0            1.043          0.994
```

Finally, a quick convergence check on the reduced-MCMC run:

```
Max RHAT, GoBurstiVAR: 1.0572 | min ESS: 62
Max RHAT, BurstiVAR:   1.0285 | min ESS: 99
```

The RHAT values are acceptable for a demonstration, but the smallest
effective sample sizes are in the double digits — a direct consequence of the
short chains, and one more reason the reduced settings are for demonstration
only. For a real analysis, run the paper's full MCMC settings and check RHAT and
effective sample sizes parameter by parameter (the full `summarizePost`
tables, written to `tutorial_fit_goburstivar.csv` and
`tutorial_fit_burstivar.csv`, contain both).

## 5. When to prefer which model

The paper's Study 2 gives the decision its shape, and the shape is
asymmetric:

- **If the Gompertz form holds**, both models are correct, and GoBurstiVAR is
  more efficient. The gain matters most exactly where data are scarce: few
  time points per burst, modest N, and interest in small cross-lagged
  effects. At T_b = 5 and N = 100 the paper found power .98 vs .84 on
  `mu_phi12 = -0.10`; with more data (N = 500 or T_b = 20) both models
  reached power 1.00 and the restriction bought nothing detectable.
- **If the Gompertz form fails**, GoBurstiVAR is misspecified and its
  estimates inherit that error. BurstiVAR, which never assumed a form,
  remains correct.

The costs are therefore not symmetric: the penalty for using BurstiVAR when
GoBurstiVAR would have been right is some efficiency; the penalty for using
GoBurstiVAR when the form is wrong is bias. **BurstiVAR is the safe
default.** Reach for GoBurstiVAR when

1. theory or prior work predicts smooth, monotone, asymptoting growth across
   bursts (learning and skill-acquisition designs are the classic case);
2. the free-intercept estimates from a BurstiVAR fit — always a sensible
   first step — actually trace such a curve; and
3. you either need the efficiency (small T_b, modest N, small effects of
   interest) or you care substantively about the growth parameters
   themselves: an asymptote or a growth rate can be the research question,
   and BurstiVAR's free means will never give you one.

A practical workflow follows directly from this tutorial: fit BurstiVAR
first, plot its estimated burst means, and only then decide whether a
parametric growth model has earned its restrictions. If it has, fitting both
and checking that the dynamics agree — as we did here — is a cheap and
persuasive robustness check.

## Files and reproduction

- `tutorial.R` — runs everything above; start it from this folder.
- `mlGVARNoCorNoME_GompertzBurst.txt` — the GoBurstiVAR JAGS model
  (published Study 2 version).
- `mlGVARNoCorNoME_BurstIntercept_5Burst.txt` — the BurstiVAR free-intercept
  JAGS model for B = 5 (published Study 2 version).
- `posteriorSummaryStats.R` — posterior summary utilities (`summarizePost`).

The full Study 2 pipeline — 100 replications per cell, paper MCMC settings,
and the summary scripts behind the paper's Study 2 tables — lives in
`../../study2_goburstivar/`.
