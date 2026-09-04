# ==============================================================================
# tutorial.R -- GoBurstiVAR tutorial (companion script to tutorial.md)
#
# New script written for this repository. The data-generation block reuses the
# generating logic of the Study 2 simulation
# (paper_simulation_archive/Study2_Gompertz/code/
#  mlGVARNoCorNoME_GompertzBurst_DataGeneration.R), scaled down to one small
# dataset. The two JAGS model files and posteriorSummaryStats.R in this folder
# are identical (apart from provenance headers) to the versions that produced
# the published Study 2 results.
#
# Run from this directory:
#   Rscript tutorial.R
#
# MCMC here is deliberately reduced for a tutorial-sized run (2 chains,
# 1000 adaptation + 1000 burn-in + 2000 sampling iterations). The paper used
# 2 chains, 5000 adaptation + 5000 burn-in + 20000 sampling iterations.
# ==============================================================================

base_path <- "."   # run from this directory

options(width = 200)

library(mvtnorm)
library(rjags)
library(coda)

source(file.path(base_path, "posteriorSummaryStats.R"))

## =============================================================================
## Step 1. Simulate one Gompertz-burst dataset (B = 5, T_b = 5, N = 40)
## =============================================================================
## Burst-level trend: trend_b = theta1 * exp(-theta2 * exp(-b * theta3)),
## b = 1,...,5. Trend is constant within a burst; a VAR(1) process runs on the
## deviations within each burst and resets at every burst start.

set.seed(1)

nP          <- 40    # persons (small, so the tutorial finishes in minutes)
nBurst      <- 5     # bursts (B)
T_per_burst <- 5     # time points per burst (T_b)
nT          <- nBurst * T_per_burst   # 25 total time points

# Gompertz parameter means (asymptote, displacement, growth rate)
theta10_1 <- 10;  theta20_1 <- 3;    theta30_1 <- 1.0   # y1
theta10_2 <- 8;   theta20_2 <- 2.5;  theta30_2 <- 0.8   # y2

# Between-person SDs of the Gompertz parameters
sigma1_1 <- 2;    sigma2_1 <- 0.5;   sigma3_1 <- 0.2    # y1
sigma1_2 <- 1.5;  sigma2_2 <- 0.4;   sigma3_2 <- 0.15   # y2

# VAR parameter means (same values as the paper's Study 2)
AR_mean1  <-  0.3    # phi_11: y1 autoregression
AR_mean2  <-  0.2    # phi_22: y2 autoregression
CR12_mean <- -0.15   # phi_21: y1 -> y2 cross-lag
CR21_mean <- -0.10   # phi_12: y2 -> y1 cross-lag (the small effect)

sigma_AR <- 0.1      # SD of AR coefficients
sigma_CR <- 0.1      # SD of cross-lag coefficients

# Innovations (fixed across persons)
sigma_innovation  <- 1.0
error_correlation <- 0.3

param_means <- c(theta10_1, theta20_1, theta30_1,
                 theta10_2, theta20_2, theta30_2,
                 AR_mean1, AR_mean2, CR12_mean, CR21_mean)
param_sds   <- c(sigma1_1, sigma2_1, sigma3_1,
                 sigma1_2, sigma2_2, sigma3_2,
                 sigma_AR, sigma_AR, sigma_CR, sigma_CR)

# Independent random effects (as in the paper)
Sigma_full <- diag(param_sds) %*% diag(10) %*% diag(param_sds)

Sigma_innov <- matrix(c(sigma_innovation^2,
                        error_correlation * sigma_innovation^2,
                        error_correlation * sigma_innovation^2,
                        sigma_innovation^2), 2, 2)

gompertz <- function(b, theta1, theta2, theta3) {
  theta1 * exp(-theta2 * exp(-b * theta3))
}

burst_starts <- seq(1, nT, by = T_per_burst)
burst_ends   <- seq(T_per_burst, nT, by = T_per_burst)

Y            <- array(NA, c(nP, nT, 2))
Y_trend      <- array(NA, c(nP, nT, 2))
Y_dev        <- array(NA, c(nP, nT, 2))
burst_trends <- array(NA, c(nP, nBurst, 2))

person_params <- matrix(NA, nP, 10)
for (i in 1:nP) {
  person_params[i, ] <- rmvnorm(1, param_means, Sigma_full)
}

theta1_y1 <- person_params[, 1]; theta2_y1 <- person_params[, 2]; theta3_y1 <- person_params[, 3]
theta1_y2 <- person_params[, 4]; theta2_y2 <- person_params[, 5]; theta3_y2 <- person_params[, 6]
AR1  <- person_params[, 7];  AR2  <- person_params[, 8]
CR12 <- person_params[, 9];  CR21 <- person_params[, 10]

for (i in 1:nP) {
  for (b in 1:nBurst) {
    burst_trends[i, b, 1] <- gompertz(b, theta1_y1[i], theta2_y1[i], theta3_y1[i])
    burst_trends[i, b, 2] <- gompertz(b, theta1_y2[i], theta2_y2[i], theta3_y2[i])
  }
  for (b in 1:nBurst) {
    t_start <- burst_starts[b]
    t_end   <- burst_ends[b]

    # First time point of the burst: deviation process resets
    Y_trend[i, t_start, 1:2] <- burst_trends[i, b, ]
    Y_dev[i, t_start, ]      <- rmvnorm(1, mean = c(0, 0), sigma = Sigma_innov)
    Y[i, t_start, ]          <- Y_trend[i, t_start, ] + Y_dev[i, t_start, ]

    # Remaining time points: VAR(1) on the deviations
    for (t in (t_start + 1):t_end) {
      Y_trend[i, t, 1:2] <- burst_trends[i, b, ]
      mean_dev <- c(AR1[i] * Y_dev[i, t - 1, 1] + CR21[i] * Y_dev[i, t - 1, 2],
                    AR2[i] * Y_dev[i, t - 1, 2] + CR12[i] * Y_dev[i, t - 1, 1])
      Y_dev[i, t, ] <- rmvnorm(1, mean = mean_dev, sigma = Sigma_innov)
      Y[i, t, ]     <- Y_trend[i, t, ] + Y_dev[i, t, ]
    }
  }
}

# Long-format copy, in case you want to inspect or plot the data
dat_long <- data.frame(
  id    = rep(1:nP, each = nT),
  time  = rep(seq(0, 9.9, length.out = nT), nP),
  burst = rep(rep(1:nBurst, each = T_per_burst), nP),
  y1    = as.vector(t(Y[, , 1])),
  y2    = as.vector(t(Y[, , 2]))
)
write.csv(dat_long, file.path(base_path, "tutorial_data_gompertzburst.csv"),
          row.names = FALSE)

cat("Simulated dataset:", nP, "persons x", nT, "time points (",
    nBurst, "bursts x", T_per_burst, "points ).\n")
cat("Population burst-level means implied by the mean Gompertz parameters:\n")
print(round(rbind(y1 = gompertz(1:5, theta10_1, theta20_1, theta30_1),
                  y2 = gompertz(1:5, theta10_2, theta20_2, theta30_2)), 3))

## =============================================================================
## Step 2. Fit the SAME dataset with the true model (GoBurstiVAR)
## =============================================================================

jags_data <- list(
  Y = Y, P = nP, W0 = diag(2), df_res = 3,
  T1 = T_per_burst,        # end of burst 1
  T2 = T_per_burst * 2,    # end of burst 2
  T3 = T_per_burst * 3,    # end of burst 3
  T4 = T_per_burst * 4,    # end of burst 4
  T5 = nT                  # end of burst 5
)

parameters <- c("Level2Mean", "Level2Sigma", "sigma_innovation")

inits <- list(
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 1),
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 501)
)

# Reduced tutorial MCMC (paper: 5000 / 5000 / 20000)
n_adapt  <- 1000
n_burnin <- 1000
n_iter   <- 2000

cat("\n===== Fitting GoBurstiVAR (Gompertz burst means; the true model) =====\n")
t0 <- proc.time()
gompModel <- jags.model(file.path(base_path, "mlGVARNoCorNoME_GompertzBurst.txt"),
                        data = jags_data, inits = inits,
                        n.chains = 2, n.adapt = n_adapt, quiet = TRUE)
update(gompModel, n.iter = n_burnin, progress.bar = "none")
gompSamples <- coda.samples(gompModel, variable.names = parameters,
                            n.iter = n_iter, progress.bar = "none")
t_gomp <- proc.time() - t0
cat("GoBurstiVAR fit took", round(t_gomp[3] / 60, 1), "minutes.\n")

fit_gomp <- summarizePost(gompSamples)

## =============================================================================
## Step 3. Fit the SAME dataset with BurstiVAR (free burst intercepts)
## =============================================================================

cat("\n===== Fitting BurstiVAR (free burst intercepts) =====\n")
t0 <- proc.time()
biModel <- jags.model(file.path(base_path, "mlGVARNoCorNoME_BurstIntercept_5Burst.txt"),
                      data = jags_data, inits = inits,
                      n.chains = 2, n.adapt = n_adapt, quiet = TRUE)
update(biModel, n.iter = n_burnin, progress.bar = "none")
biSamples <- coda.samples(biModel, variable.names = parameters,
                          n.iter = n_iter, progress.bar = "none")
t_bi <- proc.time() - t0
cat("BurstiVAR fit took", round(t_bi[3] / 60, 1), "minutes.\n")

fit_bi <- summarizePost(biSamples)

## =============================================================================
## Step 4. Compare the two fits
## =============================================================================
## Parameter index maps (see tutorial.md):
##   GoBurstiVAR:  Level2Mean[1:6] Gompertz (y1: asy, dis, rate; y2: asy, dis, rate)
##                 Level2Mean[7:10] = phi_11, phi_22, phi_21 (y1->y2), phi_12 (y2->y1)
##   BurstiVAR:    Level2Mean[1:10] burst intercepts (b1 y1, b1 y2, ..., b5 y1, b5 y2)
##                 Level2Mean[11:14] = phi_11, phi_22, phi_21 (y1->y2), phi_12 (y2->y1)

## ---- 4a. VAR fixed effects and random-effect SDs, side by side --------------
var_labels <- c("mu_phi11 (AR y1)", "mu_phi22 (AR y2)",
                "mu_phi21 (y1 -> y2)", "mu_phi12 (y2 -> y1)",
                "sigma_phi11", "sigma_phi22", "sigma_phi21", "sigma_phi12")
var_true   <- c(0.30, 0.20, -0.15, -0.10, rep(0.10, 4))

gomp_rows <- c(paste0("Level2Mean[", 7:10, "]"),  paste0("Level2Sigma[", 7:10, "]"))
bi_rows   <- c(paste0("Level2Mean[", 11:14, "]"), paste0("Level2Sigma[", 11:14, "]"))

var_compare <- data.frame(
  parameter        = var_labels,
  true             = var_true,
  GoBurstiVAR_mean = fit_gomp[gomp_rows, "mean"],
  GoBurstiVAR_PSD  = fit_gomp[gomp_rows, "PSD"],
  BurstiVAR_mean   = fit_bi[bi_rows, "mean"],
  BurstiVAR_PSD    = fit_bi[bi_rows, "PSD"],
  row.names        = NULL
)
var_compare$PSD_ratio_Bi_over_Go <-
  round(var_compare$BurstiVAR_PSD / var_compare$GoBurstiVAR_PSD, 2)

cat("\n===== VAR parameters: GoBurstiVAR vs BurstiVAR (same dataset) =====\n")
print(var_compare, digits = 3)

## ---- 4b. Gompertz parameter recovery (GoBurstiVAR only) ---------------------
gomp_labels <- c("theta1 y1 (asymptote)", "theta2 y1 (displacement)", "theta3 y1 (growth rate)",
                 "theta1 y2 (asymptote)", "theta2 y2 (displacement)", "theta3 y2 (growth rate)")
gomp_recovery <- data.frame(
  parameter = gomp_labels,
  true_mean = c(10, 3, 1.0, 8, 2.5, 0.8),
  est_mean  = fit_gomp[paste0("Level2Mean[", 1:6, "]"), "mean"],
  PSD       = fit_gomp[paste0("Level2Mean[", 1:6, "]"), "PSD"],
  true_SD   = c(2, 0.5, 0.2, 1.5, 0.4, 0.15),
  est_SD    = fit_gomp[paste0("Level2Sigma[", 1:6, "]"), "mean"],
  row.names = NULL
)
cat("\n===== Gompertz parameter recovery (GoBurstiVAR) =====\n")
print(gomp_recovery, digits = 3)

## ---- 4c. Burst-level means: free intercepts vs Gompertz curve ---------------
## The population burst mean is E[ g(b; theta_i) ] over the person-level
## parameter distribution. Because g is nonlinear, this is not exactly the
## curve at the mean parameters, so we get the true values by Monte Carlo.
set.seed(2)
M <- 200000
mc_y1 <- sapply(1:5, function(b)
  mean(gompertz(b, rnorm(M, 10, 2),   rnorm(M, 3, 0.5),   rnorm(M, 1.0, 0.2))))
set.seed(3)
mc_y2 <- sapply(1:5, function(b)
  mean(gompertz(b, rnorm(M, 8, 1.5),  rnorm(M, 2.5, 0.4), rnorm(M, 0.8, 0.15))))

# BurstiVAR: free intercept means, ordered (b1 y1, b1 y2, ..., b5 y1, b5 y2)
bi_mu_y1 <- fit_bi[paste0("Level2Mean[", seq(1, 9, by = 2), "]"), "mean"]
bi_mu_y2 <- fit_bi[paste0("Level2Mean[", seq(2, 10, by = 2), "]"), "mean"]

# GoBurstiVAR: curve implied by the posterior-mean Gompertz parameters (plug-in)
go_mu_y1 <- gompertz(1:5, fit_gomp["Level2Mean[1]", "mean"],
                          fit_gomp["Level2Mean[2]", "mean"],
                          fit_gomp["Level2Mean[3]", "mean"])
go_mu_y2 <- gompertz(1:5, fit_gomp["Level2Mean[4]", "mean"],
                          fit_gomp["Level2Mean[5]", "mean"],
                          fit_gomp["Level2Mean[6]", "mean"])

burst_means <- data.frame(
  burst            = rep(1:5, 2),
  variable         = rep(c("y1", "y2"), each = 5),
  true_MC          = c(mc_y1, mc_y2),
  BurstiVAR_free   = c(bi_mu_y1, bi_mu_y2),
  GoBurstiVAR_curve = c(go_mu_y1, go_mu_y2)
)
cat("\n===== Burst-level means: truth vs the two models =====\n")
print(burst_means, digits = 3)

## ---- 4d. Innovation covariance ----------------------------------------------
innov_compare <- data.frame(
  parameter        = c("innovation var y1", "innovation cov y1,y2", "innovation var y2"),
  true             = c(1.0, 0.3, 1.0),
  GoBurstiVAR_mean = fit_gomp[c("sigma_innovation[1,1]", "sigma_innovation[1,2]",
                                "sigma_innovation[2,2]"), "mean"],
  BurstiVAR_mean   = fit_bi[c("sigma_innovation[1,1]", "sigma_innovation[1,2]",
                              "sigma_innovation[2,2]"), "mean"],
  row.names        = NULL
)
cat("\n===== Innovation covariance matrix =====\n")
print(innov_compare, digits = 3)

## ---- 4e. Convergence check --------------------------------------------------
cat("\nMax RHAT, GoBurstiVAR:", max(fit_gomp$RHAT),
    "| min ESS:", min(fit_gomp$ESS), "\n")
cat("Max RHAT, BurstiVAR:  ", max(fit_bi$RHAT),
    "| min ESS:", min(fit_bi$ESS), "\n")

## ---- Save the full summary tables -------------------------------------------
write.csv(fit_gomp, file.path(base_path, "tutorial_fit_goburstivar.csv"))
write.csv(fit_bi,   file.path(base_path, "tutorial_fit_burstivar.csv"))
write.csv(var_compare, file.path(base_path, "tutorial_var_comparison.csv"),
          row.names = FALSE)

cat("\nDone. Summary tables written to tutorial_fit_goburstivar.csv,",
    "tutorial_fit_burstivar.csv, tutorial_var_comparison.csv\n")
