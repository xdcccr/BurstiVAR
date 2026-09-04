# ------------------------------------------------------------------------------
# tutorial.R -- BurstiVAR tutorial: simulate, fit, and read a small burst dataset
#
# New file written for the public repository. The data-generating code follows
# the Study 1 pipeline script
#   paper_simulation_archive/Study1_BurstiVAR/code/common_scripts_FlexibleBurst/
#     mlGVARNoCorNoME_BurstIntercept_Flexible_DataGeneration.R
# and the fitting code follows
#     mlGVARNoCorNoME_BurstIntercept_Flexible_Modelfit.R
# with a smaller sample (N = 50), one replication, and REDUCED MCMC settings
# suitable for a tutorial run (see Section 4 below). The JAGS model file and
# posteriorSummaryStats.R in this folder are the actual files used for the
# published Study 1 results.
#
# Run from this directory:
#   Rscript tutorial.R
# Requires: R (tested with 4.3.2), JAGS (tested with 4.3.1), and the R
# packages rjags, coda, mvtnorm.
# Expected runtime: about a minute on a typical desktop machine.
# ------------------------------------------------------------------------------

library(mvtnorm)
library(rjags)

# ==============================================================================
# Section 1. Design: a small measurement-burst study
# ==============================================================================
# 3 bursts of 5 occasions each, 50 persons, 2 variables (call them y1 and y2).
# These are the paper's Study 1 generating values, at the recommended burst
# length Tb = 5 and a deliberately small N so the tutorial runs quickly.

nBurst       <- 3    # number of bursts (B)
nT_per_burst <- 5    # occasions per burst (Tb)
nT           <- nBurst * nT_per_burst   # total occasions per person (15)
nP           <- 50   # persons (N)
ydim         <- 2    # bivariate outcome

# --- Burst-intercept means and SDs (paper's Study 1 design) ------------------
# Burst b: mean  = ( (b-1)*1.0 ,  (b-1)*0.5 + 1.0 )
#          SD    = ( 1.0 + (b-1)*0.5 , 1.0 + (b-1)*0.5 )
# So: burst 1 means (0.0, 1.0) SDs (1.0, 1.0)
#     burst 2 means (1.0, 1.5) SDs (1.5, 1.5)
#     burst 3 means (2.0, 2.0) SDs (2.0, 2.0)
Intercept_mu_list    <- list()
Intercept_sigma_list <- list()
for (b in 1:nBurst) {
  Intercept_mu_list[[b]]    <- c((b - 1) * 1.0, (b - 1) * 0.5 + 1.0)
  Intercept_sigma_list[[b]] <- c(1.0 + (b - 1) * 0.5, 1.0 + (b - 1) * 0.5)
}

# --- Within-burst dynamics (paper's Study 1 design) --------------------------
AR_mean1  <-  0.3    # autoregression of y1
AR_mean2  <-  0.2    # autoregression of y2
CR12_mean <- -0.15   # cross-lagged effect y1 -> y2
CR21_mean <-  0      # cross-lagged effect y2 -> y1 (zero in the design)

sigma_AR <- c(0.1, 0.1)   # between-person SDs of the AR coefficients
sigma_CR <- c(0.1, 0.1)   # between-person SDs of the cross-lagged coefficients

sigma_innovation  <- 1.0  # innovation SD (variance 1)
error_correlation <- 0.3  # innovation correlation (covariance 0.3)

# --- Assemble the person-parameter distribution ------------------------------
# Parameter order (this order matches the JAGS model file):
#   [1:2]  burst-1 intercepts (y1, y2)
#   [3:4]  burst-2 intercepts (y1, y2)
#   [5:6]  burst-3 intercepts (y1, y2)
#   [7]    AR y1, [8] AR y2, [9] CR y1->y2, [10] CR y2->y1
n_intercept_params <- nBurst * ydim
n_total_params     <- n_intercept_params + 4

param_means <- c(unlist(Intercept_mu_list),
                 AR_mean1, AR_mean2, CR12_mean, CR21_mean)
param_sds   <- c(unlist(Intercept_sigma_list),
                 sigma_AR[1], sigma_AR[2], sigma_CR[1], sigma_CR[2])

# Person parameters are drawn independently (identity correlation matrix),
# as in Study 1.
Sigma_full <- diag(param_sds) %*% diag(n_total_params) %*% diag(param_sds)

Sigma_innovation <- matrix(
  c(sigma_innovation^2,
    error_correlation * sigma_innovation^2,
    error_correlation * sigma_innovation^2,
    sigma_innovation^2),
  2, 2)

# ==============================================================================
# Section 2. Simulate one burst dataset
# ==============================================================================
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

    # First occasion of the burst: trend + a fresh innovation.
    # The VAR process restarts within each burst; nothing carries across the
    # burst gap except the person's parameters.
    Y_trend[i, T_start, 1:ydim] <- person_params[i, intercept_idx]
    Y_dev[i, T_start, ]   <- rmvnorm(1, mean = c(0, 0), sigma = Sigma_innovation)
    Y[i, T_start, ]       <- Y_dev[i, T_start, ] + Y_trend[i, T_start, ]

    # Remaining occasions: VAR(1) on the deviations from the burst intercept.
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

# Long format, saved for inspection (id, time, y1, y2).
dat_long <- data.frame(
  id   = rep(1:nP, each = nT),
  time = rep(0:(nT - 1), times = nP),
  y1   = as.vector(t(Y[, , 1])),
  y2   = as.vector(t(Y[, , 2]))
)
write.csv(dat_long, "tutorial_simulated_data.csv", row.names = FALSE)
cat("Simulated", nP, "persons x", nT, "occasions;",
    "saved to tutorial_simulated_data.csv\n")

# ==============================================================================
# Section 3. Prepare the data for JAGS
# ==============================================================================
# The model file mlGVARNoCorNoME_BurstIntercept_3Burst.txt is the actual file
# used for the paper's Study 1 three-burst conditions. It expects:
#   Y      : persons x occasions x 2 array
#   P      : number of persons
#   T1, T2, T3 : cumulative burst boundaries (here 5, 10, 15)
#   W0, df_res : Wishart prior for the innovation precision matrix

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

# Chain-specific RNG seeds, following the paper's (r, r + 500) convention
# with r = 1 for this single tutorial replication.
inits <- list(
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 1),
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 501)
)

# ==============================================================================
# Section 4. Fit the model (REDUCED MCMC settings)
# ==============================================================================
# The paper used 2 chains with 5000 adaptation, 5000 burn-in, and 20000
# sampling iterations. Here we use 2 chains with 1000 adaptation, 1000
# burn-in, and 2000 sampling iterations so the fit finishes in seconds.
# This is enough to see the model work on this small, well-behaved
# dataset, but it is NOT enough for reporting results: expect smaller
# effective sample sizes and less stable interval endpoints than with the
# full settings. For real analyses use the paper's settings and check
# convergence (R-hat close to 1, adequate ESS) before interpreting anything.

t0 <- proc.time()
jagsModel <- jags.model("mlGVARNoCorNoME_BurstIntercept_3Burst.txt",
                        data = jags_data, inits = inits,
                        n.chains = 2, n.adapt = 1000)
update(jagsModel, n.iter = 1000)                       # burn-in
codaSamples <- coda.samples(jagsModel,
                            variable.names = parameters,
                            n.iter = 2000)
cat("MCMC finished in", round((proc.time() - t0)[3], 1), "seconds\n")

# ==============================================================================
# Section 5. Summarize the posterior
# ==============================================================================
source("posteriorSummaryStats.R")   # the same summary utility used in the paper

resulttable <- summarizePost(codaSamples)
write.csv(resulttable, "tutorial_result_summary.csv")

# Print the summary next to the generating values.
true_values <- c(
  # Level2Mean[1..10]
  param_means,
  # Level2Sigma[1..10]
  param_sds,
  # sigma_innovation[1,1], [2,1], [1,2], [2,2]
  Sigma_innovation[1, 1], Sigma_innovation[2, 1],
  Sigma_innovation[1, 2], Sigma_innovation[2, 2]
)
comparison <- cbind(true = true_values,
                    resulttable[, c("mean", "PSD",
                                    "PCI 2.50%", "PCI 97.50%",
                                    "ESS", "RHAT")])
cat("\nPosterior summary vs generating values\n")
print(round(comparison, 3))

cat("\nFull summary table written to tutorial_result_summary.csv\n")
