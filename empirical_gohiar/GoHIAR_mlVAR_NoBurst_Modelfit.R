# ------------------------------------------------------------------------------
# Provenance header added for the public BurstiVAR repository.
# Original file: IP_1b/MOL/GoHIAR_mlVAR_NoBurst_Modelfit.R
# Modifications relative to the original:
#   1. This header was prepended.
#   2. Path configuration was made relative: data_path is now ".", so run
#      this script from this directory (empirical_gohiar/), where the JAGS
#      model files and posteriorSummaryStats.R also live.
#   3. The four load() calls now read from data_private/, a placeholder
#      directory for the private person-level inputs (datblock_Meaning.Rdata,
#      datblock_relationship.Rdata, datblock_accomplishment.Rdata,
#      datCov.Rdata), available from the original study team.
# All other content is identical to the version that produced the published
# results.
# ------------------------------------------------------------------------------

# ==============================================================================
# GoHIAR_mlVAR_NoBurst_Modelfit.R   (2026-07-17)
#
# 对照拟合: 无 burst 截距的标准 mlVAR (GoHIAR_mlVAR_NoBurst.txt), full data。
# 回应 advisor todo: burst-related changes 不建模时, VAR 推断会怎样?
#
# 数据构建与 GoHIAR_BurstiVAR_4Burst_BaselineDelta_Modelfit.R (Part 1) 完全一致:
#   同一 Y_full (N x 224 x 2), 同一 bivariate completeness 处理。
# MCMC 与全部经验拟合一致: chains=2, adapt=5000, burnin=5000, sample=20000, thin=1,
#   seeds 2026/2027 —— 保证与 BurstiVAR 拟合可比。
# ==============================================================================

rm(list = ls())

library(rjags)
load.module("dic")
library(coda)

data_path    <- "."   # run from this directory (empirical_gohiar/)
model_file   <- file.path(data_path, "GoHIAR_mlVAR_NoBurst.txt")
source_file  <- "./posteriorSummaryStats.R"

result_path <- file.path(data_path, "BurstiVAR_results_NoBurst")
if (!dir.exists(result_path)) dir.create(result_path, recursive = TRUE)

source(source_file)
setwd(data_path)
load("data_private/datblock_Meaning.Rdata")        # private input, available from the original study team
load("data_private/datblock_relationship.Rdata")   # private input, available from the original study team
load("data_private/datblock_accomplishment.Rdata") # private input, available from the original study team
load("data_private/datCov.Rdata")                  # private input, available from the original study team

cat("=== 数据加载完成 ===\n")

# --- 构建 full Y 数组 (与 Part 1 完全一致) ---
T4_full <- 224L

if (!exists("pInd")) pInd <- unique(datblock$PID)
if (!exists("N"))    N <- length(pInd)

Y_full <- array(NA, dim = c(N, T4_full, 2))

datblock$row_id <- ave(seq_len(nrow(datblock)), datblock$PID, FUN = seq_along)
datblock_relationship$row_id <- ave(seq_len(nrow(datblock_relationship)),
                                     datblock_relationship$PID, FUN = seq_along)

for (pp in 1:N) {
  pid <- pInd[pp]
  tmp_m <- datblock[datblock$PID == pid, ]
  tmp_m <- tmp_m[order(tmp_m$row_id), ]
  tmp_r <- datblock_relationship[datblock_relationship$PID == pid, ]
  tmp_r <- tmp_r[order(tmp_r$row_id), ]
  nT_person <- nrow(tmp_m)
  for (oo in 1:min(nT_person, T4_full)) {
    Y_full[pp, oo, 1] <- tmp_m$Meaning[oo]
    if (oo <= nrow(tmp_r)) {
      Y_full[pp, oo, 2] <- tmp_r$Relationship[oo]
    }
  }
}

n_partial <- 0
for (pp in 1:N) {
  for (tt in 1:T4_full) {
    has_y1 <- !is.na(Y_full[pp, tt, 1])
    has_y2 <- !is.na(Y_full[pp, tt, 2])
    if (has_y1 != has_y2) {
      Y_full[pp, tt, 1] <- NA
      Y_full[pp, tt, 2] <- NA
      n_partial <- n_partial + 1
    }
  }
}

cat("Full Y 数组维度:", dim(Y_full), "\n")
cat("部分观测 -> 全 NA 的 cell 数:", n_partial, "\n")
cat(sprintf("有效 bivariate 观测占比: %.1f%%\n\n", mean(!is.na(Y_full[,,1])) * 100))

# --- JAGS ---
W0 <- diag(rep(0.015, 2))
jags_data <- list(Y = Y_full, P = N, W0 = W0, df_res = 3L, T4 = T4_full)

parameters <- c("Level2Mean", "Level2Sigma", "sigma_innovation")

inits <- list(
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 2026),
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 2027)
)

cat("Step 1/3: jags.model (adapt=5000) ...\n")
t0 <- proc.time()
jagsModel <- jags.model(model_file, data = jags_data, inits = inits,
                        n.chains = 2, n.adapt = 5000, quiet = FALSE)
cat("  adapt 完成, 用时:", round((proc.time() - t0)[3] / 60, 1), "分钟\n")

cat("Step 2/3: burn-in (5000) ...\n")
t1 <- proc.time()
update(jagsModel, n.iter = 5000, progress.bar = "text")
cat("  burn-in 完成, 用时:", round((proc.time() - t1)[3] / 60, 1), "分钟\n")

cat("Step 3/3: sample (20000, thin=1) ...\n")
t2 <- proc.time()
codaSamples <- coda.samples(jagsModel, variable.names = parameters,
                            n.iter = 20000, thin = 1, progress.bar = "text")
cat("  sample 完成, 用时:", round((proc.time() - t2)[3] / 60, 1), "分钟\n")

save(codaSamples, file = file.path(result_path, "GoHIAR_mlVAR_NoBurst_codaSamples.RData"))
save(Y_full, file = file.path(result_path, "GoHIAR_mlVAR_NoBurst_Y.RData"))

cat("\n--- 收敛诊断 ---\n")
rhat_vals <- gelman.diag(codaSamples, multivariate = FALSE)$psrf[, "Point est."]
ess_vals  <- effectiveSize(codaSamples)
cat(sprintf("  max R-hat : %.4f\n", max(rhat_vals, na.rm = TRUE)))
cat(sprintf("  min ESS   : %.0f\n", min(ess_vals, na.rm = TRUE)))
if (max(rhat_vals, na.rm = TRUE) > 1.10) {
  bad <- names(rhat_vals[rhat_vals > 1.10])
  cat("  *** 警告: 未收敛参数:", paste(head(bad, 10), collapse = ", "), "\n")
}

resulttable <- summarizePost(codaSamples)
write.csv(resulttable, file.path(result_path, "GoHIAR_mlVAR_NoBurst_resulttable.csv"))
cat("\n结果表已保存\n")

samples_mat <- do.call(rbind, codaSamples)
var_labels <- c("Mean Meaning", "Mean Relationship",
                "AR Meaning", "AR Relationship",
                "CR Meaning->Relationship", "CR Relationship->Meaning")
cat("\n--- 关键参数 ---\n")
for (k in 1:6) {
  p <- paste0("Level2Mean[", k, "]")
  post <- samples_mat[, p]
  ll <- quantile(post, 0.025); ul <- quantile(post, 0.975)
  sig <- ifelse(ll > 0 | ul < 0, "***", "   ")
  cat(sprintf("  %-26s | Mean = %7.3f, 95%% CI = [%7.3f, %7.3f] %s\n",
              var_labels[k], mean(post), ll, ul, sig))
}
for (k in 1:6) {
  p <- paste0("Level2Sigma[", k, "]")
  post <- samples_mat[, p]
  cat(sprintf("  SD  %-22s | Mean = %7.3f, 95%% CI = [%7.3f, %7.3f]\n",
              var_labels[k], mean(post),
              quantile(post, 0.025), quantile(post, 0.975)))
}

cat("\n=== NoBurst 拟合完成 ===\n")
