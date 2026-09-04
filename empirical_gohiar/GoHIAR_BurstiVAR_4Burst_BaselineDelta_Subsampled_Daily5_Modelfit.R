# ------------------------------------------------------------------------------
# Provenance header added for the public BurstiVAR repository.
# Original file: IP_1b/MOL/GoHIAR_BurstiVAR_4Burst_BaselineDelta_Subsampled_Daily5_Modelfit.R
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
# GoHIAR_BurstiVAR_4Burst_BaselineDelta_Subsampled_Daily5_Modelfit.R  (2026-07-17)
#
# 按日采样敏感性 (advisor comment 13(4): "sample T_B daily over just 5 days?"):
#   每个 person x burst 取 5 个 occasion, 每天同一 6 小时时段, 连续 5 天
#   -> lag-1 间隔 = 24 小时 (stride = 4 slots), Tb = 5。
#
# 与 gap=1 (v2) 下采样的唯一区别: 窗口内 slot 间隔为 4 (同一时段跨日) 而非 1。
#   窗口 = {s, s+4, s+8, s+12, s+16}; s 自由滑动 (同时选择时段与起始日),
#   准则与 v2 相同: 双变量完整观测最多 -> 窗口中心靠 burst 中心 -> 最早;
#   完整观测 < 3 -> 该 person x burst 弃用。
#
# 模型/MCMC 与其他敏感性拟合完全一致 (BaselineDelta.txt; 2 chains,
#   adapt=5000, burnin=5000, sample=20000; seeds 2026/2027)。
# 注意: 本拟合结果暂不进正文, 先供检视 (预期: 24h-lag 下 AR/CR 接近零)。
# ==============================================================================

rm(list = ls())

library(rjags)
load.module("dic")
library(coda)

data_path    <- "."   # run from this directory (empirical_gohiar/)
model_file   <- file.path(data_path, "GoHIAR_BurstiVAR_4Burst_BaselineDelta.txt")
source_file  <- "./posteriorSummaryStats.R"

result_path <- file.path(data_path, "BurstiVAR_results_BaselineDelta_Subsampled")
if (!dir.exists(result_path)) dir.create(result_path, recursive = TRUE)

source(source_file)
setwd(data_path)
load("data_private/datblock_Meaning.Rdata")        # private input, available from the original study team
load("data_private/datblock_relationship.Rdata")   # private input, available from the original study team
load("data_private/datblock_accomplishment.Rdata") # private input, available from the original study team
load("data_private/datCov.Rdata")                  # private input, available from the original study team

# --- 构建 full Y 数组 (与 Part 1 完全一致) ---
T1_full <- 56L; T2_full <- 116L; T3_full <- 168L; T4_full <- 224L

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

burst_idx_full <- list(
  `1` = 1:T1_full,
  `2` = (T1_full + 1):T2_full,
  `3` = (T2_full + 1):T3_full,
  `4` = (T3_full + 1):T4_full
)

# --- 跨步窗口选择 (stride = 4, n_occ = 5) ---
best_window_start_daily <- function(obs_mask, n_occ = 5L, stride = 4L) {
  Lb   <- length(obs_mask)
  span <- (n_occ - 1L) * stride
  if (Lb < span + 1L) return(NA_integer_)
  n_windows <- Lb - span
  counts <- integer(n_windows)
  for (s in 1:n_windows) {
    counts[s] <- sum(obs_mask[seq(s, by = stride, length.out = n_occ)])
  }
  best_n <- max(counts)
  if (best_n < 3L) return(NA_integer_)
  cands <- which(counts == best_n)
  win_centers <- cands + span / 2
  burst_center <- (Lb + 1) / 2
  d <- abs(win_centers - burst_center)
  cands[which.min(d)][1]
}

n_occ <- 5L; stride <- 4L
T_total <- 4L * n_occ
Y_sub <- array(NA, dim = c(N, T_total, 2))
per_burst_n_used  <- matrix(0L, N, 4)
per_burst_dropped <- matrix(FALSE, N, 4)

for (pp in 1:N) {
  for (b in 1:4) {
    slots_b  <- burst_idx_full[[as.character(b)]]
    obs_mask <- !is.na(Y_full[pp, slots_b, 1]) & !is.na(Y_full[pp, slots_b, 2])
    ws <- best_window_start_daily(obs_mask, n_occ, stride)
    if (is.na(ws)) { per_burst_dropped[pp, b] <- TRUE; next }
    sel_rel <- seq(ws, by = stride, length.out = n_occ)
    sel_abs <- slots_b[sel_rel]
    sub_slots_b <- ((b - 1L) * n_occ + 1L):(b * n_occ)
    Y_sub[pp, sub_slots_b, 1] <- Y_full[pp, sel_abs, 1]
    Y_sub[pp, sub_slots_b, 2] <- Y_full[pp, sel_abs, 2]
    per_burst_n_used[pp, b] <- sum(obs_mask[sel_rel])
  }
}

cat("Y_sub 维度:", dim(Y_sub), "\n")
for (b in 1:4) {
  nu <- per_burst_n_used[, b]
  cat(sprintf("  Burst %d: mean obs = %.2f, range = [%d, %d], dropped = %d\n",
              b, mean(nu), min(nu), max(nu), sum(per_burst_dropped[, b])))
}
obs_cells <- sum(!is.na(Y_sub[,,1]))
cat(sprintf("有效 bivariate: %d / %d (%.1f%%)\n\n",
            obs_cells, prod(dim(Y_sub)[1:2]), obs_cells / prod(dim(Y_sub)[1:2]) * 100))

save(Y_sub, per_burst_n_used, per_burst_dropped,
     file = file.path(result_path, "GoHIAR_BurstiVAR_4Burst_BaselineDelta_Daily5_Y.RData"))

# --- JAGS (T1..T4 = 5,10,15,20; 与 Tb5 拟合同维) ---
W0 <- diag(rep(0.015, 2))
jags_data <- list(Y = Y_sub, P = N, W0 = W0, df_res = 3L,
                  T1 = 5L, T2 = 10L, T3 = 15L, T4 = 20L)

parameters <- c("mu_base", "delta", "Level2Mean", "Level2Sigma", "sigma_innovation")

inits <- list(
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 2026),
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 2027)
)

cat("Step 1/3: jags.model (adapt=5000) ...\n")
jagsModel <- jags.model(model_file, data = jags_data, inits = inits,
                        n.chains = 2, n.adapt = 5000, quiet = FALSE)

cat("Step 2/3: burn-in (5000) ...\n")
update(jagsModel, n.iter = 5000, progress.bar = "text")

cat("Step 3/3: sample (20000) ...\n")
codaSamples <- coda.samples(jagsModel, variable.names = parameters,
                            n.iter = 20000, thin = 1, progress.bar = "text")

save(codaSamples,
     file = file.path(result_path, "GoHIAR_BurstiVAR_4Burst_BaselineDelta_Daily5_codaSamples.RData"))

rhat_vals <- gelman.diag(codaSamples, multivariate = FALSE)$psrf[, "Point est."]
ess_vals  <- effectiveSize(codaSamples)
cat(sprintf("\nmax R-hat = %.4f, min ESS = %.0f\n",
            max(rhat_vals, na.rm = TRUE), min(ess_vals, na.rm = TRUE)))

resulttable <- summarizePost(codaSamples)
write.csv(resulttable,
          file.path(result_path, "GoHIAR_BurstiVAR_4Burst_BaselineDelta_Daily5_resulttable.csv"))

samples_mat <- do.call(rbind, codaSamples)
var_labels <- c("AR Meaning", "AR Relationship",
                "CR Meaning->Relationship", "CR Relationship->Meaning")
cat("\n--- VAR 参数 (24h lag) ---\n")
for (k in 9:12) {
  p <- paste0("Level2Mean[", k, "]")
  post <- samples_mat[, p]
  ll <- quantile(post, 0.025); ul <- quantile(post, 0.975)
  sig <- ifelse(ll > 0 | ul < 0, "***", "   ")
  cat(sprintf("  %-26s | Mean = %6.3f, 95%% CI = [%6.3f, %6.3f] %s\n",
              var_labels[k - 8], mean(post), ll, ul, sig))
}

cat("\n=== Daily5 拟合完成 ===\n")
