# ------------------------------------------------------------------------------
# Provenance header added for the public BurstiVAR repository.
# Original file: IP_1b/MOL/GoHIAR_BurstiVAR_4Burst_BaselineDelta_Subsampled_Modelfit.R
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
# GoHIAR_BurstiVAR_4Burst_BaselineDelta_Subsampled_Modelfit.R
#
# 目的: Part 2 of GoHIAR empirical example —— sparse within-burst sensitivity.
#       将每个 burst 内的 within-burst 数据缩减到 Tb_sub ∈ {3, 5},
#       看 BurstiVAR 能否在稀疏 within-burst 数据下仍然恢复 Part 1 的参数估计。
#
# *** 修正版 (v2): Person-specific sliding window, 严格 gap = 1 ***
#
#   对每个 person × burst:
#     1) 在该 burst 的 *原始 slot 索引上* 滑动一个长度为 Tb_sub 的窗口
#        (即窗口里的 Tb_sub 个 slot 在原始时间上严格连续 → gap=1, lag-1 = 6 小时)
#     2) 在所有可能窗口中, 挑选 *双变量观测齐全 slot 数最多* 的那个
#     3) Tie-break: 偏好窗口中心最靠近 burst 中心的; 仍平局取最早窗口
#     4) 若选中窗口的 observed 数 < 3 → 整个 person × burst 弃用 (全 NA)
#     5) 选中窗口直接搬运到 Y_sub 该 burst 段, 窗口内的 NA 也保留 (JAGS 填补)
#
#   v2 修正了 v1 的 bug: v1 在 *已过滤的 observed slot* 上取连续 Tb_sub 个,
#   导致原始时间上的间隔被拉大 (≠ gap 1), 使 AR 高估、CR 衰减.
#   v2 在 *原始 slot 索引* 上滑动窗口, 严格保 gap = 1, lag 与 Part 1 一致.
#
# Conditions:
#   Tb_sub = 3 → T1=3, T2=6,  T3=9,  T4=12
#   Tb_sub = 5 → T1=5, T2=10, T3=15, T4=20
#
# MCMC: 与 Part 1 完全一致 (chains=2, adapt=5000, burnin=5000, sample=20000, thin=1)
# JAGS 模型文件: 复用 GoHIAR_BurstiVAR_4Burst_BaselineDelta.txt 不变
# ==============================================================================

rm(list = ls())

library(rjags)
load.module("dic")
library(coda)

# ==============================================================================
# 1. 路径配置
# ==============================================================================

# ---- 本地 PC 路径 ----
data_path    <- "."   # run from this directory (empirical_gohiar/)
model_file   <- file.path(data_path, "GoHIAR_BurstiVAR_4Burst_BaselineDelta.txt")
source_file  <- "./posteriorSummaryStats.R"

# ---- Roar Cluster 路径 (取消注释以使用) ----
# data_path    <- "/storage/work/xjx5093/BurstIntercept/GoHIAR"
# model_file   <- file.path(data_path, "GoHIAR_BurstiVAR_4Burst_BaselineDelta.txt")
# source_file  <- "/storage/work/xjx5093/BurstIntercept/posteriorSummaryStats.R"

result_path <- file.path(data_path, "BurstiVAR_results_BaselineDelta_Subsampled")
if (!dir.exists(result_path)) dir.create(result_path, recursive = TRUE)

# ==============================================================================
# 2. 加载数据 & 辅助函数
# ==============================================================================

source(source_file)

setwd(data_path)
load("data_private/datblock_Meaning.Rdata")        # private input, available from the original study team
load("data_private/datblock_relationship.Rdata")   # private input, available from the original study team
load("data_private/datblock_accomplishment.Rdata") # private input, available from the original study team
load("data_private/datCov.Rdata")                  # private input, available from the original study team

cat("=== 数据加载完成 ===\n")
cat("  datblock 列名:", paste(colnames(datblock), collapse = ", "), "\n\n")

# ==============================================================================
# 3. 构建 full Y 数组 (与 Part 1 完全一致)
# ==============================================================================

cat("=== 构建 full Y 数组 (与 Part 1 一致) ===\n\n")

T1_full <- 56L
T2_full <- 116L
T3_full <- 168L
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

# Bivariate completeness
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
cat("部分观测 → 全 NA 的 cell 数:", n_partial, "\n")
cat(sprintf("Full 数据中有效 bivariate 观测占比: %.1f%%\n\n",
            mean(!is.na(Y_full[,,1])) * 100))

burst_idx_full <- list(
  `1` = 1:T1_full,
  `2` = (T1_full + 1):T2_full,
  `3` = (T2_full + 1):T3_full,
  `4` = (T3_full + 1):T4_full
)

# ==============================================================================
# 4. Person-specific sliding window (gap = 1)
# ==============================================================================

# 给定一个 burst 内的 obs_mask (该 burst 各 slot 是否双变量齐全, 长度 = burst 长度),
# 找出最优长度为 Tb_sub 的连续窗口的起始位置 (相对 burst 起点, 1-indexed).
# 窗口内的 Tb_sub 个 slot 在 *原始 slot 索引上严格相邻* (gap = 1).
# 准则:
#   1) 窗口内 observed 数最多
#   2) 平局 → 窗口中心最靠近 burst 中心
#   3) 仍平局 → 取最早窗口
# 若最优窗口的 observed 数 < 3, 返回 NA (该 person × burst 弃用).
best_window_start <- function(obs_mask, Tb_sub) {
  Lb <- length(obs_mask)
  if (Lb < Tb_sub) return(NA_integer_)

  n_windows <- Lb - Tb_sub + 1L
  csum   <- cumsum(c(0L, as.integer(obs_mask)))
  counts <- csum[(Tb_sub + 1L):(Lb + 1L)] - csum[1L:n_windows]

  best_n <- max(counts)
  if (best_n < 3L) return(NA_integer_)

  cands <- which(counts == best_n)
  win_centers <- cands + (Tb_sub - 1) / 2
  burst_center <- (Lb + 1) / 2
  d <- abs(win_centers - burst_center)
  cands[which.min(d)][1]
}

# 对给定 Tb_sub 构建 subsampled Y 数组
build_Y_sub <- function(Y_full, Tb_sub, burst_idx_full) {
  N <- dim(Y_full)[1]
  T_total <- 4L * Tb_sub
  Y_sub <- array(NA, dim = c(N, T_total, 2))

  per_burst_n_used    <- matrix(0L, N, 4)   # 选中窗口内的 observed 数
  per_burst_n_avail   <- matrix(0L, N, 4)   # burst 全段 observed 数 (full)
  per_burst_dropped   <- matrix(FALSE, N, 4)
  per_burst_win_start <- matrix(NA_integer_, N, 4)  # 窗口起点 (burst 内, 1-indexed)

  for (pp in 1:N) {
    for (b in 1:4) {
      slots_b  <- burst_idx_full[[as.character(b)]]
      obs_mask <- !is.na(Y_full[pp, slots_b, 1]) & !is.na(Y_full[pp, slots_b, 2])
      per_burst_n_avail[pp, b] <- sum(obs_mask)

      ws <- best_window_start(obs_mask, Tb_sub)
      if (is.na(ws)) {
        per_burst_dropped[pp, b] <- TRUE
        next
      }

      sel_rel <- ws:(ws + Tb_sub - 1L)        # burst 内位置 (连续)
      sel_abs <- slots_b[sel_rel]              # 绝对 slot 索引 (连续)

      sub_slots_b <- ((b - 1L) * Tb_sub + 1L):(b * Tb_sub)

      Y_sub[pp, sub_slots_b, 1] <- Y_full[pp, sel_abs, 1]
      Y_sub[pp, sub_slots_b, 2] <- Y_full[pp, sel_abs, 2]

      per_burst_n_used[pp, b]    <- sum(obs_mask[sel_rel])
      per_burst_win_start[pp, b] <- ws
    }
  }

  list(
    Y_sub               = Y_sub,
    per_burst_n_used    = per_burst_n_used,
    per_burst_n_avail   = per_burst_n_avail,
    per_burst_dropped   = per_burst_dropped,
    per_burst_win_start = per_burst_win_start,
    Tb_sub              = Tb_sub,
    T_total             = T_total
  )
}

# ==============================================================================
# 5. 单条件拟合 routine
# ==============================================================================

fit_one_condition <- function(Tb_sub) {

  cat("\n\n##############################################################\n")
  cat(sprintf("###  Fitting BurstiVAR with Tb_sub = %d (T_total = %d)  ###\n",
              Tb_sub, 4L * Tb_sub))
  cat("##############################################################\n\n")

  built <- build_Y_sub(Y_full, Tb_sub, burst_idx_full)
  Y_sub <- built$Y_sub

  total_cells <- prod(dim(Y_sub)[1:2])
  obs_cells   <- sum(!is.na(Y_sub[,,1]))
  cat(sprintf("Y_sub 维度: %s\n", paste(dim(Y_sub), collapse = " x ")))
  cat(sprintf("Y_sub 总 cell: %d, 有效 bivariate: %d (%.1f%%), 缺失: %d (%.1f%%)\n",
              total_cells, obs_cells, obs_cells / total_cells * 100,
              total_cells - obs_cells,
              (total_cells - obs_cells) / total_cells * 100))

  cat("\n每 burst 选中窗口内的 observed 数 (across persons):\n")
  for (b in 1:4) {
    nu <- built$per_burst_n_used[, b]
    cat(sprintf("  Burst %d: mean = %.2f, median = %d, range = [%d, %d], dropped persons = %d\n",
                b, mean(nu), median(nu), min(nu), max(nu),
                sum(built$per_burst_dropped[, b])))
  }

  cat("\n每 burst 窗口起始位置分布 (burst 内 1-indexed; 仅 non-dropped):\n")
  for (b in 1:4) {
    ws <- built$per_burst_win_start[, b]
    Lb <- length(burst_idx_full[[as.character(b)]])
    if (all(is.na(ws))) {
      cat(sprintf("  Burst %d (length %d): all dropped\n", b, Lb))
    } else {
      cat(sprintf("  Burst %d (length %d): mean start = %.1f, range = [%d, %d]\n",
                  b, Lb,
                  mean(ws, na.rm = TRUE),
                  min(ws, na.rm = TRUE),
                  max(ws, na.rm = TRUE)))
    }
  }

  pct_ge2 <- mean(built$per_burst_n_used >= 2) * 100
  pct_ge3 <- mean(built$per_burst_n_used >= 3) * 100
  cat(sprintf("\n%% of person x burst cells with >=2 valid obs: %.1f%%\n", pct_ge2))
  cat(sprintf("%% of person x burst cells with >=3 valid obs: %.1f%%\n", pct_ge3))

  tag <- sprintf("Tb%d", Tb_sub)
  Yarray_file      <- file.path(result_path,
    sprintf("GoHIAR_BurstiVAR_4Burst_BaselineDelta_%s_Y.RData", tag))
  resulttable_file <- file.path(result_path,
    sprintf("GoHIAR_BurstiVAR_4Burst_BaselineDelta_%s_resulttable.csv", tag))
  codasamples_file <- file.path(result_path,
    sprintf("GoHIAR_BurstiVAR_4Burst_BaselineDelta_%s_codaSamples.RData", tag))
  diff_file        <- file.path(result_path,
    sprintf("GoHIAR_BurstiVAR_4Burst_BaselineDelta_%s_SDdiff.csv", tag))

  save(Y_sub, built, file = Yarray_file)
  cat("\nY_sub 数组已保存至:", Yarray_file, "\n")

  T1 <- as.integer(Tb_sub * 1)
  T2 <- as.integer(Tb_sub * 2)
  T3 <- as.integer(Tb_sub * 3)
  T4 <- as.integer(Tb_sub * 4)

  W0_val <- 0.015
  W0 <- diag(rep(W0_val, 2))

  jags_data <- list(
    Y      = Y_sub,
    P      = N,
    W0     = W0,
    df_res = 3L,
    T1     = T1,
    T2     = T2,
    T3     = T3,
    T4     = T4
  )

  cat(sprintf("\nJAGS data: P=%d, T1=%d, T2=%d, T3=%d, T4=%d\n",
              N, T1, T2, T3, T4))

  parameters <- c("mu_base", "delta", "Level2Mean", "Level2Sigma", "sigma_innovation")

  inits <- list(
    list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 2026),
    list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 2027)
  )

  if (!file.exists(model_file)) stop("模型文件不存在: ", model_file)

  cat("\nStep 1/3: jags.model (adapt=5000) ...\n")
  t0 <- proc.time()
  jagsModel <- jags.model(
    file     = model_file,
    data     = jags_data,
    inits    = inits,
    n.chains = 2,
    n.adapt  = 5000,
    quiet    = FALSE
  )
  cat("  adapt 完成, 用时:", round((proc.time() - t0)[3] / 60, 1), "分钟\n")

  cat("\nStep 2/3: burn-in (5000) ...\n")
  t1 <- proc.time()
  update(jagsModel, n.iter = 5000, progress.bar = "text")
  cat("  burn-in 完成, 用时:", round((proc.time() - t1)[3] / 60, 1), "分钟\n")

  cat("\nStep 3/3: sample (20000, thin=1) ...\n")
  t2 <- proc.time()
  codaSamples <- coda.samples(
    model          = jagsModel,
    variable.names = parameters,
    n.iter         = 20000,
    thin           = 1,
    progress.bar   = "text"
  )
  cat("  sample 完成, 用时:", round((proc.time() - t2)[3] / 60, 1), "分钟\n")

  save(codaSamples, file = codasamples_file)
  cat("MCMC samples 已保存至:", codasamples_file, "\n")

  cat("\n--- 收敛诊断 ---\n")
  rhat_vals <- gelman.diag(codaSamples, multivariate = FALSE)$psrf[, "Point est."]
  ess_vals  <- effectiveSize(codaSamples)
  cat(sprintf("  max R-hat : %.4f\n", max(rhat_vals, na.rm = TRUE)))
  cat(sprintf("  mean R-hat: %.4f\n", mean(rhat_vals, na.rm = TRUE)))
  cat(sprintf("  # params with R-hat > 1.10: %d\n",
              sum(rhat_vals > 1.10, na.rm = TRUE)))
  cat(sprintf("  min ESS   : %.0f\n", min(ess_vals, na.rm = TRUE)))
  cat(sprintf("  mean ESS  : %.0f\n", mean(ess_vals, na.rm = TRUE)))

  if (max(rhat_vals, na.rm = TRUE) > 1.10) {
    bad <- names(rhat_vals[rhat_vals > 1.10])
    cat("  *** 警告: 未收敛参数:", paste(head(bad, 10), collapse = ", "), "\n")
  }

  resulttable <- summarizePost(codaSamples)
  write.csv(resulttable, resulttable_file)
  cat("\n结果表已保存至:", resulttable_file, "\n")

  samples_mat <- do.call(rbind, codaSamples)
  par_names   <- colnames(samples_mat)

  cat("\n--- 关键 VAR 参数 (Level2Mean[9:12]) ---\n")
  var_labels <- c("AR Meaning", "AR Relationship",
                  "CR Meaning->Relationship", "CR Relationship->Meaning")
  for (k in 9:12) {
    p <- paste0("Level2Mean[", k, "]")
    if (p %in% par_names) {
      post <- samples_mat[, p]
      ll <- quantile(post, 0.025); ul <- quantile(post, 0.975)
      sig <- ifelse(ll > 0 | ul < 0, "***", "   ")
      cat(sprintf("  %-30s | Mean = %6.3f, 95%% CI = [%6.3f, %6.3f] %s\n",
                  var_labels[k - 8], mean(post), ll, ul, sig))
    }
  }

  cat("\n--- Baseline & Delta ---\n")
  for (j in 1:2) {
    p <- paste0("mu_base[", j, "]")
    if (p %in% par_names) {
      post <- samples_mat[, p]
      cat(sprintf("  mu_base[%d] | Mean = %6.2f, 95%% CI = [%6.2f, %6.2f]\n",
                  j, mean(post), quantile(post, 0.025), quantile(post, 0.975)))
    }
  }
  for (d in 1:6) {
    p <- paste0("delta[", d, "]")
    if (p %in% par_names) {
      post <- samples_mat[, p]
      ll <- quantile(post, 0.025); ul <- quantile(post, 0.975)
      sig <- ifelse(ll > 0 | ul < 0, "***", "   ")
      cat(sprintf("  delta[%d]   | Mean = %6.2f, 95%% CI = [%6.2f, %6.2f] %s\n",
                  d, mean(post), ll, ul, sig))
    }
  }

  sd_comparisons <- list(
    list(k = 3, base = 1, label = "Burst2_y1 - Burst1_y1 (Meaning)"),
    list(k = 4, base = 2, label = "Burst2_y2 - Burst1_y2 (Relationship)"),
    list(k = 5, base = 1, label = "Burst3_y1 - Burst1_y1 (Meaning)"),
    list(k = 6, base = 2, label = "Burst3_y2 - Burst1_y2 (Relationship)"),
    list(k = 7, base = 1, label = "Burst4_y1 - Burst1_y1 (Meaning)"),
    list(k = 8, base = 2, label = "Burst4_y2 - Burst1_y2 (Relationship)")
  )
  sd_diff_results <- data.frame()
  for (comp in sd_comparisons) {
    p_k    <- paste0("Level2Sigma[", comp$k, "]")
    p_base <- paste0("Level2Sigma[", comp$base, "]")
    if (p_k %in% par_names && p_base %in% par_names) {
      diff_post <- samples_mat[, p_k] - samples_mat[, p_base]
      ll <- quantile(diff_post, 0.025); ul <- quantile(diff_post, 0.975)
      sig <- ifelse(ll > 0 | ul < 0, "***", "   ")
      sd_diff_results <- rbind(sd_diff_results, data.frame(
        Label     = comp$label,
        Mean_diff = round(mean(diff_post), 4),
        SD_diff   = round(sd(diff_post), 4),
        LL_95     = round(ll, 4),
        UL_95     = round(ul, 4),
        Sig       = sig,
        stringsAsFactors = FALSE
      ))
    }
  }
  write.csv(sd_diff_results, diff_file, row.names = FALSE)
  cat("\nSD 差值表已保存至:", diff_file, "\n")

  invisible(list(
    Tb_sub      = Tb_sub,
    rhat_max    = max(rhat_vals, na.rm = TRUE),
    ess_min     = min(ess_vals, na.rm = TRUE),
    resulttable = resulttable_file,
    coda        = codasamples_file
  ))
}

# ==============================================================================
# 6. 跑两个条件
# ==============================================================================

t_global <- proc.time()

res_Tb3 <- fit_one_condition(Tb_sub = 3L)
res_Tb5 <- fit_one_condition(Tb_sub = 5L)

cat("\n\n=== 全部完成 ===\n")
cat(sprintf("Tb=3: max R-hat = %.4f, min ESS = %.0f\n",
            res_Tb3$rhat_max, res_Tb3$ess_min))
cat(sprintf("Tb=5: max R-hat = %.4f, min ESS = %.0f\n",
            res_Tb5$rhat_max, res_Tb5$ess_min))
cat(sprintf("总用时: %.1f 分钟\n",
            (proc.time() - t_global)[3] / 60))
