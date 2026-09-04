# ------------------------------------------------------------------------------
# Provenance header added for the public BurstiVAR repository.
# Original file: IP_1b/MOL/GoHIAR_BurstiVAR_4Burst_GroupCovFull_Modelfit.R
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
# GoHIAR_BurstiVAR_4Burst_GroupCovFull_Modelfit.R
#
# 功能: 用 4-Burst BurstiVAR + *全 level-2 参数 group covariate* 版本拟合
#       Go-HIAR 经验数据 (full data, T = 224 slots).
#
#       相对 GoHIAR_BurstiVAR_4Burst_GroupCov_Modelfit.R 的变化:
#         (1) 模型文件切换: GoHIAR_BurstiVAR_4Burst_GroupCovFull.txt
#         (2) Group 现在也作为 AR (2) 和 CR (2) 的 covariate
#             → 新参数: mu_AR_base, mu_CR_base, delta_grp_AR, delta_grp_CR
#             → 派生: AR_ctrl, AR_intv, CR_ctrl, CR_intv
#         (3) Level2Mean[9:12] 不再监测 (被 mu_AR_base / mu_CR_base 替代)
#         (4) AR/CR 的 between-person 方差 tau_L2[9:12] 仍跨组共享
#         (5) delta_grp_AR / delta_grp_CR 的 prior SD = 0.3 (信息化)
#
#       y1 = Meaning of life (MOL), 0-100 scale
#       y2 = Relationship,          0-100 scale
#
# 默认组筛选: c(1, 3) = G1 Control + G3 PPI+Med
#
# MCMC: chains=2, adapt=5000, burnin=5000, iter=20000, thin=1
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
model_file   <- file.path(data_path, "GoHIAR_BurstiVAR_4Burst_GroupCovFull.txt")
source_file  <- "./posteriorSummaryStats.R"

# ---- Roar Cluster 路径 (取消注释以使用) ----
# data_path    <- "/storage/work/xjx5093/BurstIntercept/GoHIAR"
# model_file   <- file.path(data_path, "GoHIAR_BurstiVAR_4Burst_GroupCovFull.txt")
# source_file  <- "/storage/work/xjx5093/BurstIntercept/posteriorSummaryStats.R"

result_path <- file.path(data_path, "BurstiVAR_results_GroupCovFull")
if (!dir.exists(result_path)) dir.create(result_path, recursive = TRUE)

resulttable_file  <- file.path(result_path, "GoHIAR_BurstiVAR_4Burst_GroupCovFull_resulttable.csv")
descriptives_file <- file.path(result_path, "GoHIAR_BurstiVAR_4Burst_GroupCovFull_descriptives.csv")
codasamples_file  <- file.path(result_path, "GoHIAR_BurstiVAR_4Burst_GroupCovFull_codaSamples.RData")
Yarray_file       <- file.path(result_path, "GoHIAR_BurstiVAR_4Burst_GroupCovFull_Y.RData")
groupmeans_file   <- file.path(result_path, "GoHIAR_BurstiVAR_4Burst_GroupCovFull_GroupMeans.csv")
dynamics_file     <- file.path(result_path, "GoHIAR_BurstiVAR_4Burst_GroupCovFull_DynamicsByGroup.csv")

# ==============================================================================
# 2. 加载数据与辅助函数
# ==============================================================================

source(source_file)

setwd(data_path)
load("data_private/datblock_Meaning.Rdata")        # private input, available from the original study team
load("data_private/datblock_relationship.Rdata")   # private input, available from the original study team
load("data_private/datblock_accomplishment.Rdata") # private input, available from the original study team
load("data_private/datCov.Rdata")                  # private input, available from the original study team

cat("=== 数据加载完成 ===\n")
cat("  datblock 列名:", paste(colnames(datblock), collapse = ", "), "\n")
cat("  datCov  列名:", paste(colnames(datCov), collapse = ", "), "\n")
if (exists("dat")) cat("  dat     列名:", paste(colnames(dat), collapse = ", "), "\n")
cat("\n")

# ==============================================================================
# 3. 获取 Group 信息 & 筛选组
# ==============================================================================

if (exists("dat") && "Group" %in% colnames(dat)) {
  group_info <- unique(dat[, c("PID", "Group")])
  cat("Group 信息来源: dat\n")
} else if ("Group" %in% colnames(datCov)) {
  group_info <- unique(datCov[, c("PID", "Group")])
  cat("Group 信息来源: datCov\n")
} else {
  stop("找不到 Group 列, 请检查 dat 或 datCov")
}

cat("\n--- 原始样本的 Group 分布 ---\n")
print(table(group_info$Group))
cat("\n")

groups_to_keep <- c(1, 3)   # G1 Control + G3 PPI+Med
cat(sprintf("保留的组: %s\n\n", paste(groups_to_keep, collapse = ", ")))

group_info <- group_info[group_info$Group %in% groups_to_keep, ]

if (length(groups_to_keep) != 2) {
  stop("groups_to_keep 当前支持 2 组对比")
}
control_code      <- groups_to_keep[1]
intervention_code <- groups_to_keep[2]
group_info$G <- ifelse(group_info$Group == control_code, 0L, 1L)

N_ctrl <- sum(group_info$G == 0)
N_intv <- sum(group_info$G == 1)
cat(sprintf("  Control      (Group = %d, G = 0): N = %d\n", control_code, N_ctrl))
cat(sprintf("  Intervention (Group = %d, G = 1): N = %d\n", intervention_code, N_intv))
cat(sprintf("  Total kept                       : N = %d\n\n", N_ctrl + N_intv))

# ==============================================================================
# 4. 构建 Y 数组
# ==============================================================================

cat("=== 构建 Y 数组 ===\n\n")

T1 <- 56L
T2 <- 116L
T3 <- 168L
T4 <- 224L

ord <- order(group_info$G, group_info$PID)
group_info <- group_info[ord, ]
pInd <- group_info$PID
G    <- as.integer(group_info$G)
N    <- length(pInd)

cat("筛选后 N =", N, "(Control =", N_ctrl, ", Intervention =", N_intv, ")\n")
cat("Burst 边界: T1 =", T1, ", T2 =", T2, ", T3 =", T3, ", T4 =", T4, "\n\n")

Y <- array(NA, dim = c(N, T4, 2))

datblock$row_id <- ave(seq_len(nrow(datblock)),
                       datblock$PID, FUN = seq_along)
datblock_relationship$row_id <- ave(seq_len(nrow(datblock_relationship)),
                                    datblock_relationship$PID, FUN = seq_along)

obs_count <- matrix(0, N, 4)

for (pp in 1:N) {
  pid <- pInd[pp]

  tmp_m <- datblock[datblock$PID == pid, ]
  tmp_m <- tmp_m[order(tmp_m$row_id), ]
  tmp_r <- datblock_relationship[datblock_relationship$PID == pid, ]
  tmp_r <- tmp_r[order(tmp_r$row_id), ]
  nT_person <- nrow(tmp_m)

  for (oo in 1:min(nT_person, T4)) {
    Y[pp, oo, 1] <- tmp_m$Meaning[oo]
    if (oo <= nrow(tmp_r)) {
      Y[pp, oo, 2] <- tmp_r$Relationship[oo]
    }
  }

  for (b in 1:4) {
    if (b == 1) idx <- 1:T1
    if (b == 2) idx <- (T1+1):T2
    if (b == 3) idx <- (T2+1):T3
    if (b == 4) idx <- (T3+1):T4
    obs_count[pp, b] <- sum(!is.na(Y[pp, idx, 1]) & !is.na(Y[pp, idx, 2]))
  }
}

# Bivariate completeness
n_partial <- 0
for (pp in 1:N) {
  for (tt in 1:T4) {
    has_y1 <- !is.na(Y[pp, tt, 1])
    has_y2 <- !is.na(Y[pp, tt, 2])
    if (has_y1 != has_y2) {
      Y[pp, tt, 1] <- NA
      Y[pp, tt, 2] <- NA
      n_partial <- n_partial + 1
    }
  }
}

total_cells <- N * T4
obs_cells   <- sum(!is.na(Y[,,1]))
cat("部分观测 → 全 NA:", n_partial, "个 cell\n")
cat("Y 数组维度:", dim(Y), "\n")
cat(sprintf("总 cells: %d, 有效 bivariate: %d (%.1f%%)\n\n",
            total_cells, obs_cells, obs_cells / total_cells * 100))

save(Y, G, pInd, N, N_ctrl, N_intv, T1, T2, T3, T4, obs_count,
     file = Yarray_file)
cat("Y 数组 + G 已保存至:", Yarray_file, "\n\n")

# ==============================================================================
# 5. 描述性统计 (Group × Burst)
# ==============================================================================

cat("=== 描述性统计: Group × Burst bivariate-complete timepoints ===\n\n")

group_labels <- c(
  paste0("Control (Group=", control_code, ")"),
  paste0("Intervention (Group=", intervention_code, ")")
)

desc_rows <- list()

for (g in c(-1, 0, 1)) {
  if (g == -1) {
    pp_idx <- 1:N
    grp_lab <- sprintf("Pooled (N=%d)", N)
  } else {
    pp_idx <- which(G == g)
    grp_lab <- sprintf("%s (N=%d)", group_labels[g + 1], length(pp_idx))
  }

  cat(sprintf("--- %s ---\n", grp_lab))

  for (b in 1:4) {
    burst_lab <- c("P1:Pre", "P2:Intervention", "P3:Post-I", "P4:Post-II")[b]
    vals <- obs_count[pp_idx, b]
    cat(sprintf("  Burst %d (%s): Mean=%.1f, Median=%.1f, Min=%d, Max=%d, SD=%.1f\n",
                b, burst_lab, mean(vals), median(vals),
                min(vals), max(vals), sd(vals)))

    desc_rows[[length(desc_rows) + 1]] <- data.frame(
      Group  = grp_lab,
      Burst  = b,
      Phase  = burst_lab,
      N      = length(pp_idx),
      Mean   = round(mean(vals), 2),
      Median = median(vals),
      Min    = min(vals),
      Max    = max(vals),
      SD     = round(sd(vals), 2),
      stringsAsFactors = FALSE
    )
  }
  cat("\n")
}

desc_df <- do.call(rbind, desc_rows)
write.csv(desc_df, descriptives_file, row.names = FALSE)
cat("描述性表已保存至:", descriptives_file, "\n\n")

# ==============================================================================
# 6. JAGS 数据列表
# ==============================================================================

W0_val <- 0.015
W0 <- diag(rep(W0_val, 2))

jags_data <- list(
  Y      = Y,
  G      = G,
  P      = N,
  W0     = W0,
  df_res = 3L,
  T1     = T1,
  T2     = T2,
  T3     = T3,
  T4     = T4
)

cat("=== JAGS 数据列表 ===\n")
cat("  P =", N, "(Control =", N_ctrl, ", Intervention =", N_intv, ")\n")
cat("  T1 =", T1, ", T2 =", T2, ", T3 =", T3, ", T4 =", T4, "\n\n")

# ==============================================================================
# 7. 监测参数
# ==============================================================================

parameters <- c(
  # --- Intercept level ---
  "mu_base",
  "delta_main",
  "delta_grp",
  "mu_ctrl",
  "mu_intv",
  # --- Dynamics level (NEW) ---
  "mu_AR_base",
  "mu_CR_base",
  "delta_grp_AR",
  "delta_grp_CR",
  "AR_ctrl",
  "AR_intv",
  "CR_ctrl",
  "CR_intv",
  # --- Variance components ---
  "Level2Sigma",
  "sigma_innovation"
)

# ==============================================================================
# 8. 初始值
# ==============================================================================

inits <- list(
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 2026),
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 2027)
)

# ==============================================================================
# 9. 运行 JAGS
# ==============================================================================

cat("=== 开始拟合 GoHIAR_BurstiVAR_4Burst_GroupCovFull ===\n")
cat("模型文件:", model_file, "\n\n")

if (!file.exists(model_file)) {
  stop("模型文件不存在: ", model_file)
}

cat("Step 1/3: jags.model (adapt=5000)...\n")
t_start <- proc.time()
jagsModel <- jags.model(
  file     = model_file,
  data     = jags_data,
  inits    = inits,
  n.chains = 2,
  n.adapt  = 5000,
  quiet    = FALSE
)
cat("  adapt 完成, 用时:",
    round((proc.time() - t_start)[3] / 60, 1), "分钟\n\n")

cat("Step 2/3: Burn-in (5000)...\n")
t_burnin <- proc.time()
update(jagsModel, n.iter = 5000, progress.bar = "text")
cat("  Burn-in 完成, 用时:",
    round((proc.time() - t_burnin)[3] / 60, 1), "分钟\n\n")

cat("Step 3/3: 采样 (20000, thin=1)...\n")
t_sample <- proc.time()
codaSamples <- coda.samples(
  model          = jagsModel,
  variable.names = parameters,
  n.iter         = 20000,
  thin           = 1,
  progress.bar   = "text"
)
cat("  采样完成, 用时:",
    round((proc.time() - t_sample)[3] / 60, 1), "分钟\n\n")

save(codaSamples, file = codasamples_file)
cat("MCMC samples 已保存至:", codasamples_file, "\n")

# ==============================================================================
# 10. 收敛诊断
# ==============================================================================

cat("\n=== 收敛诊断 ===\n")
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
  cat("  *** 警告: 未收敛参数:",
      paste(head(bad, 10), collapse = ", "), "\n")
}

# ==============================================================================
# 11. 结果表 + 关键参数打印
# ==============================================================================

resulttable <- summarizePost(codaSamples)
write.csv(resulttable, resulttable_file)
cat("\n结果表已保存至:", resulttable_file, "\n")

samples_mat <- do.call(rbind, codaSamples)
par_names   <- colnames(samples_mat)
sig_star <- function(ll, ul) ifelse(ll > 0 | ul < 0, "***", "   ")

# --- 11a. mu_base (intercept baseline) ---
cat("\n--- mu_base: Burst 1 intercept baseline (跨组共享) ---\n")
for (k in 1:2) {
  p <- paste0("mu_base[", k, "]")
  if (p %in% par_names) {
    post <- samples_mat[, p]
    vn <- ifelse(k == 1, "Meaning", "Relationship")
    cat(sprintf("  mu_base[%d] (%-13s): Mean = %6.2f, 95%% CI = [%6.2f, %6.2f]\n",
                k, vn, mean(post),
                quantile(post, 0.025), quantile(post, 0.975)))
  }
}

# --- 11b. delta_main (Control intercept phase effects) ---
cat("\n--- delta_main: Control 组相对 Burst 1 的 intercept phase effects ---\n")
dm_labels <- c(
  "B2-B1 y1 Meaning",     "B2-B1 y2 Relationship",
  "B3-B1 y1 Meaning",     "B3-B1 y2 Relationship",
  "B4-B1 y1 Meaning",     "B4-B1 y2 Relationship"
)
for (k in 1:6) {
  p <- paste0("delta_main[", k, "]")
  if (p %in% par_names) {
    post <- samples_mat[, p]
    ll <- quantile(post, 0.025); ul <- quantile(post, 0.975)
    cat(sprintf("  delta_main[%d] %-25s: Mean = %6.2f, 95%% CI = [%6.2f, %6.2f] %s\n",
                k, dm_labels[k], mean(post), ll, ul, sig_star(ll, ul)))
  }
}

# --- 11c. delta_grp (intercept group × phase interaction) ---
cat("\n--- delta_grp: Intervention 的 *额外* intercept phase effects ---\n")
for (k in 1:6) {
  p <- paste0("delta_grp[", k, "]")
  if (p %in% par_names) {
    post <- samples_mat[, p]
    ll <- quantile(post, 0.025); ul <- quantile(post, 0.975)
    cat(sprintf("  delta_grp[%d]  %-25s: Mean = %6.2f, 95%% CI = [%6.2f, %6.2f] %s\n",
                k, dm_labels[k], mean(post), ll, ul, sig_star(ll, ul)))
  }
}

# --- 11d. mu_ctrl / mu_intv: 两组 8 burst means ---
cat("\n--- mu_ctrl: Control 组 8 个 burst means ---\n")
burst_labels_long <- c(
  "B1 P1:Pre  y1 Meaning",       "B1 P1:Pre  y2 Relationship",
  "B2 P2:Intv y1 Meaning",       "B2 P2:Intv y2 Relationship",
  "B3 P3:PI   y1 Meaning",       "B3 P3:PI   y2 Relationship",
  "B4 P4:PII  y1 Meaning",       "B4 P4:PII  y2 Relationship"
)
for (k in 1:8) {
  p <- paste0("mu_ctrl[", k, "]")
  if (p %in% par_names) {
    post <- samples_mat[, p]
    cat(sprintf("  mu_ctrl[%d] %-32s: Mean = %6.2f, 95%% CI = [%6.2f, %6.2f]\n",
                k, burst_labels_long[k],
                mean(post), quantile(post, 0.025), quantile(post, 0.975)))
  }
}
cat("\n--- mu_intv: Intervention 组 8 个 burst means ---\n")
for (k in 1:8) {
  p <- paste0("mu_intv[", k, "]")
  if (p %in% par_names) {
    post <- samples_mat[, p]
    cat(sprintf("  mu_intv[%d] %-32s: Mean = %6.2f, 95%% CI = [%6.2f, %6.2f]\n",
                k, burst_labels_long[k],
                mean(post), quantile(post, 0.025), quantile(post, 0.975)))
  }
}

# --- 11e. Dynamics: mu_AR_base, mu_CR_base (Control) ---
cat("\n--- mu_AR_base / mu_CR_base: Control 组 dynamics 均值 ---\n")
ar_labels <- c("AR y1 Meaning", "AR y2 Relationship")
cr_labels <- c("CR y1(Meaning) -> y2(Relationship)",
               "CR y2(Relationship) -> y1(Meaning)")
for (k in 1:2) {
  p <- paste0("mu_AR_base[", k, "]")
  if (p %in% par_names) {
    post <- samples_mat[, p]
    ll <- quantile(post, 0.025); ul <- quantile(post, 0.975)
    cat(sprintf("  mu_AR_base[%d] %-38s: Mean = %6.3f, 95%% CI = [%6.3f, %6.3f] %s\n",
                k, ar_labels[k], mean(post), ll, ul, sig_star(ll, ul)))
  }
}
for (k in 1:2) {
  p <- paste0("mu_CR_base[", k, "]")
  if (p %in% par_names) {
    post <- samples_mat[, p]
    ll <- quantile(post, 0.025); ul <- quantile(post, 0.975)
    cat(sprintf("  mu_CR_base[%d] %-38s: Mean = %6.3f, 95%% CI = [%6.3f, %6.3f] %s\n",
                k, cr_labels[k], mean(post), ll, ul, sig_star(ll, ul)))
  }
}

# --- 11f. delta_grp_AR / delta_grp_CR (group 对 dynamics 的效应; 主角之一) ---
cat("\n--- delta_grp_AR: Intervention 的 *额外* AR 差异 (group effect on AR) ---\n")
for (k in 1:2) {
  p <- paste0("delta_grp_AR[", k, "]")
  if (p %in% par_names) {
    post <- samples_mat[, p]
    ll <- quantile(post, 0.025); ul <- quantile(post, 0.975)
    cat(sprintf("  delta_grp_AR[%d] %-36s: Mean = %6.3f, 95%% CI = [%6.3f, %6.3f] %s\n",
                k, ar_labels[k], mean(post), ll, ul, sig_star(ll, ul)))
  }
}
cat("\n--- delta_grp_CR: Intervention 的 *额外* CR 差异 (group effect on CR) ---\n")
for (k in 1:2) {
  p <- paste0("delta_grp_CR[", k, "]")
  if (p %in% par_names) {
    post <- samples_mat[, p]
    ll <- quantile(post, 0.025); ul <- quantile(post, 0.975)
    cat(sprintf("  delta_grp_CR[%d] %-36s: Mean = %6.3f, 95%% CI = [%6.3f, %6.3f] %s\n",
                k, cr_labels[k], mean(post), ll, ul, sig_star(ll, ul)))
  }
}

# --- 11g. AR_ctrl vs AR_intv, CR_ctrl vs CR_intv (派生) ---
cat("\n--- AR by group (派生) ---\n")
for (k in 1:2) {
  pc <- samples_mat[, paste0("AR_ctrl[", k, "]")]
  pI <- samples_mat[, paste0("AR_intv[", k, "]")]
  cat(sprintf("  %-20s Ctrl: %6.3f [%6.3f, %6.3f]  |  Intv: %6.3f [%6.3f, %6.3f]\n",
              ar_labels[k],
              mean(pc), quantile(pc, 0.025), quantile(pc, 0.975),
              mean(pI), quantile(pI, 0.025), quantile(pI, 0.975)))
}
cat("\n--- CR by group (派生) ---\n")
for (k in 1:2) {
  pc <- samples_mat[, paste0("CR_ctrl[", k, "]")]
  pI <- samples_mat[, paste0("CR_intv[", k, "]")]
  cat(sprintf("  %-38s Ctrl: %6.3f [%6.3f, %6.3f]  |  Intv: %6.3f [%6.3f, %6.3f]\n",
              cr_labels[k],
              mean(pc), quantile(pc, 0.025), quantile(pc, 0.975),
              mean(pI), quantile(pI, 0.025), quantile(pI, 0.975)))
}

# --- 11h. Innovation covariance ---
cat("\n--- Innovation Covariance Matrix ---\n")
for (i in 1:2) {
  for (j in i:2) {
    p <- paste0("sigma_innovation[", i, ",", j, "]")
    if (p %in% par_names) {
      post <- samples_mat[, p]
      vname <- ifelse(i == j,
                      paste0("Var(", c("Meaning","Relationship")[i], ")"),
                      "Cov(Meaning, Relationship)")
      cat(sprintf("  sigma_innovation[%d,%d] %-30s: Mean = %7.2f, 95%% CI = [%7.2f, %7.2f]\n",
                  i, j, vname,
                  mean(post), quantile(post, 0.025), quantile(post, 0.975)))
    }
  }
}

# ==============================================================================
# 12. 导出: 两组 burst means 对比表 (画图用)
# ==============================================================================

gm_rows <- list()
for (k in 1:8) {
  if (k %% 2 == 1) {
    var_lab <- "Meaning"; burst_num <- (k + 1) / 2
  } else {
    var_lab <- "Relationship"; burst_num <- k / 2
  }
  burst_phase <- c("P1:Pre", "P2:Intv", "P3:Post-I", "P4:Post-II")[burst_num]

  pc <- samples_mat[, paste0("mu_ctrl[", k, "]")]
  pI <- samples_mat[, paste0("mu_intv[", k, "]")]

  gm_rows[[length(gm_rows) + 1]] <- data.frame(
    Burst = burst_num, Phase = burst_phase, Variable = var_lab,
    Group = "Control",
    Mean = round(mean(pc), 3),
    LL95 = round(quantile(pc, 0.025), 3),
    UL95 = round(quantile(pc, 0.975), 3),
    stringsAsFactors = FALSE
  )
  gm_rows[[length(gm_rows) + 1]] <- data.frame(
    Burst = burst_num, Phase = burst_phase, Variable = var_lab,
    Group = "Intervention",
    Mean = round(mean(pI), 3),
    LL95 = round(quantile(pI, 0.025), 3),
    UL95 = round(quantile(pI, 0.975), 3),
    stringsAsFactors = FALSE
  )
}
gm_df <- do.call(rbind, gm_rows)
write.csv(gm_df, groupmeans_file, row.names = FALSE)
cat("\n两组 burst-mean 对比表已保存至:", groupmeans_file, "\n")

# ==============================================================================
# 13. 导出: 两组 dynamics 对比表 (NEW)
# ==============================================================================

dyn_rows <- list()
dyn_labels <- list(
  AR = c("AR y1 Meaning", "AR y2 Relationship"),
  CR = c("CR y1->y2", "CR y2->y1")
)
for (type in c("AR", "CR")) {
  for (k in 1:2) {
    pc <- samples_mat[, paste0(type, "_ctrl[", k, "]")]
    pI <- samples_mat[, paste0(type, "_intv[", k, "]")]
    pdiff <- pI - pc
    ll_d <- quantile(pdiff, 0.025); ul_d <- quantile(pdiff, 0.975)

    dyn_rows[[length(dyn_rows) + 1]] <- data.frame(
      Parameter = dyn_labels[[type]][k],
      Type      = type,
      Index     = k,
      Ctrl_Mean = round(mean(pc), 4),
      Ctrl_LL95 = round(quantile(pc, 0.025), 4),
      Ctrl_UL95 = round(quantile(pc, 0.975), 4),
      Intv_Mean = round(mean(pI), 4),
      Intv_LL95 = round(quantile(pI, 0.025), 4),
      Intv_UL95 = round(quantile(pI, 0.975), 4),
      Diff_Mean = round(mean(pdiff), 4),
      Diff_LL95 = round(ll_d, 4),
      Diff_UL95 = round(ul_d, 4),
      Sig       = ifelse(ll_d > 0 | ul_d < 0, "***", "   "),
      stringsAsFactors = FALSE
    )
  }
}
dyn_df <- do.call(rbind, dyn_rows)
write.csv(dyn_df, dynamics_file, row.names = FALSE)
cat("两组 dynamics 对比表已保存至:", dynamics_file, "\n")

# ==============================================================================
# 14. 收尾
# ==============================================================================

cat("\n=== 全部完成 ===\n")
cat("产出文件:\n")
cat("  1. Y array          :", Yarray_file, "\n")
cat("  2. 描述性表 (CSV)   :", descriptives_file, "\n")
cat("  3. MCMC samples     :", codasamples_file, "\n")
cat("  4. 后验汇总表 (CSV) :", resulttable_file, "\n")
cat("  5. 两组 burst means :", groupmeans_file, "\n")
cat("  6. 两组 dynamics    :", dynamics_file, "\n\n")

cat("海报/论文关键参数:\n")
cat("  - delta_grp[1:6]    (group x phase interaction, intercepts)\n")
cat("  - delta_grp_AR[1:2] (group effect on AR) ** 新 **\n")
cat("  - delta_grp_CR[1:2] (group effect on CR) ** 新 **\n")
cat("  - mu_ctrl vs mu_intv, AR_ctrl vs AR_intv, CR_ctrl vs CR_intv -> 可视化\n")
