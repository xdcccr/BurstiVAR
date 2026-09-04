# ------------------------------------------------------------------------------
# Provenance header added for the public BurstiVAR repository.
# Original file: IP_1b/MOL/GoHIAR_BurstiVAR_4Burst_BaselineDelta_Modelfit.R
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
# GoHIAR_BurstiVAR_4Burst_BaselineDelta_Modelfit.R
#
# 功能: 用 4-Burst BurstiVAR (Baseline + Delta 参数化) 拟合 Go-HIAR 经验数据
#       y1 = Meaning of life (MOL), 0-100 scale
#       y2 = Relationship, 0-100 scale
#       Burst 1 (P1: Pre-intervention) 作为 baseline，其他 burst = baseline + delta
#
# 新增参数:
#   mu_base[1:2]: Burst 1 (Pre-intervention) 的 population mean (y1, y2)
#   delta[1:6]:   其他 burst 相对 Burst 1 的偏移
#                 delta[1:2] = Burst2 - Burst1 (Intervention effect)
#                 delta[3:4] = Burst3 - Burst1 (Post-I vs Pre)
#                 delta[5:6] = Burst4 - Burst1 (Post-II vs Pre)
#
# Burst 结构 (4 study periods × 4 blocks/day):
#   T1 = 56   (Burst 1 结束: P1 Pre-intervention)
#   T2 = 116  (Burst 2 结束: P2 Intervention)
#   T3 = 168  (Burst 3 结束: P3 Post-intervention I)
#   T4 = 224  (Burst 4 结束: P4 Post-intervention II)
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
model_file   <- file.path(data_path, "GoHIAR_BurstiVAR_4Burst_BaselineDelta.txt")
source_file  <- "./posteriorSummaryStats.R"

# ---- Roar Cluster 路径 (取消注释以使用) ----
# data_path    <- "/storage/work/xjx5093/BurstIntercept/GoHIAR"
# model_file   <- file.path(data_path, "GoHIAR_BurstiVAR_4Burst_BaselineDelta.txt")
# source_file  <- "/storage/work/xjx5093/BurstIntercept/posteriorSummaryStats.R"

# 输出路径
result_path <- file.path(data_path, "BurstiVAR_results_BaselineDelta")
if (!dir.exists(result_path)) dir.create(result_path, recursive = TRUE)

resulttable_file <- file.path(result_path, "GoHIAR_BurstiVAR_4Burst_BaselineDelta_resulttable.csv")
codasamples_file <- file.path(result_path, "GoHIAR_BurstiVAR_4Burst_BaselineDelta_codaSamples.RData")
Yarray_file      <- file.path(result_path, "GoHIAR_BurstiVAR_4Burst_BaselineDelta_Y.RData")
diff_file        <- file.path(result_path, "GoHIAR_BurstiVAR_4Burst_BaselineDelta_SDdiff.csv")

# ==============================================================================
# 2. 加载数据与辅助函数
# ==============================================================================

source(source_file)   # 加载 summarizePost()

setwd(data_path)
load("data_private/datblock_Meaning.Rdata")        # private input, available from the original study team
load("data_private/datblock_relationship.Rdata")   # private input, available from the original study team
load("data_private/datblock_accomplishment.Rdata") # private input, available from the original study team
load("data_private/datCov.Rdata")                  # private input, available from the original study team

cat("=== 数据加载完成 ===\n")
cat("  datblock 列名:", paste(colnames(datblock), collapse = ", "), "\n")
cat("  N:", ifelse(exists("N"), N, "未找到"), "\n")
cat("  pInd:", ifelse(exists("pInd"), paste0("found (", length(pInd), " persons)"), "未找到"), "\n\n")

# ==============================================================================
# 3. 数据准备：构建 Y[P, T, 2] 数组
# ==============================================================================

cat("=== 构建 Y 数组 ===\n\n")

# Burst 时间边界
T1 <- 56L     # Burst 1 结束 (P1: Pre-intervention)
T2 <- 116L    # Burst 2 结束 (P2: Intervention)
T3 <- 168L    # Burst 3 结束 (P3: Post-intervention I)
T4 <- 224L    # Burst 4 结束 (P4: Post-intervention II)

if (!exists("pInd")) pInd <- unique(datblock$PID)
if (!exists("N"))    N <- length(pInd)

cat("人数 P =", N, "\n")
cat("总时间点 T =", T4, "\n")
cat("Burst 边界: T1 =", T1, ", T2 =", T2, ", T3 =", T3, ", T4 =", T4, "\n\n")

# 初始化 Y 数组
Y <- array(NA, dim = c(N, T4, 2))

# 添加 row_id 用于对齐
datblock$row_id <- ave(seq_len(nrow(datblock)), datblock$PID, FUN = seq_along)
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

# Bivariate completeness: 部分观测 → 全 NA
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
miss_cells  <- total_cells - obs_cells

cat("部分观测设为全 NA:", n_partial, "个时间点\n")
cat("Y 数组维度:", dim(Y), "\n")
cat(sprintf("总 cells: %d, 有效: %d (%.1f%%), 缺失: %d (%.1f%%)\n",
            total_cells, obs_cells, obs_cells/total_cells*100,
            miss_cells, miss_cells/total_cells*100))

cat("\n每 burst 平均有效 bivariate 观测数:\n")
for (b in 1:4) {
  cat(sprintf("  Burst %d (P%d): M = %.1f, range = [%d, %d]\n",
              b, b, mean(obs_count[,b]), min(obs_count[,b]), max(obs_count[,b])))
}

save(Y, pInd, N, T1, T2, T3, T4, obs_count, file = Yarray_file)
cat("\nY 数组已保存至:", Yarray_file, "\n")

# ==============================================================================
# 4. JAGS 数据列表
# ==============================================================================

W0_val <- 0.015
W0 <- diag(rep(W0_val, 2))

jags_data <- list(
  Y      = Y,
  P      = N,
  W0     = W0,
  df_res = 3L,
  T1     = T1,
  T2     = T2,
  T3     = T3,
  T4     = T4
)

cat("\n=== JAGS 数据列表 ===\n")
cat("  P =", N, ", T1 =", T1, ", T2 =", T2, ", T3 =", T3, ", T4 =", T4, "\n")
cat("  W0 diagonal =", W0_val, ", df_res =", jags_data$df_res, "\n\n")

# ==============================================================================
# 5. 监测参数
# ==============================================================================

parameters <- c(
  "mu_base",           # [1:2]: Burst 1 baseline means (y1, y2)
  "delta",             # [1:6]: 其他 burst 相对 Burst 1 的偏移
  "Level2Mean",        # [1:12]: 所有 population means (验证用)
  "Level2Sigma",       # [1:12]: 所有 between-person SDs
  "sigma_innovation"   # [2,2]: innovation covariance matrix
)

# ==============================================================================
# 6. 初始值
# ==============================================================================

inits <- list(
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 2026),
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 2027)
)

# ==============================================================================
# 7. 运行 JAGS
# ==============================================================================

cat("=== 开始拟合 GoHIAR_BurstiVAR_4Burst_BaselineDelta ===\n")
cat("模型文件:", model_file, "\n\n")

if (!file.exists(model_file)) {
  stop("模型文件不存在: ", model_file,
       "\n请确保 GoHIAR_BurstiVAR_4Burst_BaselineDelta.txt 已复制到 ", data_path)
}

# Step 1: 建立模型
cat("Step 1/3: 建立 JAGS 模型 (adapt=5000)...\n")
t_start <- proc.time()

jagsModel <- jags.model(
  file     = model_file,
  data     = jags_data,
  inits    = inits,
  n.chains = 2,
  n.adapt  = 5000,
  quiet    = FALSE
)
cat("  适应阶段完成，用时:", round((proc.time() - t_start)[3]/60, 1), "分钟\n\n")

# Step 2: Burn-in
cat("Step 2/3: Burn-in (5000 iterations)...\n")
t_burnin <- proc.time()
update(jagsModel, n.iter = 5000, progress.bar = "text")
cat("  Burn-in 完成，用时:", round((proc.time() - t_burnin)[3]/60, 1), "分钟\n\n")

# Step 3: 采样
cat("Step 3/3: 正式采样 (20000 iterations, thin=1)...\n")
t_sample <- proc.time()
codaSamples <- coda.samples(
  model          = jagsModel,
  variable.names = parameters,
  n.iter         = 20000,
  thin           = 1,
  progress.bar   = "text"
)
cat("  采样完成，用时:", round((proc.time() - t_sample)[3]/60, 1), "分钟\n\n")

save(codaSamples, file = codasamples_file)
cat("MCMC samples 已保存至:", codasamples_file, "\n")

# ==============================================================================
# 8. 收敛诊断
# ==============================================================================

cat("\n=== 收敛诊断 ===\n")

rhat_vals <- gelman.diag(codaSamples, multivariate = FALSE)$psrf[, "Point est."]
cat("R-hat 统计:\n")
cat("  最大 R-hat:", round(max(rhat_vals, na.rm = TRUE), 4), "\n")
cat("  均值 R-hat:", round(mean(rhat_vals, na.rm = TRUE), 4), "\n")
cat("  R-hat > 1.05 的参数数量:", sum(rhat_vals > 1.05, na.rm = TRUE), "\n")

if (max(rhat_vals, na.rm = TRUE) > 1.10) {
  cat("  *** 警告: 部分参数可能未收敛。建议增加迭代次数。***\n")
  bad_params <- names(rhat_vals[rhat_vals > 1.10])
  cat("  问题参数:", paste(head(bad_params, 10), collapse = ", "), "\n")
}

ess_vals <- effectiveSize(codaSamples)
cat("\nESS 统计:\n")
cat("  最小 ESS:", round(min(ess_vals, na.rm = TRUE), 0), "\n")
cat("  均值 ESS:", round(mean(ess_vals, na.rm = TRUE), 0), "\n")
cat("  ESS < 200 的参数数量:", sum(ess_vals < 200, na.rm = TRUE), "\n")

# ==============================================================================
# 9. 提取结果表
# ==============================================================================

cat("\n=== 生成结果汇总表 ===\n")
resulttable <- summarizePost(codaSamples)
write.csv(resulttable, resulttable_file)
cat("结果表已保存至:", resulttable_file, "\n")

# ==============================================================================
# 10. 关键参数提取与展示
# ==============================================================================

cat("\n=== 关键参数估计结果 ===\n")

samples_mat <- do.call(rbind, codaSamples)
par_names <- colnames(samples_mat)

# --- 10a. Baseline (Burst 1, Pre-intervention) ---
cat("\n--- Baseline: Burst 1 (P1: Pre-intervention) ---\n")
for (j in 1:2) {
  var_label <- ifelse(j == 1, "Meaning (y1)", "Relationship (y2)")
  param     <- paste0("mu_base[", j, "]")
  if (param %in% par_names) {
    post  <- samples_mat[, param]
    cat(sprintf("  mu_base[%d] %-25s | Mean = %6.2f, 95%% CI = [%6.2f, %6.2f]\n",
                j, paste0("(", var_label, ")"),
                mean(post), quantile(post, 0.025), quantile(post, 0.975)))
  }
}

# --- 10b. Delta (其他 burst vs Burst 1) ---
cat("\n--- Delta: 其他 Burst 相对 Burst 1 (Pre-intervention) 的偏移 ---\n")
cat("  (正值 = 该 burst 均值高于 baseline)\n\n")
delta_labels <- c(
  "B2-B1 y1 Meaning     (Intervention effect)",
  "B2-B1 y2 Relationship (Intervention effect)",
  "B3-B1 y1 Meaning     (Post-I vs Pre)",
  "B3-B1 y2 Relationship (Post-I vs Pre)",
  "B4-B1 y1 Meaning     (Post-II vs Pre)",
  "B4-B1 y2 Relationship (Post-II vs Pre)"
)
for (d in 1:6) {
  param <- paste0("delta[", d, "]")
  if (param %in% par_names) {
    post <- samples_mat[, param]
    ll   <- quantile(post, 0.025)
    ul   <- quantile(post, 0.975)
    sig  <- ifelse(ll > 0 | ul < 0, "***", "   ")
    cat(sprintf("  delta[%d] %-50s | Mean = %6.2f, 95%% CI = [%6.2f, %6.2f] %s\n",
                d, delta_labels[d], mean(post), ll, ul, sig))
  }
}

# --- 10c. Level2Mean 重构验证 ---
cat("\n--- Level2Mean[1:8] 重构验证 ---\n")
burst_labels <- c(
  "Burst1 y1 Meaning     (baseline)",
  "Burst1 y2 Relationship (baseline)",
  "Burst2 y1 Meaning     (base+delta[1])",
  "Burst2 y2 Relationship (base+delta[2])",
  "Burst3 y1 Meaning     (base+delta[3])",
  "Burst3 y2 Relationship (base+delta[4])",
  "Burst4 y1 Meaning     (base+delta[5])",
  "Burst4 y2 Relationship (base+delta[6])"
)
for (k in 1:8) {
  param <- paste0("Level2Mean[", k, "]")
  if (param %in% par_names) {
    post <- samples_mat[, param]
    cat(sprintf("  Level2Mean[%2d] %-45s | Mean = %6.2f, 95%% CI = [%6.2f, %6.2f]\n",
                k, burst_labels[k],
                mean(post), quantile(post, 0.025), quantile(post, 0.975)))
  }
}

# --- 10d. VAR parameters ---
cat("\n--- VAR Parameters ---\n")
var_labels <- c(
  "AR y1 (Meaning inertia)",
  "AR y2 (Relationship inertia)",
  "CR y1->y2 (Meaning -> Relationship)",
  "CR y2->y1 (Relationship -> Meaning)"
)
for (k in 9:12) {
  param <- paste0("Level2Mean[", k, "]")
  if (param %in% par_names) {
    post <- samples_mat[, param]
    sig  <- ifelse(quantile(post, 0.025) > 0 | quantile(post, 0.975) < 0, "***", "   ")
    cat(sprintf("  Level2Mean[%2d] %-40s | Mean = %6.3f, 95%% CI = [%6.3f, %6.3f] %s\n",
                k, var_labels[k - 8],
                mean(post), quantile(post, 0.025), quantile(post, 0.975), sig))
  }
}

# --- 10e. Between-person SDs ---
cat("\n--- Between-Person SDs ---\n")
all_labels <- c(burst_labels, var_labels)
for (k in 1:12) {
  param <- paste0("Level2Sigma[", k, "]")
  if (param %in% par_names) {
    post <- samples_mat[, param]
    cat(sprintf("  Level2Sigma[%2d] %-45s | Mean = %6.3f, 95%% CI = [%6.3f, %6.3f]\n",
                k, all_labels[k],
                mean(post), quantile(post, 0.025), quantile(post, 0.975)))
  }
}

# --- 10f. Innovation covariance ---
cat("\n--- Innovation Covariance Matrix ---\n")
for (i in 1:2) {
  for (j in i:2) {
    param <- paste0("sigma_innovation[", i, ",", j, "]")
    if (param %in% par_names) {
      post <- samples_mat[, param]
      vname <- ifelse(i == j,
                      paste0("Var(", c("Meaning","Relationship")[i], ")"),
                      "Cov(Meaning, Relationship)")
      cat(sprintf("  sigma_innovation[%d,%d] %-30s | Mean = %7.2f, 95%% CI = [%7.2f, %7.2f]\n",
                  i, j, vname,
                  mean(post), quantile(post, 0.025), quantile(post, 0.975)))
    }
  }
}

# ==============================================================================
# 11. 后处理：Between-person SD 差值
#     (sigma_burst_k - sigma_baseline, 检验各 burst 间个体差异是否不同)
# ==============================================================================

cat("\n=== SD 差值后处理 (sigma_burst_k - sigma_baseline) ===\n")
cat("  Baseline SD = Level2Sigma[1] (y1) 和 Level2Sigma[2] (y2)\n\n")

sd_diff_results <- data.frame()

sd_comparisons <- list(
  list(k = 3, base = 1, label = "Burst2_y1 - Burst1_y1 (Meaning)"),
  list(k = 4, base = 2, label = "Burst2_y2 - Burst1_y2 (Relationship)"),
  list(k = 5, base = 1, label = "Burst3_y1 - Burst1_y1 (Meaning)"),
  list(k = 6, base = 2, label = "Burst3_y2 - Burst1_y2 (Relationship)"),
  list(k = 7, base = 1, label = "Burst4_y1 - Burst1_y1 (Meaning)"),
  list(k = 8, base = 2, label = "Burst4_y2 - Burst1_y2 (Relationship)")
)

for (comp in sd_comparisons) {
  p_k    <- paste0("Level2Sigma[", comp$k, "]")
  p_base <- paste0("Level2Sigma[", comp$base, "]")
  
  if (p_k %in% par_names && p_base %in% par_names) {
    diff_post <- samples_mat[, p_k] - samples_mat[, p_base]
    ll  <- quantile(diff_post, 0.025)
    ul  <- quantile(diff_post, 0.975)
    sig <- ifelse(ll > 0 | ul < 0, "***", "   ")
    
    row <- data.frame(
      Label     = comp$label,
      Mean_diff = round(mean(diff_post), 4),
      SD_diff   = round(sd(diff_post), 4),
      LL_95     = round(ll, 4),
      UL_95     = round(ul, 4),
      Sig       = sig,
      stringsAsFactors = FALSE
    )
    sd_diff_results <- rbind(sd_diff_results, row)
    
    cat(sprintf("  %-50s | Mean = %6.3f, 95%% CI = [%6.3f, %6.3f] %s\n",
                comp$label, mean(diff_post), ll, ul, sig))
  }
}

write.csv(sd_diff_results, diff_file, row.names = FALSE)
cat("\nSD 差值表已保存至:", diff_file, "\n")

# ==============================================================================
# 12. 完成
# ==============================================================================

cat("\n=== 运行完成 ===\n")
cat("输出文件:\n")
cat("  结果汇总表:", resulttable_file, "\n")
cat("  MCMC samples:", codasamples_file, "\n")
cat("  Y 数组:", Yarray_file, "\n")
cat("  SD 差值表:", diff_file, "\n")
cat("总用时:", round((proc.time() - t_start)[3] / 60, 1), "分钟\n")
