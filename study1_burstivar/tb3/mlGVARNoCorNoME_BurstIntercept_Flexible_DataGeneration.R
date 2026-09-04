# ----------------------------------------------------------------------------
# PUBLIC-REPO COPY. Original: IP_1b/Study1_Tb3_Corrected/mlGVARNoCorNoME_BurstIntercept_Flexible_DataGeneration.R
# Modified for this repository: the absolute work_path line was replaced with file.path(getwd(), "data")
# so the script runs from this directory (run it with the working directory
# set to study1_burstivar/tb3). All other content is identical to the version
# used for the corrected Tb = 3 re-run (2026-09-04) whose results appear in
# the manuscript.
# ----------------------------------------------------------------------------
# ============================================================================
# CORRECTED RE-RUN 2026-09-04 (Study1_Tb3_Corrected): the original 2026-05-09
# Tb=3 run (MLGVAR2605) generated BOTH variables' burst-intercept SDs as
# 1.0+(b-1)*0.5, i.e. y2 SDs (1,1.5,2) instead of the designed (1,2,3).
# This copy differs from the original ONLY in: (1) the y2 sigma line,
# (2) the work/base path, (3) [Modelfit only] num_cores 6 -> 12.
# Seeds, MCMC settings, and everything else are byte-identical.
# ============================================================================
######################## Set Environment ########################
rm(list=ls())
library(mvtnorm)
library(Matrix)
library(dplyr)

# ==============================================================================
# 灵活的数据生成脚本 - 支持任意数量的burst
# 用法: 修改 nBurst 和 nT_per_burst 参数即可
# ==============================================================================

######################## USER CONFIGURATION ########################
# ★★★ 核心参数：在此设置burst数量和每个burst的时间点数 ★★★
nBurst <- 3           # burst数量 (例如: 2, 3, 5, ...)
nT_per_burst <- 3     # 每个burst内的时间点数

# 工作路径设置
work_path <- file.path(getwd(), "data")  # REPO EDIT: run from this directory (was D:/XXYDATAanalysis/IP_1b/Study1_Tb3_Corrected/data)
print(paste0("work_path: ", work_path))

# 如果目录不存在则创建
if(!dir.exists(work_path)) {
  dir.create(work_path, recursive = TRUE)
}
setwd(work_path)

######################## Control Parameters ########################
# 条件名称（自动包含burst信息）
C <- paste0("mlGVARNoCorNoME_BurstIntercept_", nBurst, "Burst_", nT_per_burst, "T")

# 实验条件
nT_list <- c(nBurst * nT_per_burst)  # 总时间点 = burst数 × 每burst时间点
nP_list <- c(100, 500)               # 被试数
N_repl <- 100                        # Monte Carlo复制数

ydim <- 2    # 维度数 (固定为2: reading and math)
npad <- 0    # 初始丢弃的时间点数

######################## Model Parameters ########################

# ==============================================================================
# 函数：生成每个burst的Intercept参数
# ==============================================================================
generate_burst_parameters <- function(nBurst, ydim = 2) {
  # 为每个burst生成intercept均值和标准差
  # 这里使用一个简单的线性增长模式，你可以根据需要修改
  
  Intercept_mu_list <- list()
  Intercept_sigma_list <- list()
  
  for(b in 1:nBurst) {
    # 均值: 随burst递增
    Intercept_mu_list[[b]] <- c((b-1) * 1.0, (b-1) * 0.5 + 1.0)
    
    # 标准差: 可以固定或递增
    Intercept_sigma_list[[b]] <- c(1.0 + (b-1) * 0.5, b * 1.0)  # CORRECTED: y1 SDs (1,1.5,2), y2 SDs (1,2,3) per the Table 1 design
  }
  
  return(list(mu = Intercept_mu_list, sigma = Intercept_sigma_list))
}

# 生成burst参数
burst_params <- generate_burst_parameters(nBurst, ydim)

# VAR Parameters
AR_mean1 <- 0.3      # AR for reading
AR_mean2 <- 0.2      # AR for math
CR12_mean <- -0.15   # Cross effect: reading -> math
CR21_mean <- 0       # Cross effect: math -> reading (设为0)

# Standard deviations for VAR parameters
sigma_AR <- c(0.1, 0.1)    # SD for AR coefficients
sigma_CR <- c(0.1, 0.1)    # SD for cross-regression coefficients

# Innovation parameters
sigma_innovation <- 1.0    # innovation SD
error_correlation <- 0.3   # innovation correlation

# ==============================================================================
# 构建完整参数向量
# ==============================================================================
# 参数顺序: [Burst1_Y1, Burst1_Y2, Burst2_Y1, Burst2_Y2, ..., AR1, AR2, CR12, CR21]
n_intercept_params <- nBurst * ydim  # intercept参数数量
n_var_params <- 4                     # VAR参数数量 (AR1, AR2, CR12, CR21)
n_total_params <- n_intercept_params + n_var_params

# 构建均值向量
param_means <- c()
for(b in 1:nBurst) {
  param_means <- c(param_means, burst_params$mu[[b]])
}
param_means <- c(param_means, AR_mean1, AR_mean2, CR12_mean, CR21_mean)

# 构建标准差向量
param_sds <- c()
for(b in 1:nBurst) {
  param_sds <- c(param_sds, burst_params$sigma[[b]])
}
param_sds <- c(param_sds, sigma_AR[1], sigma_AR[2], sigma_CR[1], sigma_CR[2])

######################## Correlation Structure ########################
# 构建相关矩阵 (当前设为单位矩阵，即无相关)
R_full <- diag(n_total_params)

# 转换为协方差矩阵
Sigma_full <- diag(param_sds) %*% R_full %*% diag(param_sds)

# Innovation covariance matrix
Sigma_innovation <- matrix(
  c(sigma_innovation^2, 
    error_correlation * sigma_innovation^2,
    error_correlation * sigma_innovation^2, 
    sigma_innovation^2), 
  2, 2
)

# ==============================================================================
# 打印配置信息
# ==============================================================================
print("========================================")
print(paste0("配置: ", nBurst, " bursts, ", nT_per_burst, " time points per burst"))
print(paste0("总时间点: ", nBurst * nT_per_burst))
print(paste0("总参数数: ", n_total_params, " (", n_intercept_params, " intercepts + ", n_var_params, " VAR)"))
print("Intercept均值:")
for(b in 1:nBurst) {
  print(paste0("  Burst ", b, ": Y1=", burst_params$mu[[b]][1], ", Y2=", burst_params$mu[[b]][2]))
}
print(paste0("VAR参数: AR1=", AR_mean1, ", AR2=", AR_mean2, ", CR12=", CR12_mean, ", CR21=", CR21_mean))
print("========================================")

########################################################################
######################## MC Experiments begins ########################
########################################################################

for(nT in nT_list) {
  for(nP in nP_list) {
    
    rs <- 1:N_repl
    
    # 计算每个burst的边界
    burst_boundaries <- seq(0, nT, by = nT_per_burst)
    
    Total_t_begin <- proc.time()
    
    for(k in 1:N_repl) {
      r <- rs[k]
      set.seed(r)
      
      dat_name <- paste0("Data_", C, "_nT", nT, "_nP", nP, "_r", r)
      dat_filename <- paste0(dat_name, ".csv")
      print(paste0("生成: ", dat_filename))
      
      ######################## Data Generation ########################
      # Storage arrays
      Y <- array(NA, c(nP, nT, ydim))
      Y_trend <- array(NA, c(nP, nT, ydim))
      Y_dev <- array(NA, c(nP, nT, ydim))
      person_params <- matrix(NA, nP, n_total_params)
      
      # Generate person-specific parameters
      for(i in 1:nP) {
        person_params[i,] <- rmvnorm(1, param_means, Sigma_full)
      }
      
      # Extract VAR parameters
      AR1 <- person_params[, n_intercept_params + 1]
      AR2 <- person_params[, n_intercept_params + 2]
      CR12 <- person_params[, n_intercept_params + 3]
      CR21 <- person_params[, n_intercept_params + 4]
      
      # Time sequence
      time <- seq(0, nT - 1, length.out = nT)
      
      # Generate data for each person
      for(i in 1:nP) {
        
        # Loop through each burst
        for(b in 1:nBurst) {
          # Burst起止时间点
          T_start <- burst_boundaries[b] + 1
          T_end <- burst_boundaries[b + 1]
          
          # 该burst的intercept参数索引
          intercept_idx <- ((b - 1) * ydim + 1):(b * ydim)
          
          # Burst的第一个时间点
          Y_trend[i, T_start, 1:ydim] <- person_params[i, intercept_idx]
          Y_dev[i, T_start, ] <- rmvnorm(1, mean = c(0, 0), sigma = Sigma_innovation)
          Y[i, T_start, ] <- Y_dev[i, T_start, ] + Y_trend[i, T_start, ]
          
          # Burst内的后续时间点 (如果有的话)
          if(T_end > T_start) {
            for(t in (T_start + 1):T_end) {
              # Trend保持不变
              Y_trend[i, t, 1:ydim] <- person_params[i, intercept_idx]
              
              # VAR process
              mean_dev <- c(
                AR1[i] * Y_dev[i, t-1, 1] + CR21[i] * Y_dev[i, t-1, 2],  # Y1
                AR2[i] * Y_dev[i, t-1, 2] + CR12[i] * Y_dev[i, t-1, 1]   # Y2
              )
              
              # Add innovations
              Y_dev[i, t, ] <- rmvnorm(1, mean = mean_dev, sigma = Sigma_innovation)
              
              # Combine trend and deviation
              Y[i, t, ] <- Y_trend[i, t, ] + Y_dev[i, t, ]
            }
          }
        }  # End burst loop
      }  # End person loop
      
      ######################## Save Data ########################
      # Reshape to long format
      dat_long <- matrix(NA, nrow = nP * nT, ncol = 4)
      colnames(dat_long) <- c("id", "time", "y1", "y2")
      
      for(i in 1:nP) {
        for(t in 1:nT) {
          row_idx <- (i - 1) * nT + t
          dat_long[row_idx, "id"] <- i
          dat_long[row_idx, "time"] <- time[t]
          dat_long[row_idx, "y1"] <- Y[i, t, 1]
          dat_long[row_idx, "y2"] <- Y[i, t, 2]
        }
      }
      
      dat_long <- as.data.frame(dat_long)
      
      # Save data
      write.csv(dat_long, paste0(work_path, "/", dat_name, ".csv"), row.names = FALSE)
      
      print(paste0("完成复制 ", k, " for N=", nP, ", T=", nT))
      
    }  # End replication loop
    
    Total_t_end <- proc.time()
    Total_t_elapsed <- Total_t_end - Total_t_begin
    print(paste0("条件 N=", nP, ", T=", nT, " 完成, 耗时: ", round(Total_t_elapsed[3], 2), "秒"))
    
  }  # End nP loop
}  # End nT loop

########################################################################
######################## MC Experiments ends ########################
########################################################################

# ==============================================================================
# 保存参数配置 (用于模型拟合时参考)
# ==============================================================================
config <- list(
  nBurst = nBurst,
  nT_per_burst = nT_per_burst,
  nT = nBurst * nT_per_burst,
  ydim = ydim,
  n_intercept_params = n_intercept_params,
  n_var_params = n_var_params,
  n_total_params = n_total_params,
  burst_boundaries = burst_boundaries,
  param_means = param_means,
  param_sds = param_sds,
  Sigma_innovation = Sigma_innovation
)

save(config, file = paste0(work_path, "/../config_", C, ".RData"))
print(paste0("配置已保存到: config_", C, ".RData"))

print("")
print("===== 数据生成全部完成 =====")
