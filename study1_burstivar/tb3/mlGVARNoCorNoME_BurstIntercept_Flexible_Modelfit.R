# ----------------------------------------------------------------------------
# PUBLIC-REPO COPY. Original: IP_1b/Study1_Tb3_Corrected/mlGVARNoCorNoME_BurstIntercept_Flexible_Modelfit.R
# Modified for this repository: the absolute base_path line was replaced with getwd()
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
# ==============================================================================
# mlGVARNoCorNoME_BurstIntercept_Flexible_Modelfit.R
# 功能：灵活的模型拟合脚本，支持任意数量的burst
# 特点：只需修改配置区域即可适配不同burst数
# ==============================================================================
rm(list=ls())
library(parallel)

# ==============================================================================
# ★★★ USER CONFIGURATION - 在此修改配置 ★★★
# ==============================================================================
# Burst配置
nBurst <- 3           # burst数量
nT_per_burst <- 3     # 每个burst内的时间点数

# 路径配置
base_path <- getwd()  # REPO EDIT: run from this directory (was D:/XXYDATAanalysis/IP_1b/Study1_Tb3_Corrected)

# 实验条件
nT_list <- c(nBurst * nT_per_burst)  # 总时间点
nP_list <- c(100, 500)               # 被试数
N_repl <- 100                        # Monte Carlo复制数

# 并行配置
num_cores <- 12

# ==============================================================================
# 自动计算的参数 (无需修改)
# ==============================================================================
ydim <- 2  # 维度数 (固定)
n_intercept_params <- nBurst * ydim
n_var_params <- 4
n_total_params <- n_intercept_params + n_var_params

# 条件名称
C_data <- paste0("mlGVARNoCorNoME_BurstIntercept_", nBurst, "Burst_", nT_per_burst, "T")
C_fit <- paste0("mlGVARNoCorNoME_BurstIntercept_", nBurst, "Burst")

# ==============================================================================
# 辅助函数：找出缺失的replications
# ==============================================================================
find_missing <- function(base_path, C_data, nT, nP, N_repl) {
  result_path <- file.path(base_path, "result")
  
  missing <- c()
  for(r in 1:N_repl) {
    result_file <- file.path(result_path, paste0(C_data, "_resulttable_nT", nT, "_nP", nP, "_r", r, ".csv"))
    if(!file.exists(result_file)) {
      missing <- c(missing, r)
    }
  }
  return(missing)
}

# ==============================================================================
# 辅助函数：生成JAGS数据的时间边界
# ==============================================================================
get_jags_time_boundaries <- function(nBurst, nT_per_burst) {
  boundaries <- list()
  for(b in 1:nBurst) {
    boundaries[[paste0("T", b)]] <- as.integer(b * nT_per_burst)
  }
  return(boundaries)
}

# ==============================================================================
# Worker函数：运行单个replication
# ==============================================================================
run_one_replication <- function(r, base_path, C_data, C_fit, nT, nP, nBurst, nT_per_burst) {
  
  library(rjags)
  library(coda)
  library(mvtnorm)
  
  # 路径设置
  data_path <- file.path(base_path, "data")
  result_path <- file.path(base_path, "result")
  if(!dir.exists(result_path)) dir.create(result_path, recursive = TRUE)
  
  # 加载辅助函数
  source(file.path(base_path, "posteriorSummaryStats.R"))
  
  # 文件路径
  model_file <- file.path(base_path, paste0(C_fit, ".txt"))
  data_file <- file.path(data_path, paste0("Data_", C_data, "_nT", nT, "_nP", nP, "_r", r, ".csv"))
  resulttable_file <- file.path(result_path, paste0(C_data, "_resulttable_nT", nT, "_nP", nP, "_r", r, ".csv"))
  codasamples_path <- file.path(result_path, paste0(C_data, "_codaSamples_nT", nT, "_nP", nP, "_r", r, ".RData"))
  
  # 检查数据文件
  if(!file.exists(data_file)) {
    return(paste0("Error: Data not found r=", r))
  }
  
  # 检查模型文件
  if(!file.exists(model_file)) {
    return(paste0("Error: Model file not found: ", model_file))
  }
  
  tryCatch({
    # 读取数据
    dat_long <- read.csv(data_file)
    Y <- array(NA, dim=c(nP, nT, 2))
    Y[,,1] <- matrix(dat_long$y1, nrow=nP, byrow=TRUE)
    Y[,,2] <- matrix(dat_long$y2, nrow=nP, byrow=TRUE)
    
    # 构建JAGS数据
    jags_data <- list(
      Y = Y, 
      P = nP, 
      W0 = diag(2), 
      df_res = 3
    )
    
    # 添加时间边界
    for(b in 1:nBurst) {
      jags_data[[paste0("T", b)]] <- as.integer(b * nT_per_burst)
    }
    
    # 监测参数
    parameters <- c("Level2Mean", "Level2Sigma", "sigma_innovation")
    
    # 初始值
    inits <- list(
      list(.RNG.name="base::Mersenne-Twister", .RNG.seed=r),
      list(.RNG.name="base::Mersenne-Twister", .RNG.seed=r+500)
    )
    
    # 运行JAGS
    jagsModel <- jags.model(model_file, data=jags_data, inits=inits, 
                            n.chains=2, n.adapt=5000, quiet=TRUE)
    update(jagsModel, n.iter=5000, progress.bar="none")
    codaSamples <- coda.samples(jagsModel, variable.names=parameters, 
                                n.iter=20000, progress.bar="none")
    
    # 保存结果
    resulttable <- summarizePost(codaSamples)
    write.csv(resulttable, resulttable_file)
    save(codaSamples, file = codasamples_path)
    
    return(paste0("Success: r=", r))
    
  }, error = function(e) {
    return(paste0("Fail: r=", r, " - ", e$message))
  })
}

# ==============================================================================
# 打印配置信息
# ==============================================================================
print("========================================")
print("模型拟合配置")
print("========================================")
print(paste0("Burst数量: ", nBurst))
print(paste0("每Burst时间点: ", nT_per_burst))
print(paste0("总时间点: ", nBurst * nT_per_burst))
print(paste0("参数总数: ", n_total_params, " (", n_intercept_params, " intercepts + ", n_var_params, " VAR)"))
print(paste0("数据条件名: ", C_data))
print(paste0("模型文件: ", C_fit, ".txt"))
print(paste0("并行核心数: ", num_cores))
print("========================================")

# ==============================================================================
# 检查缺失状态
# ==============================================================================
print("")
print("========== 检查缺失状态 ==========")
for(nT in nT_list) {
  for(nP in nP_list) {
    missing <- find_missing(base_path, C_data, nT, nP, N_repl)
    print(paste0("nT=", nT, ", nP=", nP, ": 缺失 ", length(missing), " 个"))
    if(length(missing) > 0 && length(missing) <= 30) {
      print(paste0("  缺失的r: ", paste(missing, collapse=",")))
    }
  }
}
print("===================================")

# ==============================================================================
# 运行缺失的任务
# ==============================================================================
for(nT in nT_list) {
  for(nP in nP_list) {
    
    repl_list <- find_missing(base_path, C_data, nT, nP, N_repl)
    
    if(length(repl_list) == 0) {
      print(paste0("nT=", nT, ", nP=", nP, ": 已全部完成，跳过"))
      next
    }
    
    print("")
    print(paste0("######## nT=", nT, ", nP=", nP, " ########"))
    print(paste0("运行 ", length(repl_list), " 个replications"))
    
    # 创建并行集群
    cl <- makeCluster(num_cores)
    clusterExport(cl, varlist = c("base_path", "C_data", "C_fit", "nT", "nP", 
                                   "nBurst", "nT_per_burst", "run_one_replication"))
    
    t_start <- proc.time()
    results <- parLapplyLB(cl, repl_list, function(r) {
      run_one_replication(r, base_path, C_data, C_fit, nT, nP, nBurst, nT_per_burst)
    })
    stopCluster(cl)
    
    t_elapsed <- proc.time() - t_start
    
    # 汇总结果
    success_count <- sum(grepl("Success", results))
    fail_count <- length(results) - success_count
    print(paste0("完成: ", success_count, " 成功, ", fail_count, " 失败"))
    print(paste0("耗时: ", round(t_elapsed[3]/60, 2), " 分钟"))
    
    if(fail_count > 0) {
      print("失败详情:")
      print(unlist(results[!grepl("Success", results)]))
    }
  }
}

print("")
print("===== 全部完成 =====")

# ==============================================================================
# 真值表 (用于后续分析参考)
# ==============================================================================
# 打印真值信息，方便核对
print("")
print("========== 真值参数参考 ==========")
print("Intercept参数:")
for(b in 1:nBurst) {
  idx_y1 <- (b - 1) * 2 + 1
  idx_y2 <- b * 2
  print(paste0("  Burst ", b, ": Level2Mean[", idx_y1, "] (Y1), Level2Mean[", idx_y2, "] (Y2)"))
}
print(paste0("VAR参数: Level2Mean[", n_intercept_params + 1, "]=AR1, Level2Mean[", n_intercept_params + 2, "]=AR2"))
print(paste0("         Level2Mean[", n_intercept_params + 3, "]=CR12, Level2Mean[", n_intercept_params + 4, "]=CR21"))
print("===================================")
