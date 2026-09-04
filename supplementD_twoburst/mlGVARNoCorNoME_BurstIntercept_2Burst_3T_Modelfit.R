# ------------------------------------------------------------------
# BurstiVAR public repository copy (Supplement D, two-burst check)
# Original path: paper_simulation_archive/SupplementD_TwoBurst/code/mlGVARNoCorNoME_BurstIntercept_2Burst_3T_Modelfit.R
# Modified for this repository: the absolute base_path was replaced with ".";
#   run this script with this folder as the R working directory.
# All other content is identical to the version that produced the
# published results.
# ------------------------------------------------------------------

# ==============================================================================
# mlGVARNoCorNoME_BurstIntercept_2Burst_3T_Modelfit.R
# 功能：2-Burst, 3-TimePoints 专用模型拟合脚本
# ==============================================================================
rm(list=ls())
library(parallel)

# ==============================================================================
# USER CONFIGURATION
# ==============================================================================
base_path <- "."  # run this script from this directory

# 固定配置
nBurst <- 2
nT_per_burst <- 3
nT <- 6  # = nBurst * nT_per_burst

# 实验条件
nT_list <- c(6)
nP_list <- c(100, 500)
N_repl <- 100

# 并行配置
num_cores <- 2

# 条件名称
C_data <- "mlGVARNoCorNoME_BurstIntercept_2Burst_3T"
C_fit <- "mlGVARNoCorNoME_BurstIntercept_2Burst"

# ==============================================================================
# 辅助函数
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
# Worker函数
# ==============================================================================
run_one_replication <- function(r, base_path, C_data, C_fit, nT, nP) {
  
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
  
  # 检查文件
  if(!file.exists(data_file)) {
    return(paste0("Error: Data not found r=", r))
  }
  if(!file.exists(model_file)) {
    return(paste0("Error: Model file not found: ", model_file))
  }
  
  tryCatch({
    # 读取数据
    dat_long <- read.csv(data_file)
    Y <- array(NA, dim=c(nP, nT, 2))
    Y[,,1] <- matrix(dat_long$y1, nrow=nP, byrow=TRUE)
    Y[,,2] <- matrix(dat_long$y2, nrow=nP, byrow=TRUE)
    
    # JAGS数据
    # 注意: 2-burst模型需要 T1 和 T2
    jags_data <- list(
      Y = Y, 
      P = nP, 
      W0 = diag(2), 
      df_res = 3,
      T1 = 3L,   # Burst 1 结束点
      T2 = 6L    # Burst 2 结束点 (= 总时间点)
    )
    
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
# 打印配置
# ==============================================================================
print("========================================")
print("2-Burst, 3-TimePoints 模型拟合")
print("========================================")
print(paste0("总时间点: ", nT))
print(paste0("Burst边界: T1=3, T2=6"))
print(paste0("参数: 8个 (4 intercepts + 4 VAR)"))
print(paste0("数据条件: ", C_data))
print(paste0("模型文件: ", C_fit, ".txt"))
print("========================================")

# ==============================================================================
# 检查状态
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
# 运行
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
    
    cl <- makeCluster(num_cores)
    clusterExport(cl, varlist = c("base_path", "C_data", "C_fit", "nT", "nP", 
                                   "run_one_replication"))
    
    t_start <- proc.time()
    results <- parLapplyLB(cl, repl_list, function(r) {
      run_one_replication(r, base_path, C_data, C_fit, nT, nP)
    })
    stopCluster(cl)
    
    t_elapsed <- proc.time() - t_start
    
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
# 真值参考
# ==============================================================================
print("")
print("========== 真值参数参考 ==========")
print("Level2Mean[1] = Intercept_Burst1_Y1 = 0")
print("Level2Mean[2] = Intercept_Burst1_Y2 = 1")
print("Level2Mean[3] = Intercept_Burst2_Y1 = 1")
print("Level2Mean[4] = Intercept_Burst2_Y2 = 1.5")
print("Level2Mean[5] = AR1 = 0.3")
print("Level2Mean[6] = AR2 = 0.2")
print("Level2Mean[7] = CR12 = -0.15")
print("Level2Mean[8] = CR21 = 0")
print("===================================")
