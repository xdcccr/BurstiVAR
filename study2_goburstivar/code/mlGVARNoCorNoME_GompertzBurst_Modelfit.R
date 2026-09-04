# --------------------------------------------------------------------
# BurstiVAR public repository -- Study 2 (GoBurstiVAR).
# Original file: paper_simulation_archive/Study2_Gompertz/code/mlGVARNoCorNoME_GompertzBurst_Modelfit.R
# Modification for this repository: base_path set to "." (was an
# absolute local path); run this script from this directory. It expects
# mlGVARNoCorNoME_GompertzBurst.txt, posteriorSummaryStats.R and ./data
# here, and writes to ./result. The nT_list/nP_list headers are
# preserved as last edited; set them per cell (see README).
# All other content is identical to the version that produced the
# published results.
# --------------------------------------------------------------------

# ==============================================================================
# mlGVARNoCorNoME_GompertzBurst_Modelfit.R
# 功能：用GompertzBurst模型（真实模型）拟合GompertzBurst生成的数据
# 说明：Gompertz函数基于burst number (1,2,3,4,5)，而非连续时间
# ==============================================================================
rm(list=ls())
library(parallel)

# --- 1. 配置区域 ---
# Roar Cluster路径
# base_path <- "/storage/work/xjx5093/BurstIntercept/mlGVARNoCorNoME_GompertzBurst"

# 本地PC路径
base_path <- "."  # repo: run this script from this directory

nT_list <- c(100)
nP_list <- c(500)
N_repl <- 100

# 核心数设置
num_cores <- 6  

# --- 2. 找出缺失的 replications ---
find_missing <- function(base_path, nT, nP, N_repl) {
  # 数据来源和拟合模型都是GompertzBurst（真实模型）
  C_result <- "mlGVARNoCorNoME_GompertzBurst"
  result_path <- file.path(base_path, "result")
  
  missing <- c()
  for(r in 1:N_repl) {
    result_file <- file.path(result_path, paste0(C_result, "_resulttable_nT", nT, "_nP", nP, "_r", r, ".csv"))
    if(!file.exists(result_file)) {
      missing <- c(missing, r)
    }
  }
  return(missing)
}

# --- 3. 检查状态 ---
print("========== 检查缺失状态 ==========")
for(nT in nT_list) {
  for(nP in nP_list) {
    missing <- find_missing(base_path, nT, nP, N_repl)
    print(paste0("nT=", nT, ", nP=", nP, ": 缺失 ", length(missing), " 个"))
    if(length(missing) > 0 && length(missing) <= 30) {
      print(paste0("  缺失的r: ", paste(missing, collapse=",")))
    }
  }
}
print("===================================")

# --- 4. Worker 函数 ---
run_one_replication <- function(r, base_path, nT, nP) {
  
  library(rjags)
  library(coda)
  library(mvtnorm)
  
  # 数据来源：GompertzBurst数据
  C_data <- "mlGVARNoCorNoME_GompertzBurst"
  # 拟合模型：GompertzBurst（真实模型）
  C_fit <- "mlGVARNoCorNoME_GompertzBurst"
  # 结果命名
  C_result <- "mlGVARNoCorNoME_GompertzBurst"
  
  data_path <- file.path(base_path, "data")
  result_path <- file.path(base_path, "result")
  if(!dir.exists(result_path)) dir.create(result_path, recursive = TRUE)
  
  source(file.path(base_path, "posteriorSummaryStats.R"))
  
  model_file <- file.path(base_path, paste0(C_fit, ".txt"))
  data_file <- file.path(data_path, paste0("Data_", C_data, "_nT", nT, "_nP", nP, "_r", r, ".csv"))
  resulttable_file <- file.path(result_path, paste0(C_result, "_resulttable_nT", nT, "_nP", nP, "_r", r, ".csv"))
  codasamples_path <- file.path(result_path, paste0(C_result, "_codaSamples_nT", nT, "_nP", nP, "_r", r, ".RData"))
  
  if(!file.exists(data_file)) {
    return(paste0("Error: Data not found r=", r))
  }
  
  tryCatch({
    # 读取数据
    dat_long <- read.csv(data_file)
    Y <- array(NA, dim=c(nP, nT, 2))
    Y[,,1] <- matrix(dat_long$y1, nrow=nP, byrow=TRUE)
    Y[,,2] <- matrix(dat_long$y2, nrow=nP, byrow=TRUE)
    
    # 5 burst的时间边界
    nBurst <- 5
    T_per_burst <- as.integer(nT / nBurst)
    
    # 准备JAGS数据
    jags_data <- list(
      Y = Y, 
      P = nP, 
      W0 = diag(2), 
      df_res = 3,
      T1 = T_per_burst,           # End of Burst 1
      T2 = T_per_burst * 2,       # End of Burst 2
      T3 = T_per_burst * 3,       # End of Burst 3
      T4 = T_per_burst * 4,       # End of Burst 4
      T5 = nT                     # End of Burst 5
    )
    
    # 监测参数
    parameters <- c("Level2Mean", "Level2Sigma", "sigma_innovation")
    
    # 初始化
    inits <- list(
      list(.RNG.name="base::Mersenne-Twister", .RNG.seed=r),
      list(.RNG.name="base::Mersenne-Twister", .RNG.seed=r+500)
    )
    
    # MCMC设置
    n_adapt <- 5000
    n_burnin <- 5000
    n_iter <- 20000
    
    # 运行JAGS
    jagsModel <- jags.model(model_file, data=jags_data, inits=inits, 
                            n.chains=2, n.adapt=n_adapt, quiet=TRUE)
    update(jagsModel, n.iter=n_burnin, progress.bar="none")
    codaSamples <- coda.samples(jagsModel, variable.names=parameters, 
                                n.iter=n_iter, progress.bar="none")
    
    # 保存结果
    resulttable <- summarizePost(codaSamples)
    write.csv(resulttable, resulttable_file)
    save(codaSamples, file = codasamples_path)
    
    return(paste0("Success: r=", r))
    
  }, error = function(e) {
    return(paste0("Fail: r=", r, " - ", e$message))
  })
}

# --- 5. 只运行缺失的任务 ---
for(nT in nT_list) {
  for(nP in nP_list) {
    
    repl_list <- find_missing(base_path, nT, nP, N_repl)
    
    if(length(repl_list) == 0) {
      print(paste0("nT=", nT, ", nP=", nP, ": 已全部完成，跳过"))
      next
    }
    
    print("")
    print(paste0("######## nT=", nT, ", nP=", nP, " ########"))
    print(paste0("运行 ", length(repl_list), " 个缺失的 replications"))
    
    cl <- makeCluster(num_cores)
    clusterExport(cl, varlist = c("base_path", "nT", "nP", "run_one_replication"))
    
    t_start <- proc.time()
    results <- parLapplyLB(cl, repl_list, function(r) {
      run_one_replication(r, base_path, nT, nP)
    })
    stopCluster(cl)
    
    t_elapsed <- proc.time() - t_start
    
    # 汇总
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
