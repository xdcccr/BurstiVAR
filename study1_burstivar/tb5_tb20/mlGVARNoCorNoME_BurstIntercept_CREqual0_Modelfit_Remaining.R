# ----------------------------------------------------------------------------
# PUBLIC-REPO COPY. Original: IP_1b/MLGVAR2601/mlGVARNoCorNoME_BurstIntercept_CREqual0/mlGVARNoCorNoME_BurstIntercept_CREqual0_Modelfit_Remaining.R
# Modified for this repository: the absolute base_path line was replaced with getwd()
# so the script runs from this directory (run it with the working directory
# set to study1_burstivar/tb5_tb20). All other content is identical to the
# version that produced the published results.
# ----------------------------------------------------------------------------
# ==============================================================================
# mlGVARNoCorNoME_BurstIntercept_Modelfit_Remaining.R
# 功能：只运行缺失的 replications
# ==============================================================================
rm(list=ls())
library(parallel)

# --- 1. 配置区域 ---
base_path <- getwd()  # REPO EDIT: run from this directory (was D:/XXYDATAanalysis/IP_1b/MLGVAR2601/mlGVARNoCorNoME_BurstIntercept_CREqual0)

nT_list <- c(15,60)
nP_list <- c(100,500)
N_repl <- 100

# ★ 减少核心数以避免内存不足 ★
num_cores <- 4  

# --- 2. 找出缺失的 replications ---
find_missing <- function(base_path, nT, nP, N_repl) {
  C <- "mlGVARNoCorNoME_BurstIntercept_CREqual0"
  result_path <- file.path(base_path, "result")
  
  missing <- c()
  for(r in 1:N_repl) {
    result_file <- file.path(result_path, paste0(C, "_resulttable_nT", nT, "_nP", nP, "_r", r, ".csv"))
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

# --- 4. Worker 函数 (与原脚本相同) ---
run_one_replication <- function(r, base_path, nT, nP) {
  
  library(rjags)
  library(coda)
  library(mvtnorm)
  
  C_data <- "mlGVARNoCorNoME_BurstIntercept_CREqual0"
  C_fit <- "mlGVARNoCorNoME_BurstIntercept"
  data_path <- file.path(base_path, "data")
  result_path <- file.path(base_path, "result")
  if(!dir.exists(result_path)) dir.create(result_path, recursive = TRUE)
  
  source(file.path(base_path, "posteriorSummaryStats.R"))
  
  model_file <- file.path(base_path, paste0(C_fit, ".txt"))
  data_file <- file.path(data_path, paste0("Data_", C_data, "_nT", nT, "_nP", nP, "_r", r, ".csv"))
  resulttable_file <- file.path(result_path, paste0(C_data, "_resulttable_nT", nT, "_nP", nP, "_r", r, ".csv"))
  codasamples_path <- file.path(result_path, paste0(C_data, "_codaSamples_nT", nT, "_nP", nP, "_r", r, ".RData"))
  
  if(!file.exists(data_file)) {
    return(paste0("Error: Data not found r=", r))
  }
  
  tryCatch({
    dat_long <- read.csv(data_file)
    Y <- array(NA, dim=c(nP, nT, 2))
    Y[,,1] <- matrix(dat_long$y1, nrow=nP, byrow=TRUE)
    Y[,,2] <- matrix(dat_long$y2, nrow=nP, byrow=TRUE)
    
    jags_data <- list(
      Y = Y, P = nP, W0 = diag(2), df_res = 3,
      T1 = as.integer(nT / 3),    
      T2 = as.integer(nT * 2 / 3), 
      T3 = as.integer(nT)  
    )
    
    parameters <- c("Level2Mean", "Level2Sigma", "sigma_innovation")
    
    inits <- list(
      list(.RNG.name="base::Mersenne-Twister", .RNG.seed=r),
      list(.RNG.name="base::Mersenne-Twister", .RNG.seed=r+500)
    )
    
    jagsModel <- jags.model(model_file, data=jags_data, inits=inits, 
                            n.chains=2, n.adapt=5000, quiet=TRUE)
    update(jagsModel, n.iter=5000, progress.bar="none")
    codaSamples <- coda.samples(jagsModel, variable.names=parameters, 
                                n.iter=20000, progress.bar="none")
    
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
