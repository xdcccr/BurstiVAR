# ------------------------------------------------------------------
# BurstiVAR public repository copy (Supplement D, two-burst check)
# Original path: paper_simulation_archive/SupplementD_TwoBurst/code/mlGVARNoCorNoME_BurstIntercept_2Burst_3T_DataGeneration.R
# Modified for this repository: the absolute work_path was replaced with "./data";
#   run this script with this folder as the R working directory.
# All other content is identical to the version that produced the
# published results.
# ------------------------------------------------------------------

######################## Set Environment ########################
rm(list=ls())
library(mvtnorm)
library(Matrix)
library(dplyr)

# ==============================================================================
# 2-Burst, 3-TimePoints 专用数据生成脚本
# 配置: nBurst=2, nT_per_burst=3, 总共6个时间点
# ==============================================================================

######################## USER CONFIGURATION ########################
# 工作路径设置
work_path <- "./data"  # run this script from this directory
print(paste0("work_path: ", work_path))

# 如果目录不存在则创建
if(!dir.exists(work_path)) {
  dir.create(work_path, recursive = TRUE)
}
setwd(work_path)

######################## Control Parameters ########################
C <- "mlGVARNoCorNoME_BurstIntercept_2Burst_3T"

# 固定配置
nBurst <- 2
nT_per_burst <- 3
nT <- nBurst * nT_per_burst  # = 6

# 实验条件
for(nP in c(100, 500)){  # 被试数
  
  ydim <- 2        # 维度数 (reading and math)
  N_repl <- 100    # Monte Carlo复制数
  rs <- 1:N_repl
  
  ######################## Model Parameters ########################
  # Burst 1 Intercept
  Intercept_mu_Burst1 <- c(0, 1)        # Y1=0, Y2=1
  Intercept_sigma_Burst1 <- c(1, 1)      # SD: 1, 1
  
  # Burst 2 Intercept
  Intercept_mu_Burst2 <- c(1, 1.5)       # Y1=1, Y2=1.5
  Intercept_sigma_Burst2 <- c(1.5, 2)    # SD: 1.5, 2
  
  # VAR Parameters
  AR_mean1 <- 0.3      # AR for reading
  AR_mean2 <- 0.2      # AR for math
  CR12_mean <- -0.15   # Cross effect: reading -> math
  CR21_mean <- 0       # Cross effect: math -> reading
  
  # Standard deviations for VAR parameters
  sigma_AR <- c(0.1, 0.1)
  sigma_CR <- c(0.1, 0.1)
  
  # Innovation parameters
  sigma_innovation <- 1.0
  error_correlation <- 0.3
  
  # 完整参数向量 (8个参数)
  # [Burst1_Y1, Burst1_Y2, Burst2_Y1, Burst2_Y2, AR1, AR2, CR12, CR21]
  param_means <- c(
    Intercept_mu_Burst1,   # 1-2
    Intercept_mu_Burst2,   # 3-4
    AR_mean1, AR_mean2,    # 5-6
    CR12_mean, CR21_mean   # 7-8
  )
  
  param_sds <- c(
    Intercept_sigma_Burst1,  # 1-2
    Intercept_sigma_Burst2,  # 3-4
    sigma_AR[1], sigma_AR[2],  # 5-6
    sigma_CR[1], sigma_CR[2]   # 7-8
  )
  
  ######################## Correlation Structure ########################
  # 8x8 单位矩阵 (无相关)
  R_full <- diag(8)
  
  # 协方差矩阵
  Sigma_full <- diag(param_sds) %*% R_full %*% diag(param_sds)
  
  # Innovation covariance matrix
  Sigma_innovation <- matrix(
    c(sigma_innovation^2, 
      error_correlation * sigma_innovation^2,
      error_correlation * sigma_innovation^2, 
      sigma_innovation^2), 
    2, 2
  )
  
  ########################################################################
  ######################## MC Experiments begins ########################
  ########################################################################
  
  # 时间边界
  T1 <- nT_per_burst       # = 3
  T2 <- nT_per_burst * 2   # = 6
  
  Total_t_begin <- proc.time()
  
  for(k in 1:N_repl){
    r <- rs[k]
    set.seed(r)
    
    dat_name <- paste0("Data_", C, "_nT", nT, "_nP", nP, "_r", r)
    dat_filename <- paste0(dat_name, ".csv")
    print(paste0("生成: ", dat_filename))
    
    ######################## Data Generation ########################
    # Storage arrays
    Y <- array(NA, c(nP, nT, 2))
    Y_trend <- array(NA, c(nP, nT, 2))
    Y_dev <- array(NA, c(nP, nT, 2))
    person_params <- matrix(NA, nP, 8)
    
    # Generate person-specific parameters
    for(i in 1:nP){
      person_params[i,] <- rmvnorm(1, param_means, Sigma_full)
    }
    
    # Extract parameters
    # Intercept parameters
    Intercept1_Burst1 <- person_params[,1]
    Intercept2_Burst1 <- person_params[,2]
    Intercept1_Burst2 <- person_params[,3]
    Intercept2_Burst2 <- person_params[,4]
    
    # VAR parameters
    AR1 <- person_params[,5]
    AR2 <- person_params[,6]
    CR12 <- person_params[,7]
    CR21 <- person_params[,8]
    
    time <- seq(0, nT - 1, length.out = nT) 
    
    # Generate data for each person
    for(i in 1:nP){
      
      # ==================== Burst 1 ====================
      # First time point
      Y_trend[i,1,1:2] <- person_params[i,1:2]
      Y_dev[i,1,] <- rmvnorm(1, mean=c(0,0), sigma=Sigma_innovation)
      Y[i,1,] <- Y_dev[i,1,] + Y_trend[i,1,]
      
      # t=2 to T1
      for(t in 2:T1){
        Y_trend[i,t,1:2] <- person_params[i,1:2]
        
        mean_dev <- c(
          AR1[i] * Y_dev[i,t-1,1] + CR21[i] * Y_dev[i,t-1,2],
          AR2[i] * Y_dev[i,t-1,2] + CR12[i] * Y_dev[i,t-1,1]
        )
        
        Y_dev[i,t,] <- rmvnorm(1, mean=mean_dev, sigma=Sigma_innovation)
        Y[i,t,] <- Y_trend[i,t,] + Y_dev[i,t,]
      }
      
      # ==================== Burst 2 ====================
      # First time point of Burst 2
      Y_trend[i,T1+1,1:2] <- person_params[i,3:4]
      Y_dev[i,T1+1,] <- rmvnorm(1, mean=c(0,0), sigma=Sigma_innovation)
      Y[i,T1+1,] <- Y_dev[i,T1+1,] + Y_trend[i,T1+1,]
      
      # t=T1+2 to T2
      for(t in (T1+2):T2){
        Y_trend[i,t,1:2] <- person_params[i,3:4]
        
        mean_dev <- c(
          AR1[i] * Y_dev[i,t-1,1] + CR21[i] * Y_dev[i,t-1,2],
          AR2[i] * Y_dev[i,t-1,2] + CR12[i] * Y_dev[i,t-1,1]
        )
        
        Y_dev[i,t,] <- rmvnorm(1, mean=mean_dev, sigma=Sigma_innovation)
        Y[i,t,] <- Y_trend[i,t,] + Y_dev[i,t,]
      }
      
    }  # End person loop
    
    ######################## Save Data ########################
    # Reshape to long format
    dat_long <- matrix(NA, nrow=nP*nT, ncol=4)
    colnames(dat_long) <- c("id", "time", "y1", "y2")
    
    for(i in 1:nP){
      for(t in 1:nT){
        row_idx <- (i-1)*nT + t
        dat_long[row_idx, "id"] <- i
        dat_long[row_idx, "time"] <- time[t]
        dat_long[row_idx, "y1"] <- Y[i,t,1]
        dat_long[row_idx, "y2"] <- Y[i,t,2]
      }
    }
    
    dat_long <- as.data.frame(dat_long)
    
    write.csv(dat_long, paste0(work_path, "/", dat_name, ".csv"), row.names=FALSE)
    
    print(paste0("完成复制 ", k, " for N=", nP, ", T=", nT))
    
  }  # End replication loop
  
  Total_t_end <- proc.time()
  Total_t_elapsed <- Total_t_end - Total_t_begin
  print(paste0("条件 N=", nP, ", T=", nT, " 完成, 耗时: ", round(Total_t_elapsed[3], 2), "秒"))
  
}  # End nP loop

########################################################################
######################## MC Experiments ends ########################
########################################################################

# ==============================================================================
# 保存真值参数 (用于后续结果分析)
# ==============================================================================
true_params <- data.frame(
  Parameter = c("Level2Mean[1]", "Level2Mean[2]", "Level2Mean[3]", "Level2Mean[4]",
                "Level2Mean[5]", "Level2Mean[6]", "Level2Mean[7]", "Level2Mean[8]",
                "Level2Sigma[1]", "Level2Sigma[2]", "Level2Sigma[3]", "Level2Sigma[4]",
                "Level2Sigma[5]", "Level2Sigma[6]", "Level2Sigma[7]", "Level2Sigma[8]",
                "sigma_innovation[1,1]", "sigma_innovation[2,2]", "sigma_innovation[1,2]"),
  Label = c("Intercept_Burst1_Y1", "Intercept_Burst1_Y2", 
            "Intercept_Burst2_Y1", "Intercept_Burst2_Y2",
            "AR1", "AR2", "CR12", "CR21",
            "SD_Intercept_Burst1_Y1", "SD_Intercept_Burst1_Y2",
            "SD_Intercept_Burst2_Y1", "SD_Intercept_Burst2_Y2",
            "SD_AR1", "SD_AR2", "SD_CR12", "SD_CR21",
            "Var_Innovation_Y1", "Var_Innovation_Y2", "Cov_Innovation"),
  TrueValue = c(Intercept_mu_Burst1[1], Intercept_mu_Burst1[2],
                Intercept_mu_Burst2[1], Intercept_mu_Burst2[2],
                AR_mean1, AR_mean2, CR12_mean, CR21_mean,
                Intercept_sigma_Burst1[1], Intercept_sigma_Burst1[2],
                Intercept_sigma_Burst2[1], Intercept_sigma_Burst2[2],
                sigma_AR[1], sigma_AR[2], sigma_CR[1], sigma_CR[2],
                sigma_innovation^2, sigma_innovation^2, error_correlation*sigma_innovation^2)
)

write.csv(true_params, paste0(work_path, "/../TrueParameters_", C, ".csv"), row.names = FALSE)
print("")
print("真值参数已保存")
print(true_params)

print("")
print("===== 数据生成全部完成 =====")
