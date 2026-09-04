# --------------------------------------------------------------------
# BurstiVAR public repository -- Study 2 (GoBurstiVAR).
# Original file: paper_simulation_archive/Study2_Gompertz/code/Generate_Summ_GompertzBurst_BurstIntercept.R
# Modification for this repository: base_path set to "." (was an
# absolute local path); run this script from this directory. It reads
# per-replication result tables from ./result and writes the MCfile and
# MCfileSumm CSVs there. The nT_list/nP_list headers are preserved as
# last edited; set them per cell (see README).
# All other content is identical to the version that produced the
# published results.
# --------------------------------------------------------------------

########################################################################
# Script: Generate_Summ_GompertzBurst_BurstIntercept.R
# Description: Aggregates results from fitting BurstIntercept model 
#              to GompertzBurst-generated data.
#              - For mismatched parameters (burst intercepts, Level2Sigma for intercepts),
#                true values are set to NA but other statistics are still reported.
#              - For matched parameters (AR, CR, sigma_innovation), 
#                full statistics including coverage and power are computed.
########################################################################

rm(list=ls())

required_packages <- c("dplyr", "tidyverse", "data.table")
for(pkg in required_packages){
  if(!require(pkg, character.only = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

######################## 1. Environment & Path Setting ########################
# Roar Cluster path
# base_path <- "/storage/work/xjx5093/BurstIntercept/mlGVARNoCorNoME_GompertzBurst"

# Local PC path
base_path <- "."  # repo: run this script from this directory

result_path <- paste0(base_path, "/result")

print(paste0("Result path: ", result_path))
if(!dir.exists(result_path)) stop("Result directory not found! Please check 'base_path'.")
setwd(result_path)

######################## 2. Helper Functions ########################
Mean = function(x){ mean(x, na.rm=TRUE) }
SD = function(x){ sd(x, na.rm=TRUE) }
# Relative Bias (avoid division by zero, return NA if true value is NA or near zero)
relBias = function(x, truex){ 
  if(is.na(truex) || abs(truex) < 1e-6) return(NA) 
  mean((x - truex) / truex, na.rm=T) 
}
# RMSE (return NA if true value is NA)
RMSE = function(x, truex){ 
  if(any(is.na(truex))) return(NA)
  sqrt(Mean((x - truex)^2)) 
}

# Gompertz function for computing expected burst intercepts
gompertz <- function(b, theta1, theta2, theta3) {
  theta1 * exp(-theta2 * exp(-b * theta3))
}

# --- Robust Value Extractor ---
find_col_index <- function(col_names, keywords) {
  for (key in keywords) {
    idx <- grep(key, col_names, ignore.case = TRUE)
    if (length(idx) > 0) return(idx[1])
  }
  return(NA)
}

######################## 3. Simulation Conditions ########################
nT_list <- c(25)       
nP_list <- c(500)      
N_repl <- 100          

# Data generation model
C_DataGen <- "mlGVARNoCorNoME_GompertzBurst"
# Fitting model (result file naming)
C_Fit <- "mlGVARNoCorNoME_BurstIntercept_5Burst"

######################## 4. True Parameters from Gompertz Data Generation ########################
# Gompertz parameters for y1
theta10_1 = 10    # asymptote
theta20_1 = 3     # displacement  
theta30_1 = 1.0   # growth rate

# Gompertz parameters for y2
theta10_2 = 8     # asymptote
theta20_2 = 2.5   # displacement
theta30_2 = 0.8   # growth rate

# VAR parameters (these have direct correspondence)
mu_AR1 = 0.3      # AR for y1
mu_AR2 = 0.2      # AR for y2
mu_CR12 = -0.15   # Cross effect: y1 -> y2 (corresponds to CR[pp,1] in JAGS model)
mu_CR21 = -0.1    # Cross effect: y2 -> y1 (corresponds to CR[pp,2] in JAGS model)

# SD for VAR parameters (these have direct correspondence)
sigma_AR = 0.1    # SD for AR coefficients
sigma_CR = 0.1    # SD for cross-regression coefficients

# Innovation parameters
true_sigma_val = 1.0^2           # innovation variance
true_cov_val   = 0.3 * 1.0^2     # innovation covariance

# Calculate "expected" burst intercepts using Gompertz function
# These are the population mean values at each burst
burst_intercepts_y1 <- sapply(1:5, function(b) gompertz(b, theta10_1, theta20_1, theta30_1))
burst_intercepts_y2 <- sapply(1:5, function(b) gompertz(b, theta10_2, theta20_2, theta30_2))

cat("\n========== Expected Burst Intercepts from Gompertz ==========\n")
cat("Y1 (Gompertz: theta1=10, theta2=3, theta3=1.0):\n")
for(b in 1:5) cat(paste0("  Burst ", b, ": ", round(burst_intercepts_y1[b], 4), "\n"))
cat("Y2 (Gompertz: theta1=8, theta2=2.5, theta3=0.8):\n")
for(b in 1:5) cat(paste0("  Burst ", b, ": ", round(burst_intercepts_y2[b], 4), "\n"))
cat("===============================================================\n\n")

######################## 5. Build True Values Map ########################
# BurstIntercept model parameters:
# Level2Mean[1-10]: Burst intercepts (5 bursts x 2 variables)
#   [1,2] = Burst1 (y1, y2)
#   [3,4] = Burst2 (y1, y2)
#   [5,6] = Burst3 (y1, y2)
#   [7,8] = Burst4 (y1, y2)
#   [9,10] = Burst5 (y1, y2)
# Level2Mean[11-14]: VAR parameters
#   [11] = AR1 (y1 autoregression)
#   [12] = AR2 (y2 autoregression)
#   [13] = CR1 (y1->y2, corresponds to CR12_mean)
#   [14] = CR2 (y2->y1, corresponds to CR21_mean)

True_Values_Map <- c(
  # Burst intercepts - NA (model mismatch: Gompertz vs constant intercepts)
  "Level2Mean[1]"  = NA,   # Burst1 y1
  "Level2Mean[2]"  = NA,   # Burst1 y2
  "Level2Mean[3]"  = NA,   # Burst2 y1
  "Level2Mean[4]"  = NA,   # Burst2 y2
  "Level2Mean[5]"  = NA,   # Burst3 y1
  "Level2Mean[6]"  = NA,   # Burst3 y2
  "Level2Mean[7]"  = NA,   # Burst4 y1
  "Level2Mean[8]"  = NA,   # Burst4 y2
  "Level2Mean[9]"  = NA,   # Burst5 y1
  "Level2Mean[10]" = NA,   # Burst5 y2
  
  # VAR parameters - direct correspondence
  "Level2Mean[11]" = mu_AR1,                   # AR1
  "Level2Mean[12]" = mu_AR2,                   # AR2
  "Level2Mean[13]" = mu_CR12,                  # CR: y1->y2
  "Level2Mean[14]" = mu_CR21,                  # CR: y2->y1
  
  # Level2Sigma for burst intercepts - NA (structure mismatch)
  # In Gompertz model, individual differences come from Gompertz parameters,
  # not directly from intercept random effects
  "Level2Sigma[1]"  = NA,   # Burst1 y1 SD
  "Level2Sigma[2]"  = NA,   # Burst1 y2 SD
  "Level2Sigma[3]"  = NA,   # Burst2 y1 SD
  "Level2Sigma[4]"  = NA,   # Burst2 y2 SD
  "Level2Sigma[5]"  = NA,   # Burst3 y1 SD
  "Level2Sigma[6]"  = NA,   # Burst3 y2 SD
  "Level2Sigma[7]"  = NA,   # Burst4 y1 SD
  "Level2Sigma[8]"  = NA,   # Burst4 y2 SD
  "Level2Sigma[9]"  = NA,   # Burst5 y1 SD
  "Level2Sigma[10]" = NA,   # Burst5 y2 SD
  
  # Level2Sigma for VAR parameters - direct correspondence
  "Level2Sigma[11]" = sigma_AR,                # AR1 SD
  "Level2Sigma[12]" = sigma_AR,                # AR2 SD
  "Level2Sigma[13]" = sigma_CR,                # CR12 SD
  "Level2Sigma[14]" = sigma_CR,                # CR21 SD
  
  # Innovation covariance matrix - direct correspondence
  "sigma_innovation[1,1]" = true_sigma_val,
  "sigma_innovation[2,2]" = true_sigma_val,
  "sigma_innovation[1,2]" = true_cov_val,
  "sigma_innovation[2,1]" = true_cov_val
)

# Parameter labels for better readability in output
Param_Labels <- c(
  "Level2Mean[1]"  = "Burst1_y1",
  "Level2Mean[2]"  = "Burst1_y2",
  "Level2Mean[3]"  = "Burst2_y1",
  "Level2Mean[4]"  = "Burst2_y2",
  "Level2Mean[5]"  = "Burst3_y1",
  "Level2Mean[6]"  = "Burst3_y2",
  "Level2Mean[7]"  = "Burst4_y1",
  "Level2Mean[8]"  = "Burst4_y2",
  "Level2Mean[9]"  = "Burst5_y1",
  "Level2Mean[10]" = "Burst5_y2",
  "Level2Mean[11]" = "AR1",
  "Level2Mean[12]" = "AR2",
  "Level2Mean[13]" = "CR_y1toy2",
  "Level2Mean[14]" = "CR_y2toy1",
  "Level2Sigma[1]"  = "SD_Burst1_y1",
  "Level2Sigma[2]"  = "SD_Burst1_y2",
  "Level2Sigma[3]"  = "SD_Burst2_y1",
  "Level2Sigma[4]"  = "SD_Burst2_y2",
  "Level2Sigma[5]"  = "SD_Burst3_y1",
  "Level2Sigma[6]"  = "SD_Burst3_y2",
  "Level2Sigma[7]"  = "SD_Burst4_y1",
  "Level2Sigma[8]"  = "SD_Burst4_y2",
  "Level2Sigma[9]"  = "SD_Burst5_y1",
  "Level2Sigma[10]" = "SD_Burst5_y2",
  "Level2Sigma[11]" = "SD_AR1",
  "Level2Sigma[12]" = "SD_AR2",
  "Level2Sigma[13]" = "SD_CR_y1toy2",
  "Level2Sigma[14]" = "SD_CR_y2toy1",
  "sigma_innovation[1,1]" = "InnoVar_y1",
  "sigma_innovation[2,2]" = "InnoVar_y2",
  "sigma_innovation[1,2]" = "InnoCov",
  "sigma_innovation[2,1]" = "InnoCov_sym"
)

######################## 6. Main Loop ########################

for(nT in nT_list){
  for(nP in nP_list){
    
    cat(paste0("\n========== Processing: nT=", nT, ", nP=", nP, " ==========\n"))
    
    target_params <- names(True_Values_Map)
    # MCfile columns: r, then for each parameter: Tru, Est, se, LL, UL, cFlag, pFlag
    summ_cols_base <- c("Tru", "Est", "se", "LL", "UL", "cFlag", "pFlag")
    all_colnames <- c("r", as.vector(outer(summ_cols_base, target_params, paste, sep="_")))
    
    MCfile <- matrix(NA, nrow = N_repl, ncol = length(all_colnames))
    colnames(MCfile) <- all_colnames
    
    # ---------------------------------------------------------
    # First Pass: Detect Column Names from the first valid file
    # ---------------------------------------------------------
    col_map <- list()
    first_file <- list.files(result_path, pattern = paste0(C_Fit, ".*_nT", nT, "_nP", nP, "_r1\\.csv$"), full.names = T)
    if(length(first_file) > 1) first_file <- first_file[1]
    if(length(first_file) == 0) first_file <- list.files(result_path, pattern = "\\.csv$", full.names = T)[1]
    
    if(length(first_file) == 1 && !is.na(first_file) && file.exists(first_file)){
      temp_df <- read.csv(first_file)
      cnames <- colnames(temp_df)
      
      col_map$mean <- find_col_index(cnames, c("Mean", "mean", "Estimate", "Est"))
      col_map$se   <- find_col_index(cnames, c("SD", "PSD", "se", "StdDev"))
      col_map$ll   <- find_col_index(cnames, c("2.5", "Low", "LL", "min"))
      col_map$ul   <- find_col_index(cnames, c("97.5", "High", "UL", "max"))
      
      cat("Detected Column Mapping (Index):\n")
      print(unlist(col_map))
      cat("Corresponding Names:\n")
      print(cnames[unlist(col_map)])
    } else {
      cat("Warning: No result files found for this condition. Skipping...\n")
      next
    }
    
    # ---------------------------------------------------------
    # Processing Loop
    # ---------------------------------------------------------
    valid_reps <- 0
    
    for(r in 1:N_repl){
      file_name <- paste0(C_Fit, "_resulttable_nT", nT, "_nP", nP, "_r", r, ".csv")
      full_path <- file.path(result_path, file_name)
      
      if(file.exists(full_path)){
        res_table <- read.csv(full_path, row.names = 1) 
        rownames(res_table) <- trimws(rownames(res_table))
        
        MCfile[r, "r"] <- r
        
        for(par in target_params){
          if(par %in% rownames(res_table)){
            tru_val <- True_Values_Map[par]  # Can be NA for mismatched parameters
            
            # Safe extraction function
            get_val <- function(df, row, col_idx){
              val <- df[row, col_idx - 1]
              if(length(val) == 0 || is.null(val)) return(NA)
              return(as.numeric(val))
            }
            
            est_val <- get_val(res_table, par, col_map$mean)
            se_val  <- get_val(res_table, par, col_map$se)
            ll_val  <- get_val(res_table, par, col_map$ll)
            ul_val  <- get_val(res_table, par, col_map$ul)
            
            # Fallback: direct name matching
            if(is.na(est_val)){
              cnames_curr <- colnames(res_table)
              est_val <- res_table[par, grep("mean|Mean|Est", cnames_curr)[1]]
              se_val  <- res_table[par, grep("SD|PSD|se", cnames_curr)[1]]
              ll_val  <- res_table[par, grep("2.5|Low", cnames_curr)[1]]
              ul_val  <- res_table[par, grep("97.5|High", cnames_curr)[1]]
            }
            
            # Ensure scalar
            if(length(est_val) != 1) est_val <- NA
            if(length(se_val) != 1) se_val <- NA
            if(length(ll_val) != 1) ll_val <- NA
            if(length(ul_val) != 1) ul_val <- NA
            
            # Coverage flag (only if true value is not NA)
            cov_flag <- NA
            if(!is.na(tru_val) & !is.na(ll_val) & !is.na(ul_val)){
              cov_flag <- ifelse(tru_val >= ll_val & tru_val <= ul_val, 1, 0)
            }
            
            # Power / Type I Error Flag
            # For parameters with NA true value, still compute whether CI excludes zero
            power_flag <- NA
            if(!is.na(ll_val) & !is.na(ul_val)){
              is_variance_param <- grepl("Level2Sigma", par) | 
                                   par == "sigma_innovation[1,1]" | 
                                   par == "sigma_innovation[2,2]"
              
              if(is_variance_param){
                power_flag <- ifelse(ll_val > 0, 1, 0)
              } else {
                power_flag <- ifelse(ll_val > 0 | ul_val < 0, 1, 0)
              }
            }
            
            MCfile[r, paste0("Tru_", par)] <- tru_val
            MCfile[r, paste0("Est_", par)] <- est_val
            MCfile[r, paste0("se_", par)]  <- se_val
            MCfile[r, paste0("LL_", par)]  <- ll_val
            MCfile[r, paste0("UL_", par)]  <- ul_val
            MCfile[r, paste0("cFlag_", par)] <- cov_flag
            MCfile[r, paste0("pFlag_", par)] <- power_flag
          }
        }
        valid_reps <- valid_reps + 1
      }
    }
    
    cat(paste0("Found ", valid_reps, " valid replication files.\n"))
    
    if(valid_reps == 0){
      cat("No valid files found. Skipping this condition.\n")
      next
    }
    
    # Save Raw MCfile
    MCfile_name <- paste0(C_Fit, "_on_", C_DataGen, "_MCfile_nT", nT, "_nP", nP, ".csv")
    write.csv(MCfile, MCfile_name, row.names = FALSE)
    cat(paste0("MCfile saved to: ", MCfile_name, "\n"))
    
    ######################## 7. Compute Summary Statistics ########################
    
    # Columns: 
    # - TruePar: true parameter value (NA for mismatched)
    # - MeanPar_hat: mean of estimates
    # - RMSE, rBias: only computed if TruePar is not NA
    # - SD: SD of estimates
    # - Mean_SE_hat: mean of posterior SDs
    # - RDSE: relative deviation of SE
    # - 95%CI_LL, 95%CI_UL: quantiles of estimates
    # - Coverage: only computed if TruePar is not NA
    # - Power: for parameters with true value != 0
    # - TypeIError: for parameters with true value == 0
    # - DetectRate: proportion of CIs excluding zero (for all parameters)
    # - MissingPer: proportion of missing estimates
    
    summ_stat_cols <- c("Label", "TruePar", "MeanPar_hat", "RMSE", "rBias", "SD", 
                        "Mean_SE_hat", "RDSE", "95%CI_LL", "95%CI_UL", "Coverage", 
                        "Power", "TypeIError", "DetectRate", "MissingPer")
    
    MCfile_summ <- matrix(NA, nrow = length(target_params), ncol = length(summ_stat_cols))
    rownames(MCfile_summ) <- target_params
    colnames(MCfile_summ) <- summ_stat_cols
    
    MC_df <- as.data.frame(MCfile)
    
    for(par in target_params){
      est_vec <- as.numeric(MC_df[[paste0("Est_", par)]])
      tru_vec <- as.numeric(MC_df[[paste0("Tru_", par)]])
      se_vec  <- as.numeric(MC_df[[paste0("se_", par)]])
      flag_vec <- as.numeric(MC_df[[paste0("cFlag_", par)]])
      pflag_vec <- as.numeric(MC_df[[paste0("pFlag_", par)]])
      
      # Skip if all missing
      if(all(is.na(est_vec))) next
      
      # Get true parameter value (can be NA)
      true_par_val <- True_Values_Map[par]
      
      # Label
      MCfile_summ[par, "Label"] <- Param_Labels[par]
      
      # True value
      MCfile_summ[par, "TruePar"] <- true_par_val
      
      # Basic statistics (always computable)
      MCfile_summ[par, "MeanPar_hat"] <- Mean(est_vec)
      MCfile_summ[par, "SD"]          <- SD(est_vec) 
      MCfile_summ[par, "Mean_SE_hat"] <- Mean(se_vec)
      
      # RMSE and rBias - only if true value is not NA
      if(!is.na(true_par_val)){
        MCfile_summ[par, "RMSE"]  <- RMSE(est_vec, rep(true_par_val, length(est_vec)))
        MCfile_summ[par, "rBias"] <- relBias(Mean(est_vec), true_par_val)
      }
      
      # RDSE
      sd_est <- SD(est_vec)
      if(!is.na(sd_est) && sd_est > 0){
        MCfile_summ[par, "RDSE"] <- (Mean(se_vec) - sd_est) / sd_est
      }
      
      # CI of estimates
      MCfile_summ[par, "95%CI_LL"] <- quantile(est_vec, probs = 0.025, na.rm = T)
      MCfile_summ[par, "95%CI_UL"] <- quantile(est_vec, probs = 0.975, na.rm = T)
      
      # Coverage - only if true value is not NA
      if(!is.na(true_par_val)){
        valid_flags <- sum(!is.na(flag_vec))
        if(valid_flags > 0){
          MCfile_summ[par, "Coverage"] <- sum(flag_vec, na.rm = T) / valid_flags
        }
      }
      
      # Detection rate (always computable)
      valid_pflags <- sum(!is.na(pflag_vec))
      if(valid_pflags > 0){
        detection_rate <- sum(pflag_vec, na.rm = T) / valid_pflags
        MCfile_summ[par, "DetectRate"] <- detection_rate
        
        # Assign to Power or TypeIError based on true value
        if(!is.na(true_par_val)){
          if(abs(true_par_val) < 1e-6){
            # True value is essentially 0 -> Type I Error
            MCfile_summ[par, "TypeIError"] <- detection_rate
          } else {
            # True value is not 0 -> Power
            MCfile_summ[par, "Power"] <- detection_rate
          }
        }
      }
      
      # Missing proportion
      MCfile_summ[par, "MissingPer"] <- sum(is.na(est_vec)) / N_repl
    }
    
    # Convert to data frame for better display
    MCfile_summ_df <- as.data.frame(MCfile_summ)
    
    # Save summary
    file_summ_name <- paste0(C_Fit, "_on_", C_DataGen, "_MCfileSumm_nT", nT, "_nP", nP, ".csv")
    write.csv(MCfile_summ_df, file_summ_name)
    cat(paste0("Summary saved to: ", file_summ_name, "\n"))
    
    # Print summary for key parameters
    cat("\n--- Key Parameter Summary ---\n")
    key_params <- c("Level2Mean[11]", "Level2Mean[12]", "Level2Mean[13]", "Level2Mean[14]",
                    "sigma_innovation[1,1]", "sigma_innovation[1,2]")
    for(kp in key_params){
      if(kp %in% rownames(MCfile_summ_df)){
        cat(paste0(Param_Labels[kp], " (", kp, "): ",
                   "True=", round(as.numeric(MCfile_summ_df[kp, "TruePar"]), 4), ", ",
                   "Est=", round(as.numeric(MCfile_summ_df[kp, "MeanPar_hat"]), 4), ", ",
                   "Coverage=", round(as.numeric(MCfile_summ_df[kp, "Coverage"]), 3), "\n"))
      }
    }
  }
}

cat("\n========== Analysis Completed Successfully! ==========\n")
