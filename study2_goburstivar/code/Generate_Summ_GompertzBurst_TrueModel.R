# --------------------------------------------------------------------
# BurstiVAR public repository -- Study 2 (GoBurstiVAR).
# Original file: paper_simulation_archive/Study2_Gompertz/code/Generate_Summ_GompertzBurst_TrueModel.R
# Modification for this repository: base_path set to "." (was an
# absolute local path); run this script from this directory. It reads
# per-replication result tables from ./result and writes the MCfile and
# MCfileSumm CSVs there. The nT_list/nP_list headers are preserved as
# last edited; set them per cell (see README).
# All other content is identical to the version that produced the
# published results.
# --------------------------------------------------------------------

########################################################################
# Script: Generate_Summ_GompertzBurst_TrueModel.R
# Description: Aggregates results from fitting GompertzBurst model (TRUE model)
#              to GompertzBurst-generated data.
#              All parameters have corresponding true values.
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
# Relative Bias
relBias = function(x, truex){ 
  if(is.na(truex) || abs(truex) < 1e-6) return(NA) 
  mean((x - truex) / truex, na.rm=T) 
}
# RMSE
RMSE = function(x, truex){ 
  if(any(is.na(truex))) return(NA)
  sqrt(Mean((x - truex)^2)) 
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
nT_list <- c(25)#c(25, 100)       # Updated to match data generation conditions
nP_list <- c(500)      
N_repl <- 100          

# Model name (TRUE model fitting TRUE data)
C_Model <- "mlGVARNoCorNoME_GompertzBurst"

######################## 4. True Parameters from Data Generation ########################
# Gompertz parameters for y1
theta10_1 = 10    # asymptote
theta20_1 = 3     # displacement  
theta30_1 = 1.0   # growth rate
sigma1_1 = 2      # SD for asymptote
sigma2_1 = 0.5    # SD for displacement
sigma3_1 = 0.2    # SD for growth rate

# Gompertz parameters for y2
theta10_2 = 8     # asymptote
theta20_2 = 2.5   # displacement
theta30_2 = 0.8   # growth rate
sigma1_2 = 1.5    # SD for asymptote
sigma2_2 = 0.4    # SD for displacement
sigma3_2 = 0.15   # SD for growth rate

# VAR parameters
mu_AR1 = 0.3      # AR for y1
mu_AR2 = 0.2      # AR for y2
mu_CR12 = -0.15   # Cross effect: y1 -> y2
mu_CR21 = -0.1    # Cross effect: y2 -> y1
sigma_AR = 0.1    # SD for AR coefficients
sigma_CR = 0.1    # SD for cross-regression coefficients

# Innovation parameters
true_sigma_val = 1.0^2           # innovation variance
true_cov_val   = 0.3 * 1.0^2     # innovation covariance

######################## 5. Build True Values Map ########################
# GompertzBurst model parameters:
# Level2Mean[1-6]: Gompertz parameters (3 per variable)
#   [1] = theta1_y1 (asymptote)
#   [2] = theta2_y1 (displacement)
#   [3] = theta3_y1 (growth rate)
#   [4] = theta1_y2 (asymptote)
#   [5] = theta2_y2 (displacement)
#   [6] = theta3_y2 (growth rate)
# Level2Mean[7-10]: VAR parameters
#   [7] = AR1
#   [8] = AR2
#   [9] = CR y1->y2
#   [10] = CR y2->y1

True_Values_Map <- c(
  # Gompertz parameters - all have true values
  "Level2Mean[1]"  = theta10_1,    # y1 asymptote
  "Level2Mean[2]"  = theta20_1,    # y1 displacement
  "Level2Mean[3]"  = theta30_1,    # y1 growth rate
  "Level2Mean[4]"  = theta10_2,    # y2 asymptote
  "Level2Mean[5]"  = theta20_2,    # y2 displacement
  "Level2Mean[6]"  = theta30_2,    # y2 growth rate
  
  # VAR parameters
  "Level2Mean[7]"  = mu_AR1,       # AR1
  "Level2Mean[8]"  = mu_AR2,       # AR2
  "Level2Mean[9]"  = mu_CR12,      # CR: y1->y2
  "Level2Mean[10]" = mu_CR21,      # CR: y2->y1
  
  # Level2Sigma for Gompertz parameters
  "Level2Sigma[1]"  = sigma1_1,    # y1 asymptote SD
  "Level2Sigma[2]"  = sigma2_1,    # y1 displacement SD
  "Level2Sigma[3]"  = sigma3_1,    # y1 growth rate SD
  "Level2Sigma[4]"  = sigma1_2,    # y2 asymptote SD
  "Level2Sigma[5]"  = sigma2_2,    # y2 displacement SD
  "Level2Sigma[6]"  = sigma3_2,    # y2 growth rate SD
  
  # Level2Sigma for VAR parameters
  "Level2Sigma[7]"  = sigma_AR,    # AR1 SD
  "Level2Sigma[8]"  = sigma_AR,    # AR2 SD
  "Level2Sigma[9]"  = sigma_CR,    # CR12 SD
  "Level2Sigma[10]" = sigma_CR,    # CR21 SD
  
  # Innovation covariance matrix
  "sigma_innovation[1,1]" = true_sigma_val,
  "sigma_innovation[2,2]" = true_sigma_val,
  "sigma_innovation[1,2]" = true_cov_val,
  "sigma_innovation[2,1]" = true_cov_val
)

# Parameter labels for better readability
Param_Labels <- c(
  "Level2Mean[1]"  = "Theta1_y1_Asymptote",
  "Level2Mean[2]"  = "Theta2_y1_Displacement",
  "Level2Mean[3]"  = "Theta3_y1_GrowthRate",
  "Level2Mean[4]"  = "Theta1_y2_Asymptote",
  "Level2Mean[5]"  = "Theta2_y2_Displacement",
  "Level2Mean[6]"  = "Theta3_y2_GrowthRate",
  "Level2Mean[7]"  = "AR1",
  "Level2Mean[8]"  = "AR2",
  "Level2Mean[9]"  = "CR_y1toy2",
  "Level2Mean[10]" = "CR_y2toy1",
  "Level2Sigma[1]"  = "SD_Theta1_y1",
  "Level2Sigma[2]"  = "SD_Theta2_y1",
  "Level2Sigma[3]"  = "SD_Theta3_y1",
  "Level2Sigma[4]"  = "SD_Theta1_y2",
  "Level2Sigma[5]"  = "SD_Theta2_y2",
  "Level2Sigma[6]"  = "SD_Theta3_y2",
  "Level2Sigma[7]"  = "SD_AR1",
  "Level2Sigma[8]"  = "SD_AR2",
  "Level2Sigma[9]"  = "SD_CR_y1toy2",
  "Level2Sigma[10]" = "SD_CR_y2toy1",
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
    summ_cols_base <- c("Tru", "Est", "se", "LL", "UL", "cFlag", "pFlag")
    all_colnames <- c("r", as.vector(outer(summ_cols_base, target_params, paste, sep="_")))
    
    MCfile <- matrix(NA, nrow = N_repl, ncol = length(all_colnames))
    colnames(MCfile) <- all_colnames
    
    # ---------------------------------------------------------
    # First Pass: Detect Column Names
    # ---------------------------------------------------------
    col_map <- list()
    first_file <- list.files(result_path, pattern = paste0(C_Model, "_resulttable_nT", nT, "_nP", nP, "_r1\\.csv$"), full.names = T)
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
    } else {
      cat("Warning: No result files found for this condition. Skipping...\n")
      next
    }
    
    # ---------------------------------------------------------
    # Processing Loop
    # ---------------------------------------------------------
    valid_reps <- 0
    
    for(r in 1:N_repl){
      file_name <- paste0(C_Model, "_resulttable_nT", nT, "_nP", nP, "_r", r, ".csv")
      full_path <- file.path(result_path, file_name)
      
      if(file.exists(full_path)){
        res_table <- read.csv(full_path, row.names = 1) 
        rownames(res_table) <- trimws(rownames(res_table))
        
        MCfile[r, "r"] <- r
        
        for(par in target_params){
          if(par %in% rownames(res_table)){
            tru_val <- True_Values_Map[par]
            
            get_val <- function(df, row, col_idx){
              val <- df[row, col_idx - 1]
              if(length(val) == 0 || is.null(val)) return(NA)
              return(as.numeric(val))
            }
            
            est_val <- get_val(res_table, par, col_map$mean)
            se_val  <- get_val(res_table, par, col_map$se)
            ll_val  <- get_val(res_table, par, col_map$ll)
            ul_val  <- get_val(res_table, par, col_map$ul)
            
            # Fallback
            if(is.na(est_val)){
              cnames_curr <- colnames(res_table)
              est_val <- res_table[par, grep("mean|Mean|Est", cnames_curr)[1]]
              se_val  <- res_table[par, grep("SD|PSD|se", cnames_curr)[1]]
              ll_val  <- res_table[par, grep("2.5|Low", cnames_curr)[1]]
              ul_val  <- res_table[par, grep("97.5|High", cnames_curr)[1]]
            }
            
            if(length(est_val) != 1) est_val <- NA
            if(length(se_val) != 1) se_val <- NA
            if(length(ll_val) != 1) ll_val <- NA
            if(length(ul_val) != 1) ul_val <- NA
            
            # Coverage
            cov_flag <- NA
            if(!is.na(tru_val) & !is.na(ll_val) & !is.na(ul_val)){
              cov_flag <- ifelse(tru_val >= ll_val & tru_val <= ul_val, 1, 0)
            }
            
            # Power / Type I Error Flag
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
    
    # Save MCfile
    MCfile_name <- paste0(C_Model, "_TrueModel_MCfile_nT", nT, "_nP", nP, ".csv")
    write.csv(MCfile, MCfile_name, row.names = FALSE)
    cat(paste0("MCfile saved to: ", MCfile_name, "\n"))
    
    ######################## 7. Compute Summary Statistics ########################
    
    summ_stat_cols <- c("Label", "TruePar", "MeanPar_hat", "RMSE", "rBias", "SD", 
                        "Mean_SE_hat", "RDSE", "95%CI_LL", "95%CI_UL", "Coverage", 
                        "Power", "TypeIError", "MissingPer")
    
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
      
      if(all(is.na(est_vec))) next
      
      true_par_val <- True_Values_Map[par]
      
      MCfile_summ[par, "Label"] <- Param_Labels[par]
      MCfile_summ[par, "TruePar"] <- true_par_val
      MCfile_summ[par, "MeanPar_hat"] <- Mean(est_vec)
      MCfile_summ[par, "SD"] <- SD(est_vec) 
      MCfile_summ[par, "Mean_SE_hat"] <- Mean(se_vec)
      
      if(!is.na(true_par_val)){
        MCfile_summ[par, "RMSE"] <- RMSE(est_vec, rep(true_par_val, length(est_vec)))
        MCfile_summ[par, "rBias"] <- relBias(Mean(est_vec), true_par_val)
      }
      
      sd_est <- SD(est_vec)
      if(!is.na(sd_est) && sd_est > 0){
        MCfile_summ[par, "RDSE"] <- (Mean(se_vec) - sd_est) / sd_est
      }
      
      MCfile_summ[par, "95%CI_LL"] <- quantile(est_vec, probs = 0.025, na.rm = T)
      MCfile_summ[par, "95%CI_UL"] <- quantile(est_vec, probs = 0.975, na.rm = T)
      
      if(!is.na(true_par_val)){
        valid_flags <- sum(!is.na(flag_vec))
        if(valid_flags > 0){
          MCfile_summ[par, "Coverage"] <- sum(flag_vec, na.rm = T) / valid_flags
        }
      }
      
      valid_pflags <- sum(!is.na(pflag_vec))
      if(valid_pflags > 0){
        detection_rate <- sum(pflag_vec, na.rm = T) / valid_pflags
        
        if(!is.na(true_par_val)){
          if(abs(true_par_val) < 1e-6){
            MCfile_summ[par, "TypeIError"] <- detection_rate
          } else {
            MCfile_summ[par, "Power"] <- detection_rate
          }
        }
      }
      
      MCfile_summ[par, "MissingPer"] <- sum(is.na(est_vec)) / N_repl
    }
    
    MCfile_summ_df <- as.data.frame(MCfile_summ)
    
    file_summ_name <- paste0(C_Model, "_TrueModel_MCfileSumm_nT", nT, "_nP", nP, ".csv")
    write.csv(MCfile_summ_df, file_summ_name)
    cat(paste0("Summary saved to: ", file_summ_name, "\n"))
    
    # Print summary
    cat("\n--- Key Parameter Summary ---\n")
    key_params <- c("Level2Mean[1]", "Level2Mean[3]", "Level2Mean[7]", "Level2Mean[9]",
                    "sigma_innovation[1,1]", "sigma_innovation[1,2]")
    for(kp in key_params){
      if(kp %in% rownames(MCfile_summ_df)){
        cat(paste0(Param_Labels[kp], ": ",
                   "True=", round(as.numeric(MCfile_summ_df[kp, "TruePar"]), 4), ", ",
                   "Est=", round(as.numeric(MCfile_summ_df[kp, "MeanPar_hat"]), 4), ", ",
                   "Coverage=", round(as.numeric(MCfile_summ_df[kp, "Coverage"]), 3), "\n"))
      }
    }
  }
}

cat("\n========== Analysis Completed Successfully! ==========\n")
