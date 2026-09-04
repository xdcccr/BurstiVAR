# ------------------------------------------------------------------
# BurstiVAR public repository copy (Supplement D, two-burst check)
# Original path: paper_simulation_archive/SupplementD_TwoBurst/code/Generate_Summ_BurstIntercept_2Burst_3T.R
# Modified for this repository: the absolute base_path was replaced with ".";
#   run this script with this folder as the R working directory.
# All other content is identical to the version that produced the
# published results.
# ------------------------------------------------------------------

########################################################################
# Script: Generate_Summ_BurstIntercept_2Burst_3T.R
# Description: Aggregates results from fitting 2-Burst BurstIntercept model.
#              Computes RMSE, Bias, Coverage, Power, Type I Error, etc.
########################################################################

rm(list=ls())

required_packages <- c("dplyr", "tidyverse", "data.table")
for(pkg in required_packages){
  if(!require(pkg, character.only = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

######################## 1. Environment & Path Setting ########################
# Roar Cluster path
# base_path <- "/storage/work/xjx5093/BurstIntercept/mlGVARNoCorNoME_BurstIntercept_2Burst_3T"

# Local PC path
base_path <- "."  # run this script from this directory

result_path <- paste0(base_path, "/result")

print(paste0("Result path: ", result_path))
if(!dir.exists(result_path)) stop("Result directory not found! Please check 'base_path'.")
setwd(result_path)

######################## 2. Helper Functions ########################
Mean = function(x){ mean(x, na.rm=TRUE) }
SD = function(x){ sd(x, na.rm=TRUE) }

# Relative Bias (avoid division by zero)
relBias = function(x, truex){ 
  if(is.na(truex) || abs(truex) < 1e-6) return(NA) 
  mean((x - truex) / truex, na.rm=T) 
}

# RMSE
RMSE = function(x, truex){ 
  if(any(is.na(truex))) return(NA)
  sqrt(Mean((x - truex)^2)) 
}

# --- Robust Column Index Finder ---
find_col_index <- function(col_names, keywords) {
  for (key in keywords) {
    idx <- grep(key, col_names, ignore.case = TRUE)
    if (length(idx) > 0) return(idx[1])
  }
  return(NA)
}

######################## 3. Simulation Conditions ########################
nT_list <- c(6)         # 2 bursts × 3 time points = 6
nP_list <- c(100, 500)  # Sample sizes
N_repl <- 100           # Number of replications

# Condition name
C <- "mlGVARNoCorNoME_BurstIntercept_2Burst_3T"

######################## 4. True Parameters ########################
# From the data generation script:

# Burst 1 Intercept
Intercept_mu_Burst1 <- c(0, 1)        # Y1=0, Y2=1
Intercept_sigma_Burst1 <- c(1, 1)     # SD: 1, 1

# Burst 2 Intercept
Intercept_mu_Burst2 <- c(1, 1.5)      # Y1=1, Y2=1.5
Intercept_sigma_Burst2 <- c(1.5, 2)   # SD: 1.5, 2

# VAR Parameters
mu_AR1 <- 0.3       # AR for reading
mu_AR2 <- 0.2       # AR for math
mu_CR12 <- -0.15    # Cross effect: reading -> math
mu_CR21 <- 0        # Cross effect: math -> reading

# SD for VAR parameters
sigma_AR <- c(0.1, 0.1)
sigma_CR <- c(0.1, 0.1)

# Innovation parameters
sigma_innovation <- 1.0
error_correlation <- 0.3
true_sigma_val <- sigma_innovation^2                       # = 1.0
true_cov_val <- error_correlation * sigma_innovation^2     # = 0.3

######################## 5. Build True Values Map ########################
# 2-Burst BurstIntercept model parameters:
# Level2Mean[1-4]: Burst intercepts (2 bursts × 2 variables)
#   [1] = Burst1 y1, [2] = Burst1 y2
#   [3] = Burst2 y1, [4] = Burst2 y2
# Level2Mean[5-8]: VAR parameters
#   [5] = AR1, [6] = AR2, [7] = CR12 (y1->y2), [8] = CR21 (y2->y1)

True_Values_Map <- c(
  # Burst intercepts
  "Level2Mean[1]" = Intercept_mu_Burst1[1],   # Burst1 y1 = 0
  "Level2Mean[2]" = Intercept_mu_Burst1[2],   # Burst1 y2 = 1
  "Level2Mean[3]" = Intercept_mu_Burst2[1],   # Burst2 y1 = 1
  "Level2Mean[4]" = Intercept_mu_Burst2[2],   # Burst2 y2 = 1.5
  
  # VAR parameters
  "Level2Mean[5]" = mu_AR1,    # AR1 = 0.3
  "Level2Mean[6]" = mu_AR2,    # AR2 = 0.2
  "Level2Mean[7]" = mu_CR12,   # CR12 = -0.15
  "Level2Mean[8]" = mu_CR21,   # CR21 = 0
  
  # Level2Sigma for burst intercepts
  "Level2Sigma[1]" = Intercept_sigma_Burst1[1],   # SD Burst1 y1 = 1
  "Level2Sigma[2]" = Intercept_sigma_Burst1[2],   # SD Burst1 y2 = 1
  "Level2Sigma[3]" = Intercept_sigma_Burst2[1],   # SD Burst2 y1 = 1.5
  "Level2Sigma[4]" = Intercept_sigma_Burst2[2],   # SD Burst2 y2 = 2
  
  # Level2Sigma for VAR parameters
  "Level2Sigma[5]" = sigma_AR[1],    # SD AR1 = 0.1
  "Level2Sigma[6]" = sigma_AR[2],    # SD AR2 = 0.1
  "Level2Sigma[7]" = sigma_CR[1],    # SD CR12 = 0.1
  "Level2Sigma[8]" = sigma_CR[2],    # SD CR21 = 0.1
  
  # Innovation covariance matrix
  "sigma_innovation[1,1]" = true_sigma_val,   # = 1.0
  "sigma_innovation[2,2]" = true_sigma_val,   # = 1.0
  "sigma_innovation[1,2]" = true_cov_val,     # = 0.3
  "sigma_innovation[2,1]" = true_cov_val      # = 0.3
)

# Parameter labels for better readability
Param_Labels <- c(
  "Level2Mean[1]" = "Burst1_y1",
  "Level2Mean[2]" = "Burst1_y2",
  "Level2Mean[3]" = "Burst2_y1",
  "Level2Mean[4]" = "Burst2_y2",
  "Level2Mean[5]" = "AR1",
  "Level2Mean[6]" = "AR2",
  "Level2Mean[7]" = "CR_y1toy2",
  "Level2Mean[8]" = "CR_y2toy1",
  "Level2Sigma[1]" = "SD_Burst1_y1",
  "Level2Sigma[2]" = "SD_Burst1_y2",
  "Level2Sigma[3]" = "SD_Burst2_y1",
  "Level2Sigma[4]" = "SD_Burst2_y2",
  "Level2Sigma[5]" = "SD_AR1",
  "Level2Sigma[6]" = "SD_AR2",
  "Level2Sigma[7]" = "SD_CR_y1toy2",
  "Level2Sigma[8]" = "SD_CR_y2toy1",
  "sigma_innovation[1,1]" = "InnoVar_y1",
  "sigma_innovation[2,2]" = "InnoVar_y2",
  "sigma_innovation[1,2]" = "InnoCov",
  "sigma_innovation[2,1]" = "InnoCov_dup"
)

# Print true values
cat("\n========== True Parameter Values ==========\n")
cat("Burst Intercepts:\n")
cat(paste0("  Level2Mean[1] (Burst1_y1) = ", True_Values_Map["Level2Mean[1]"], "\n"))
cat(paste0("  Level2Mean[2] (Burst1_y2) = ", True_Values_Map["Level2Mean[2]"], "\n"))
cat(paste0("  Level2Mean[3] (Burst2_y1) = ", True_Values_Map["Level2Mean[3]"], "\n"))
cat(paste0("  Level2Mean[4] (Burst2_y2) = ", True_Values_Map["Level2Mean[4]"], "\n"))
cat("VAR Parameters:\n")
cat(paste0("  Level2Mean[5] (AR1) = ", True_Values_Map["Level2Mean[5]"], "\n"))
cat(paste0("  Level2Mean[6] (AR2) = ", True_Values_Map["Level2Mean[6]"], "\n"))
cat(paste0("  Level2Mean[7] (CR12) = ", True_Values_Map["Level2Mean[7]"], "\n"))
cat(paste0("  Level2Mean[8] (CR21) = ", True_Values_Map["Level2Mean[8]"], "\n"))
cat("Innovation Covariance:\n")
cat(paste0("  sigma_innovation[1,1] = ", True_Values_Map["sigma_innovation[1,1]"], "\n"))
cat(paste0("  sigma_innovation[1,2] = ", True_Values_Map["sigma_innovation[1,2]"], "\n"))
cat("=============================================\n\n")

######################## 6. Target Parameters ########################
target_params <- names(True_Values_Map)

######################## 7. Main Loop ########################
for(nT in nT_list){
  for(nP in nP_list){
    
    cat(paste0("\n============ Processing nT=", nT, ", nP=", nP, " ============\n"))
    
    # Initialize MCfile matrix
    n_cols <- 1 + length(target_params) * 7  # r + 7 columns per parameter
    MCfile <- matrix(NA, nrow = N_repl, ncol = n_cols)
    
    col_names <- c("r")
    for(par in target_params){
      col_names <- c(col_names, 
                     paste0("Tru_", par),
                     paste0("Est_", par),
                     paste0("se_", par),
                     paste0("LL_", par),
                     paste0("UL_", par),
                     paste0("cFlag_", par),
                     paste0("pFlag_", par))
    }
    colnames(MCfile) <- col_names
    
    valid_reps <- 0
    
    # Loop through replications
    for(r in 1:N_repl){
      MCfile[r, "r"] <- r
      
      # File path for this replication
      res_file <- paste0(C, "_resulttable_nT", nT, "_nP", nP, "_r", r, ".csv")
      
      if(!file.exists(res_file)){
        cat(paste0("Missing file: ", res_file, "\n"))
        next
      }
      
      # Read result table
      res_table <- tryCatch({
        read.csv(res_file, row.names = 1, check.names = FALSE)
      }, error = function(e){
        cat(paste0("Error reading: ", res_file, "\n"))
        return(NULL)
      })
      
      if(is.null(res_table)) next
      
      # Detect column mapping
      cnames <- colnames(res_table)
      col_map <- list(
        mean = find_col_index(cnames, c("mean", "Mean", "Est")),
        se   = find_col_index(cnames, c("PSD", "SD", "se")),
        ll   = find_col_index(cnames, c("2.50%", "PCI 2.50%", "quantileLow", "Low")),
        ul   = find_col_index(cnames, c("97.50%", "PCI 97.50%", "quantileHigh", "High"))
      )
      
      # Extract values for each parameter
      for(par in target_params){
        if(par %in% rownames(res_table)){
          tru_val <- True_Values_Map[par]
          
          # Safe extraction
          get_val <- function(df, row, col_idx){
            if(is.na(col_idx)) return(NA)
            val <- df[row, col_idx]
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
            mean_idx <- grep("mean|Mean|Est", cnames_curr)[1]
            se_idx   <- grep("SD|PSD|se", cnames_curr)[1]
            ll_idx   <- grep("2.5|Low", cnames_curr)[1]
            ul_idx   <- grep("97.5|High", cnames_curr)[1]
            
            if(!is.na(mean_idx)) est_val <- as.numeric(res_table[par, mean_idx])
            if(!is.na(se_idx))   se_val  <- as.numeric(res_table[par, se_idx])
            if(!is.na(ll_idx))   ll_val  <- as.numeric(res_table[par, ll_idx])
            if(!is.na(ul_idx))   ul_val  <- as.numeric(res_table[par, ul_idx])
          }
          
          # Ensure scalar
          if(length(est_val) != 1) est_val <- NA
          if(length(se_val) != 1)  se_val <- NA
          if(length(ll_val) != 1)  ll_val <- NA
          if(length(ul_val) != 1)  ul_val <- NA
          
          # Coverage flag
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
              # For variance parameters: significant if CI > 0
              power_flag <- ifelse(ll_val > 0, 1, 0)
            } else {
              # For other parameters: significant if CI excludes 0
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
    
    cat(paste0("Found ", valid_reps, " valid replication files.\n"))
    
    if(valid_reps == 0){
      cat("No valid files found. Skipping this condition.\n")
      next
    }
    
    # Save Raw MCfile
    MCfile_name <- paste0(C, "_MCfile_nT", nT, "_nP", nP, ".csv")
    write.csv(MCfile, MCfile_name, row.names = FALSE)
    cat(paste0("MCfile saved to: ", MCfile_name, "\n"))
    
    ######################## 8. Compute Summary Statistics ########################
    
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
      
      # Skip if all missing
      if(all(is.na(est_vec))) next
      
      # Get true parameter value
      true_par_val <- True_Values_Map[par]
      
      # Label
      MCfile_summ[par, "Label"] <- Param_Labels[par]
      
      # True value
      MCfile_summ[par, "TruePar"] <- true_par_val
      
      # Basic statistics
      MCfile_summ[par, "MeanPar_hat"] <- Mean(est_vec)
      MCfile_summ[par, "SD"]          <- SD(est_vec) 
      MCfile_summ[par, "Mean_SE_hat"] <- Mean(se_vec)
      
      # RMSE and rBias
      if(!is.na(true_par_val)){
        MCfile_summ[par, "RMSE"]  <- RMSE(est_vec, rep(true_par_val, length(est_vec)))
        MCfile_summ[par, "rBias"] <- relBias(Mean(est_vec), true_par_val)
      }
      
      # RDSE (Relative Deviation of SE)
      sd_est <- SD(est_vec)
      if(!is.na(sd_est) && sd_est > 0){
        MCfile_summ[par, "RDSE"] <- (Mean(se_vec) - sd_est) / sd_est
      }
      
      # CI of estimates
      MCfile_summ[par, "95%CI_LL"] <- quantile(est_vec, probs = 0.025, na.rm = T)
      MCfile_summ[par, "95%CI_UL"] <- quantile(est_vec, probs = 0.975, na.rm = T)
      
      # Coverage
      if(!is.na(true_par_val)){
        valid_flags <- sum(!is.na(flag_vec))
        if(valid_flags > 0){
          MCfile_summ[par, "Coverage"] <- sum(flag_vec, na.rm = T) / valid_flags
        }
      }
      
      # Power / Type I Error
      valid_pflags <- sum(!is.na(pflag_vec))
      if(valid_pflags > 0){
        detection_rate <- sum(pflag_vec, na.rm = T) / valid_pflags
        
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
    
    # Convert to data frame
    MCfile_summ_df <- as.data.frame(MCfile_summ)
    
    # Save summary
    file_summ_name <- paste0(C, "_MCfileSumm_nT", nT, "_nP", nP, ".csv")
    write.csv(MCfile_summ_df, file_summ_name)
    cat(paste0("Summary saved to: ", file_summ_name, "\n"))
    
    # Print summary for key parameters
    cat("\n--- Key Parameter Summary ---\n")
    key_params <- c("Level2Mean[1]", "Level2Mean[2]", "Level2Mean[3]", "Level2Mean[4]",
                    "Level2Mean[5]", "Level2Mean[6]", "Level2Mean[7]", "Level2Mean[8]",
                    "sigma_innovation[1,1]", "sigma_innovation[1,2]")
    
    for(kp in key_params){
      if(kp %in% rownames(MCfile_summ_df)){
        label <- Param_Labels[kp]
        true_val <- as.numeric(MCfile_summ_df[kp, "TruePar"])
        est_val <- as.numeric(MCfile_summ_df[kp, "MeanPar_hat"])
        coverage <- as.numeric(MCfile_summ_df[kp, "Coverage"])
        power <- as.numeric(MCfile_summ_df[kp, "Power"])
        type1 <- as.numeric(MCfile_summ_df[kp, "TypeIError"])
        
        cat(paste0(label, " (", kp, "): ",
                   "True=", round(true_val, 4), ", ",
                   "Est=", round(est_val, 4), ", ",
                   "Coverage=", ifelse(is.na(coverage), "NA", round(coverage, 3))))
        
        if(!is.na(power)){
          cat(paste0(", Power=", round(power, 3)))
        }
        if(!is.na(type1)){
          cat(paste0(", TypeIErr=", round(type1, 3)))
        }
        cat("\n")
      }
    }
    
    cat("\n")
  }
}

cat("\n========== Analysis Completed Successfully! ==========\n")
