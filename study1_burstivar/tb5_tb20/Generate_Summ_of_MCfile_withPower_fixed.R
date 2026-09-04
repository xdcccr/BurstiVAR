# ----------------------------------------------------------------------------
# PUBLIC-REPO COPY. Original: IP_1b/MLGVAR2601/mlGVARNoCorNoME_BurstIntercept_CREqual0/Generate_Summ_of_MCfile_withPower_fixed.R
# Modified for this repository: the absolute base_path line was replaced with getwd()
# so the script runs from this directory (run it with the working directory
# set to study1_burstivar/tb5_tb20). All other content is identical to the
# version that produced the published results.
# ----------------------------------------------------------------------------
########################################################################
# Script: Generate_Summ_Robust.R
# Description: Robustly aggregates individual replication result tables.
#              Auto-detects column names to prevent "replacement length zero" errors.
#              Added: Power calculation for detecting significant effects
#              Added: Type I Error calculation for parameters with true value = 0
########################################################################

rm(list=ls())

required_packages <- c("dplyr", "tidyverse", "data.table")
# Check and load packages
for(pkg in required_packages){
  if(!require(pkg, character.only = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

######################## 1. Environment & Path Setting ########################
base_path <- getwd()  # REPO EDIT: run from this directory (was D:/XXYDATAanalysis/IP_1b/MLGVAR2601/mlGVARNoCorNoME_BurstIntercept_CREqual0)
result_path <- paste0(base_path, "/result")

print(paste0("Result path: ", result_path))
if(!dir.exists(result_path)) stop("Result directory not found! Please check 'base_path'.")
setwd(result_path)

######################## 2. Helper Functions ########################
Mean = function(x){ mean(x, na.rm=TRUE) }
SD = function(x){ sd(x, na.rm=TRUE) }
# Relative Bias (avoid division by zero)
relBias = function(x, truex){ 
  if(abs(truex) < 1e-6) return(NA) 
  mean((x - truex) / truex, na.rm=T) 
}
# RMSE
RMSE = function(x, truex){ sqrt(Mean((x - truex)^2)) }

# --- Robust Value Extractor ---
# Finds the correct column index based on keywords
find_col_index <- function(col_names, keywords) {
  for (key in keywords) {
    # Case-insensitive partial match
    idx <- grep(key, col_names, ignore.case = TRUE)
    if (length(idx) > 0) return(idx[1]) # Return first match
  }
  return(NA) # Not found
}

######################## 3. Simulation Conditions ########################
nT_list <- c(15,60)       
nP_list <- c(100,500)      
N_repl <- 100          
Condition_Name <- "mlGVARNoCorNoME_BurstIntercept_CREqual0"

######################## 4. True Parameters Setting (Matches DataGen) ########################
mu_Burst1_y1 = 0; mu_Burst1_y2 = 1
mu_Burst2_y1 = 1; mu_Burst2_y2 = 1.5
mu_Burst3_y1 = 2; mu_Burst3_y2 = 2
mu_AR1 = 0.3; mu_AR2 = 0.2
mu_CR12 = -0.15; mu_CR21 = 0#-0.1
true_sigma_val = 1.0^2
true_cov_val   = 0.3 * 1.0^2


sigma_Burst1_y1 = 1
sigma_Burst1_y2 = 1
sigma_Burst2_y1 = 1.5
sigma_Burst2_y2 = 2
sigma_Burst3_y1 = 2
sigma_Burst3_y2 = 3
sigma_AR1 = 0.1      
sigma_AR2 = 0.1      
sigma_CR12 = 0.1       
sigma_CR21 = 0.1  


True_Values_Map <- c(
  "Level2Mean[1]"  = mu_Burst1_y1, 
  "Level2Mean[2]"  = mu_Burst1_y2, 
  "Level2Mean[3]"  = mu_Burst2_y1, 
  "Level2Mean[4]"  = mu_Burst2_y2, 
  "Level2Mean[5]"  = mu_Burst3_y1, 
  "Level2Mean[6]"  = mu_Burst3_y2, 
  "Level2Mean[7]"  = mu_AR1,       
  "Level2Mean[8]"  = mu_AR2,       
  "Level2Mean[9]"  = mu_CR12,      
  "Level2Mean[10]" = mu_CR21,
  
  "Level2Sigma[1]"  = sigma_Burst1_y1, 
  "Level2Sigma[2]"  = sigma_Burst1_y2, 
  "Level2Sigma[3]"  = sigma_Burst2_y1, 
  "Level2Sigma[4]"  = sigma_Burst2_y2, 
  "Level2Sigma[5]"  = sigma_Burst3_y1, 
  "Level2Sigma[6]"  = sigma_Burst3_y2, 
  "Level2Sigma[7]"  = sigma_AR1,       
  "Level2Sigma[8]"  = sigma_AR2,       
  "Level2Sigma[9]"  = sigma_CR12,      
  "Level2Sigma[10]" = sigma_CR21,  
  
  "sigma_innovation[1,1]" = true_sigma_val,
  "sigma_innovation[2,2]" = true_sigma_val,
  "sigma_innovation[1,2]" = true_cov_val,
  "sigma_innovation[2,1]" = true_cov_val
)

######################## 5. Main Loop ########################

for(nT in nT_list){
  for(nP in nP_list){
    
    cat(paste0("\nProcessing: nT=", nT, ", nP=", nP, "\n"))
    
    target_params <- names(True_Values_Map)
    # Added "pFlag" for power/Type I error calculation
    summ_cols_base <- c("Tru", "Est", "se", "LL", "UL", "cFlag", "pFlag")
    all_colnames <- c("r", as.vector(outer(summ_cols_base, target_params, paste, sep="_")))
    
    MCfile <- matrix(NA, nrow = N_repl, ncol = length(all_colnames))
    colnames(MCfile) <- all_colnames
    
    # ---------------------------------------------------------
    # First Pass: Detect Column Names from the first valid file
    # ---------------------------------------------------------
    col_map <- list()
    # Use precise regex: _r1.csv$ to match exactly "_r1.csv" at the end
    first_file <- list.files(result_path, pattern = paste0(Condition_Name, ".*_r1\\.csv$"), full.names = T)
    # Take only the first match if multiple files found
    if(length(first_file) > 1) first_file <- first_file[1]
    if(length(first_file) == 0) first_file <- list.files(result_path, pattern = "\\.csv$", full.names = T)[1]
    
    if(length(first_file) == 1 && !is.na(first_file) && file.exists(first_file)){
      # Try reading to find headers
      temp_df <- read.csv(first_file)
      # Check if first column is rownames (often "X", "Unnamed", or empty)
      cnames <- colnames(temp_df)
      
      # Determine mapping based on detected column names
      col_map$mean <- find_col_index(cnames, c("Mean", "mean", "Estimate", "Est"))
      col_map$se   <- find_col_index(cnames, c("SD", "PSD", "se", "StdDev"))
      col_map$ll   <- find_col_index(cnames, c("2.5", "Low", "LL", "min"))
      col_map$ul   <- find_col_index(cnames, c("97.5", "High", "UL", "max"))
      
      cat("Detected Column Mapping (Index):\n")
      print(unlist(col_map))
      cat("Corresponding Names:\n")
      print(cnames[unlist(col_map)])
    } else {
      stop("No result files found to initialize column mapping.")
    }
    
    # ---------------------------------------------------------
    # Processing Loop
    # ---------------------------------------------------------
    valid_reps <- 0
    
    for(r in 1:N_repl){
      file_name <- paste0(Condition_Name, "_resulttable_nT", nT, "_nP", nP, "_r", r, ".csv")
      full_path <- file.path(result_path, file_name)
      
      if(file.exists(full_path)){
        # Read file, assume first column is Row Names
        res_table <- read.csv(full_path, row.names = 1) 
        
        # Ensure row names are clean (remove leading spaces if any)
        rownames(res_table) <- trimws(rownames(res_table))
        
        MCfile[r, "r"] <- r
        
        for(par in target_params){
          if(par %in% rownames(res_table)){
            tru_val <- True_Values_Map[par]
            
            # Use detected column indices to extract values robustly
            # Note: read.csv with row.names=1 reduces column count by 1 in the data frame object
            # We need to adjust indices if we used the raw CSV headers for detection
            # Easier approach: Use names found in detection
            
            # Safe extraction function
            get_val <- function(df, row, col_idx){
              val <- df[row, col_idx - 1] # -1 because row.names took col 1
              if(length(val) == 0 || is.null(val)) return(NA)
              return(as.numeric(val))
            }
            
            est_val <- get_val(res_table, par, col_map$mean)
            se_val  <- get_val(res_table, par, col_map$se)
            ll_val  <- get_val(res_table, par, col_map$ll)
            ul_val  <- get_val(res_table, par, col_map$ul)
            
            # Fallback: if indices failed, try name matching directly on the loaded DF
            if(is.na(est_val)){
              cnames_curr <- colnames(res_table)
              est_val <- res_table[par, grep("mean|Mean|Est", cnames_curr)[1]]
              se_val  <- res_table[par, grep("SD|PSD|se", cnames_curr)[1]]
              ll_val  <- res_table[par, grep("2.5|Low", cnames_curr)[1]]
              ul_val  <- res_table[par, grep("97.5|High", cnames_curr)[1]]
            }
            
            # Safety check: ensure scalar (length 1)
            if(length(est_val) != 1) est_val <- NA
            if(length(se_val) != 1) se_val <- NA
            if(length(ll_val) != 1) ll_val <- NA
            if(length(ul_val) != 1) ul_val <- NA
            
            # Coverage
            cov_flag <- NA
            if(!is.na(tru_val) & !is.na(ll_val) & !is.na(ul_val)){
              cov_flag <- ifelse(tru_val >= ll_val & tru_val <= ul_val, 1, 0)
            }
            
            # ----------------------------------------------------------
            # Power / Type I Error Flag Calculation (NEW)
            # ----------------------------------------------------------
            # This flag indicates whether the 95% CI excluded zero
            # Interpretation depends on true value:
            #   - If true value != 0: This is POWER (correctly detecting a real effect)
            #   - If true value == 0: This is TYPE I ERROR (falsely detecting a non-existent effect)
            #
            # For Level2Sigma (variance parameters): LL > 0 indicates significant individual differences
            # For Level2Mean (mean parameters): LL > 0 OR UL < 0 indicates significant effect
            # For sigma_innovation covariance: LL > 0 OR UL < 0 indicates significant covariance
            
            power_flag <- NA
            if(!is.na(ll_val) & !is.na(ul_val)){
              # Check if parameter is a variance parameter (Level2Sigma or diagonal of sigma_innovation)
              is_variance_param <- grepl("Level2Sigma", par) | 
                                   par == "sigma_innovation[1,1]" | 
                                   par == "sigma_innovation[2,2]"
              
              if(is_variance_param){
                # For variance parameters: significant if LL > 0
                power_flag <- ifelse(ll_val > 0, 1, 0)
              } else {
                # For mean/covariance parameters: significant if CI excludes 0
                power_flag <- ifelse(ll_val > 0 | ul_val < 0, 1, 0)
              }
            }
            # ----------------------------------------------------------
            
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
    
    # Save Raw MCfile
    write.csv(MCfile, paste0(Condition_Name, "_MCfile_nT", nT, "_nP", nP, ".csv"), row.names = FALSE)
    
    ######################## 6. Compute Summary Statistics ########################
    
    # Added "Power" and "TypeIError" columns
    # Power: for parameters with true value != 0
    # TypeIError: for parameters with true value == 0
    summ_stat_cols <- c("TruePar", "MeanPar_hat", "RMSE", "rBias", "SD", 
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
      pflag_vec <- as.numeric(MC_df[[paste0("pFlag_", par)]])  # Power/TypeI flag vector
      
      # Skip if all missing
      if(all(is.na(est_vec))) next
      
      # Get true parameter value
      true_par_val <- mean(tru_vec, na.rm=T)
      
      # Calculate Stats
      MCfile_summ[par, "TruePar"]     <- true_par_val
      MCfile_summ[par, "MeanPar_hat"] <- Mean(est_vec)
      MCfile_summ[par, "SD"]          <- SD(est_vec) 
      MCfile_summ[par, "Mean_SE_hat"] <- Mean(se_vec)
      MCfile_summ[par, "RMSE"]        <- RMSE(est_vec, tru_vec)
      MCfile_summ[par, "rBias"]       <- relBias(Mean(est_vec), mean(tru_vec, na.rm=T))
      
      # RDSE
      MCfile_summ[par, "RDSE"] <- (Mean(se_vec) - SD(est_vec)) / SD(est_vec)
      
      # CI
      MCfile_summ[par, "95%CI_LL"] <- quantile(est_vec, probs = 0.025, na.rm = T)
      MCfile_summ[par, "95%CI_UL"] <- quantile(est_vec, probs = 0.975, na.rm = T)
      MCfile_summ[par, "Coverage"] <- sum(flag_vec, na.rm = T) / sum(!is.na(flag_vec)) # Coverage rate among valid
      
      # ----------------------------------------------------------
      # Power vs Type I Error Calculation (NEW)
      # ----------------------------------------------------------
      # Calculate the detection rate (proportion of CIs excluding zero)
      detection_rate <- sum(pflag_vec, na.rm = T) / sum(!is.na(pflag_vec))
      
      # Assign to appropriate column based on true value
      # Use small threshold (1e-6) to handle floating point comparison
      if(abs(true_par_val) < 1e-6){
        # True value is essentially 0 -> This is Type I Error
        MCfile_summ[par, "Power"] <- NA
        MCfile_summ[par, "TypeIError"] <- detection_rate
      } else {
        # True value is not 0 -> This is Power
        MCfile_summ[par, "Power"] <- detection_rate
        MCfile_summ[par, "TypeIError"] <- NA
      }
      # ----------------------------------------------------------
      
      MCfile_summ[par, "MissingPer"] <- sum(is.na(est_vec)) / N_repl
    }
    
    file_summ_name <- paste0(Condition_Name, "_MCfileSumm_nT", nT, "_nP", nP, ".csv")
    write.csv(MCfile_summ, file_summ_name)
    cat(paste0("Summary saved to: ", file_summ_name, "\n"))
  }
}

print("Analysis Completed Successfully!")
