# --------------------------------------------------------------------
# BurstiVAR public repository -- Study 2 (GoBurstiVAR).
# Original file: paper_simulation_archive/Study2_Gompertz/code/mlGVARNoCorNoME_GompertzBurst_DataGeneration.R
# Modification for this repository: base_path set to "." (was an
# absolute local path); run this script from this directory. Simulated
# data are written to ./data. The nT/nP loop headers are preserved as
# last edited; set them to the target cell before running (see README).
# All other content is identical to the version that produced the
# published results.
# --------------------------------------------------------------------

######################## Set Environment ########################
# Note: mlGVARNoCorNoME_GompertzBurst
# Trend: Each burst has a constant value, but values across bursts follow Gompertz curve
#        trend_b = theta1 * exp(-theta2 * exp(-b * theta3)), where b = burst number (1,2,3,4,5)
# VAR: Standard VAR(1) process for deviations

rm(list=ls())
library(mvtnorm)
library(Matrix)
library(dplyr)

# Set paths - modify for your environment
# base_path <- "/storage/work/xjx5093/BurstIntercept/mlGVARNoCorNoME_GompertzBurst"  # Roar Cluster
base_path <- "."   # Local PC (repo: run this script from this directory)
work_path <- paste0(base_path, "/data")
print(paste0("work_path: ", work_path))

# Create directory if it doesn't exist
if(!dir.exists(work_path)){
  dir.create(work_path, recursive = TRUE)
}
setwd(work_path)

######################## Control Parameters ########################
C = "mlGVARNoCorNoME_GompertzBurst"

for(nT in c(100)){  # Total time points (must be divisible by nBurst)
  for(nP in c(500)){  # Number of people
    ydim = 2       # number of dimensions (y1 and y2)
    npad = 0       # initial number of time points to discard
    N_repl = 100   # number of replications
    rs = 1:N_repl
    nBurst = 5     # number of bursts
    
    ######################## Model Parameters ########################
    # 1. Gompertz Parameters for Burst-level Trend
    # The trend at burst b is: theta1 * exp(-theta2 * exp(-b * theta3))
    # where b = 1, 2, 3, 4, 5 (burst number)
    
    # Variable 1 (y1) Gompertz parameters
    theta10_1 = 10    # asymptote (upper limit of growth)
    theta20_1 = 3     # displacement (controls where growth starts)
    theta30_1 = 1.0   # growth rate (how fast it approaches asymptote)
    
    # Variable 2 (y2) Gompertz parameters
    theta10_2 = 8     # asymptote
    theta20_2 = 2.5   # displacement
    theta30_2 = 0.8   # growth rate
    
    # Standard deviations for Gompertz parameters (between-person variability)
    sigma1_1 = 2      # SD for y1 asymptote
    sigma2_1 = 0.5    # SD for y1 displacement
    sigma3_1 = 0.2    # SD for y1 growth rate
    
    sigma1_2 = 1.5    # SD for y2 asymptote
    sigma2_2 = 0.4    # SD for y2 displacement
    sigma3_2 = 0.15   # SD for y2 growth rate
    
    # 2. VAR Parameters
    # Mean parameters
    AR_mean1 = 0.3    # AR for y1
    AR_mean2 = 0.2    # AR for y2
    CR12_mean = -0.15 # Cross effect: y1 -> y2
    CR21_mean = -0.1  # Cross effect: y2 -> y1
    
    # Standard deviations for VAR parameters
    sigma_AR = 0.1    # SD for AR coefficients
    sigma_CR = 0.1    # SD for cross-regression coefficients
    
    # Innovation parameters
    sigma_innovation = 1.0  # innovation SD
    error_correlation = 0.3 # innovation correlation
    
    # Complete parameter vector (10 parameters total)
    # Parameters 1-3: y1 Gompertz (theta1, theta2, theta3)
    # Parameters 4-6: y2 Gompertz (theta1, theta2, theta3)
    # Parameters 7-8: AR coefficients
    # Parameters 9-10: CR coefficients
    param_means = c(
      theta10_1, theta20_1, theta30_1,  # y1 Gompertz means (1-3)
      theta10_2, theta20_2, theta30_2,  # y2 Gompertz means (4-6)
      AR_mean1, AR_mean2,               # AR means (7-8)
      CR12_mean, CR21_mean              # CR means (9-10)
    )
    
    # Standard deviation vector
    param_sds = c(
      sigma1_1, sigma2_1, sigma3_1,     # y1 Gompertz SDs (1-3)
      sigma1_2, sigma2_2, sigma3_2,     # y2 Gompertz SDs (4-6)
      sigma_AR, sigma_AR,               # AR SDs (7-8)
      sigma_CR, sigma_CR                # CR SDs (9-10)
    )
    
    ######################## Correlation Structure ########################
    # Build the base correlation matrix (10x10) - assuming independence
    R_full = diag(10)
    
    # Uncomment below to add correlations if needed:
    # # Within-Gompertz correlations for y1 (1:3, 1:3)
    # r12_g = 0.3  # correlation between asymptote and displacement
    # r13_g = 0.2  # correlation between asymptote and growth rate
    # r23_g = 0.1  # correlation between displacement and growth rate
    # R_full[1,2] = R_full[2,1] = r12_g
    # R_full[1,3] = R_full[3,1] = r13_g
    # R_full[2,3] = R_full[3,2] = r23_g
    
    # Convert correlation matrix to covariance matrix
    Sigma_full = diag(param_sds) %*% R_full %*% diag(param_sds)
    
    # Innovation covariance matrix (fixed across subjects)
    Sigma_innovation = matrix(
      c(sigma_innovation^2, 
        error_correlation*sigma_innovation^2,
        error_correlation*sigma_innovation^2, 
        sigma_innovation^2), 
      2, 2
    )
    
    ######################## Gompertz Function ########################
    # Gompertz curve: f(b) = theta1 * exp(-theta2 * exp(-b * theta3))
    # b is the burst number (1, 2, 3, 4, 5)
    gompertz <- function(b, theta1, theta2, theta3) {
      theta1 * exp(-theta2 * exp(-b * theta3))
    }
    
    ########################################################################
    ######################## MC Experiments begins ########################
    ########################################################################
    Total_t_begin <- proc.time()
    
    # Calculate time points per burst
    T_per_burst = nT / nBurst
    if(T_per_burst != floor(T_per_burst)){
      stop("nT must be divisible by nBurst!")
    }
    
    # Define burst boundaries
    burst_starts = seq(1, nT, by = T_per_burst)
    burst_ends = seq(T_per_burst, nT, by = T_per_burst)
    
    cat("Burst structure:\n")
    for(b in 1:nBurst){
      cat(paste0("  Burst ", b, ": time points ", burst_starts[b], " to ", burst_ends[b], "\n"))
    }
    
    for(k in 1:N_repl){
      r = rs[k]
      set.seed(r)
      
      dat_name <- paste0("Data_",C,"_nT",nT,"_nP",nP,"_r",r)
      dat_filename = paste0(dat_name,".csv")
      print(paste0("dat_filename: ", dat_filename))
      
      ######################## Data Generation ########################
      # Storage arrays
      Y = array(NA, c(nP, nT, 2))          # Final observed data [nP persons, nT times, 2 variables]
      Y_trend = array(NA, c(nP, nT, 2))    # Trend component
      Y_dev = array(NA, c(nP, nT, 2))      # Deviation from trend (VAR component)
      person_params = matrix(NA, nP, 10)   # Person-specific parameters
      burst_trends = array(NA, c(nP, nBurst, 2))  # Store burst-level trends for each person
      
      # Generate person-specific parameters
      for(i in 1:nP){
        person_params[i,] = rmvnorm(1, param_means, Sigma_full)
      }
      
      # Extract Gompertz parameters for easier reference
      theta1_y1 = person_params[,1]  # y1 asymptote
      theta2_y1 = person_params[,2]  # y1 displacement
      theta3_y1 = person_params[,3]  # y1 growth rate
      
      theta1_y2 = person_params[,4]  # y2 asymptote
      theta2_y2 = person_params[,5]  # y2 displacement
      theta3_y2 = person_params[,6]  # y2 growth rate
      
      # VAR parameters
      AR1 = person_params[,7]   # y1 AR
      AR2 = person_params[,8]   # y2 AR
      CR12 = person_params[,9]  # y1 -> y2
      CR21 = person_params[,10] # y2 -> y1
      
      # Time vector
      time = seq(0, 9.9, length.out = nT)
      
      # Generate data for each person
      for(i in 1:nP){
        
        # Calculate burst-level trends using Gompertz curve
        for(b in 1:nBurst){
          burst_trends[i, b, 1] = gompertz(b, theta1_y1[i], theta2_y1[i], theta3_y1[i])
          burst_trends[i, b, 2] = gompertz(b, theta1_y2[i], theta2_y2[i], theta3_y2[i])
        }
        
        # Loop through bursts
        for(b in 1:nBurst){
          t_start = burst_starts[b]
          t_end = burst_ends[b]
          
          # First time point of this burst
          Y_trend[i, t_start, 1:2] = burst_trends[i, b, ]
          Y_dev[i, t_start, ] = rmvnorm(1, mean = c(0,0), sigma = Sigma_innovation)
          Y[i, t_start, ] = Y_trend[i, t_start, ] + Y_dev[i, t_start, ]
          
          # Subsequent time points within this burst
          if(t_start < t_end){
            for(t in (t_start + 1):t_end){
              # Trend is constant within burst
              Y_trend[i, t, 1:2] = burst_trends[i, b, ]
              
              # VAR process for deviations
              mean_dev = c(
                AR1[i] * Y_dev[i, t-1, 1] + CR21[i] * Y_dev[i, t-1, 2],  # y1
                AR2[i] * Y_dev[i, t-1, 2] + CR12[i] * Y_dev[i, t-1, 1]   # y2
              )
              
              # Add innovations
              Y_dev[i, t, ] = rmvnorm(1, mean = mean_dev, sigma = Sigma_innovation)
              
              # Combine trend and deviation
              Y[i, t, ] = Y_trend[i, t, ] + Y_dev[i, t, ]
            }
          }
        } # End burst loop
      } # End person loop
      
      ######################## Save Data ########################
      # Reshape to long format
      dat_long <- matrix(NA, nrow = nP*nT, ncol = 5)  # [id, time, burst, y1, y2]
      colnames(dat_long) <- c("id", "time", "burst", "y1", "y2")
      
      # Fill data
      for(i in 1:nP){
        for(t in 1:nT){
          row_idx = (i-1)*nT + t
          dat_long[row_idx, "id"] = i
          dat_long[row_idx, "time"] = time[t]
          dat_long[row_idx, "burst"] = ceiling(t / T_per_burst)  # Burst number
          dat_long[row_idx, "y1"] = Y[i,t,1]
          dat_long[row_idx, "y2"] = Y[i,t,2]
        }
      }
      
      dat_long = as.data.frame(dat_long)
      
      # Save data as CSV
      write.csv(dat_long, paste0(work_path, "/", dat_name, ".csv"), 
                row.names = FALSE)
      
      print(paste0("Completed replication ", k, " for N=", nP, ", T=", nT))
      
    }  # End replication loop
    
    print(paste0("Completed condition N=", nP, ", T=", nT))
    
  }  # End nP loop
}  # End nT loop

########################################################################
######################## MC Experiments ends ########################
########################################################################
Total_t_end <- proc.time()
Total_t_elapsed <- Total_t_end - Total_t_begin
show(Total_t_elapsed)

########################################################################
######################## Display Example Gompertz Curves ##############
########################################################################
# # Show what the Gompertz curves look like with the specified parameters
# cat("\n\n========== Example Gompertz Curves ==========\n")
# cat("Burst-level trend values (population mean parameters):\n\n")
# 
# burst_nums = 1:nBurst
# 
# # Calculate population mean trends
# pop_trend_y1 = sapply(burst_nums, function(b) gompertz(b, theta10_1, theta20_1, theta30_1))
# pop_trend_y2 = sapply(burst_nums, function(b) gompertz(b, theta10_2, theta20_2, theta30_2))
# 
# cat("Variable 1 (y1):\n")
# cat(paste0("  theta1=", theta10_1, " (asymptote), theta2=", theta20_1, 
#            " (displacement), theta3=", theta30_1, " (growth rate)\n"))
# for(b in 1:nBurst){
#   cat(paste0("  Burst ", b, ": trend = ", round(pop_trend_y1[b], 3), "\n"))
# }
# 
# cat("\nVariable 2 (y2):\n")
# cat(paste0("  theta1=", theta10_2, " (asymptote), theta2=", theta20_2, 
#            " (displacement), theta3=", theta30_2, " (growth rate)\n"))
# for(b in 1:nBurst){
#   cat(paste0("  Burst ", b, ": trend = ", round(pop_trend_y2[b], 3), "\n"))
# }
# 
# cat("\n==============================================\n")