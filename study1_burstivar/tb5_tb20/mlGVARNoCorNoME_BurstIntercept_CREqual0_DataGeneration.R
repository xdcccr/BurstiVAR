# ----------------------------------------------------------------------------
# PUBLIC-REPO COPY. Original: IP_1b/MLGVAR2601/mlGVARNoCorNoME_BurstIntercept_CREqual0/mlGVARNoCorNoME_BurstIntercept_CREqual0_DataGeneration.R
# Modified for this repository: the absolute work_path line was replaced with file.path(getwd(), "data")
# so the script runs from this directory (run it with the working directory
# set to study1_burstivar/tb5_tb20). All other content is identical to the
# version that produced the published results.
# ----------------------------------------------------------------------------
######################## Set Environment ########################
rm(list=ls())
library(mvtnorm)
library(Matrix)
library(dplyr)

work_path <- file.path(getwd(), "data")  # REPO EDIT: run from this directory; create the data subfolder first (was D:/XXYDATAanalysis/IP_1b/MLGVAR2601/mlGVARNoCorNoME_BurstIntercept_CREqual0/data)
print(paste0("work_path: ", work_path))
setwd(work_path)

######################## Control Parameters ########################
C = "mlGVARNoCorNoME_BurstIntercept_CREqual0"

for(nT in c(15,60)){  # Using only nT=100 for initial validation
  for(nP in c(100,500)){  # Using only nP=500 for initial validation
    ydim = 2    # number of dimensions (reading and math)
    npad = 0    # initial number of time points to discard
    N_repl = 100  # number of replications
    rs = 1:N_repl
    
    ######################## Model Parameters ########################
    # 1. Gompertz Parameters
    # Reading parameters (slower growth)
    # theta10_1 = 0#40   # asymptote
    # theta20_1 = 0#2    # displacement for growth onset
    # theta30_1 = 0#0.5  # growth rate
    # 
    # # Math parameters (even slower growth)
    # theta10_2 = 0#30   # lower asymptote
    # theta20_2 = 0#2.5  # later displacement
    # theta30_2 = 0#0.6  # slower growth rate
    # 
    # # Standard deviations for Gompertz parameters
    # sigma1_1 = 0#6     # SD for asymptote
    # sigma2_1 = 0#0.4   # SD for displacement
    # sigma3_1 = 0#0.05  # SD for growth rate
    # 
    # # Same SDs for math trends
    # sigma1_2 = sigma1_1
    # sigma2_2 = sigma2_1
    # sigma3_2 = sigma3_1
    
    Intercept_mu_Burst1 <- c(0,1)
    Intercept_mu_Burst2 <- c(1,1.5)
    Intercept_mu_Burst3 <- c(2,2)
    
    Intercept_sigma_Burst1 <- c(1,1)
    Intercept_sigma_Burst2 <- c(1.5,2)
    Intercept_sigma_Burst3 <- c(2,3)
    # Intercept_sigma_Burst1 <- c(1,1.5)
    # Intercept_sigma_Burst2 <- c(1,1.5)
    # Intercept_sigma_Burst3 <- c(1,1.5)
    # 
    # 2. VAR Parameters
    # Mean parameters
    AR_mean1 = 0.3    # AR for reading
    AR_mean2 = 0.2    # AR for math
    CR12_mean = -0.15 # Cross effect: reading -> math
    CR21_mean = 0#-0.1  # Cross effect: math -> reading
  
    
    # Standard deviations for VAR parameters
    sigma_AR = c(0.1, 0.1)    # SD for AR coefficients
    sigma_CR = c(0.1, 0.1)    # SD for cross-regression coefficients
    
    # Innovation parameters
    sigma_innovation = 1.0  # innovation SD
    error_correlation = 0.3 # innovation correlation
    
    # Measurement Error parameters
    # sigma_ME = 1.0  # ME SD
    # ME_correlation = 0 # ME correlation
    # 
    # Complete parameter vector (10 parameters total)
    param_means = c(
      # Gompertz parameters
      # theta10_1, theta20_1, theta30_1,  # reading trend means (1-3)
      # theta10_2, theta20_2, theta30_2,  # math trend means (4-6)
      Intercept_mu_Burst1,
      Intercept_mu_Burst2,
      Intercept_mu_Burst3,
      # VAR parameters
      AR_mean1, AR_mean2,              # AR means (7-8)
      CR12_mean, CR21_mean             # CR means (9-10)
    )
    
    # Standard deviation vector
    param_sds = c(
      # Gompertz SDs
      # sigma1_1, sigma2_1, sigma3_1,     # reading trend SDs (1-3)
      # sigma1_2, sigma2_2, sigma3_2,     # math trend SDs (4-6)
      Intercept_sigma_Burst1,
      Intercept_sigma_Burst2,
      Intercept_sigma_Burst3,
      # VAR SDs
      sigma_AR[1], sigma_AR[2],               # AR SDs (7-8)
      sigma_CR[1], sigma_CR[2]                # CR SDs (9-10)
    )
    ######################## Correlation Structure ########################
    # Build the base correlation matrix (10x10)
    R_full = diag(10)
    # 
    # # 1. Within-Gompertz correlations for reading (1:3, 1:3)
    # r12_g = 0  # correlation between asymptote and displacement
    # r13_g = 0  # correlation between asymptote and growth rate
    # r23_g = 0  # correlation between displacement and growth rate
    # 
    # R_full[1,2] = R_full[2,1] = r12_g
    # R_full[1,3] = R_full[3,1] = r13_g
    # R_full[2,3] = R_full[3,2] = r23_g
    # 
    # # 2. Within-Gompertz correlations for math (4:6, 4:6)
    # R_full[4,5] = R_full[5,4] = r12_g
    # R_full[4,6] = R_full[6,4] = r13_g
    # R_full[5,6] = R_full[6,5] = r23_g
    # 
    # # 3. Between-subject Gompertz correlations (1:3, 4:6)
    # r_theta1 = 0  # correlation between reading and math asymptotes
    # r_theta2 = 0  # correlation between reading and math displacements
    # r_theta3 = 0  # correlation between reading and math growth rates
    # 
    # R_full[1,4] = R_full[4,1] = r_theta1
    # R_full[2,5] = R_full[5,2] = r_theta2
    # R_full[3,6] = R_full[6,3] = r_theta3
    # 
    # # 4. Within-VAR correlations (7:10, 7:10)
    # r_AR = 0      # correlation between AR parameters
    # r_CR = 0     # correlation between CR parameters
    # 
    # R_full[7,8] = R_full[8,7] = r_AR    # AR1-AR2
    # R_full[9,10] = R_full[10,9] = r_CR   # CR12-CR21
    # 
    # # 5. Gompertz-VAR correlations 
    # # Only within same subject, asymptote with AR
    # r_trend_AR = 0  # correlation between trend level and AR
    # 
    # # Reading: asymptote-AR correlation
    # R_full[1,7] = R_full[7,1] = r_trend_AR
    # # Math: asymptote-AR correlation
    # R_full[4,8] = R_full[8,4] = r_trend_AR
    
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
    
    # ME covariance matrix (fixed across subjects)
    # Sigma_ME = matrix(
    #   c(sigma_ME^2, 
    #     ME_correlation*sigma_ME^2,
    #     ME_correlation*sigma_ME^2, 
    #     sigma_ME^2), 
    #   2, 2
    # )
    
    # # Innovation covariance matrix for the initial condition
    # IIV_init <- diag(1, nrow=2, ncol=2)
    # Initial Condition
    # mu0_y1 = 1
    # mu0_y2 = 1
    # Sigma0 = diag(1,ydim)
    
    # Intercept_mu <- c(0,1)
    # Intercept_sigma <- c(1,1)
    # Intercept_Sigma <- matrix(c(Intercept_sigma[1],0,
    #                             0,Intercept_sigma[2]),nrow=2,byrow=T)
    # Slope_mu <- c(1,1)
    # Slope_sigma <- c(0.5,0.5)
    # Slope_Sigma <- matrix(c(Slope_sigma[1],0,
    #                         0,Slope_sigma[2]),nrow=2,byrow=T)
    
    ########################################################################
    ######################## MC Experiments begins ########################
    ########################################################################
    Total_t_begin <- proc.time()
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
      person_params = matrix(NA, nP, 10)    # Person-specific parameters
      #person_params_Intercept = matrix(NA, nP, 6) 
      #person_params_Slope = matrix(NA, nP, 2) 
      # Generate person-specific parameters
      for(i in 1:nP){
        person_params[i,] = rmvnorm(1, param_means, Sigma_full)
        #person_params_Intercept[i,] = rmvnorm(1, Intercept_mu, Intercept_Sigma) # covariance matrix
        #person_params_Slope[i,] = rmvnorm(1, Slope_mu, Slope_Sigma) # covariance matrix
      }
      
      # Extract parameters for easier reference
      # # Gompertz parameters
      # theta1_1 = person_params[,1]  # reading asymptote
      # theta2_1 = person_params[,2]  # reading displacement
      # theta3_1 = person_params[,3]  # reading growth rate
      # 
      # theta1_2 = person_params[,4]  # math asymptote
      # theta2_2 = person_params[,5]  # math displacement
      # theta3_2 = person_params[,6]  # math growth rate
      
      # Intercept parameters
      Intercept1_Burst1 = person_params[,1]
      Intercept2_Burst1 = person_params[,2]
      Intercept1_Burst2 = person_params[,3]
      Intercept2_Burst2 = person_params[,4]
      Intercept1_Burst3 = person_params[,5]
      Intercept2_Burst3 = person_params[,6]
      
      # VAR parameters
      AR1 = person_params[,7]   # reading AR
      AR2 = person_params[,8]   # math AR
      CR12 = person_params[,9]  # reading -> math
      CR21 = person_params[,10] # math -> reading
      
      time = seq(0, 9.9, length.out=nT) 
      T1 = nT*1/3
      T2 = nT*2/3
      T3 = nT*3/3
      # Generate data for each person
      # Generate data for each person
      for(i in 1:nP){
        
        #1. Generate first burst1 observations
        Y_trend[i,1,1:2] = person_params[i,1:2] #Intercept1_Burst1, Intercept2_Burst1
        Y_dev[i,1,] = rmvnorm(1, mean=c(0,0), sigma=Sigma_innovation)
        Y[i,1,] = Y_dev[i,1,]+Y_trend[i,1,]
        
        for(t in 2:T1){ # Burst 1
          # 2. Generate Gompertz trends
          #Y_trend[i,t,1] = 0#theta1_1[i] * exp(-theta2_1[i] * exp(-time[t] * theta3_1[i])) # reading
          #Y_trend[i,t,2] = 0#theta1_2[i] * exp(-theta2_2[i] * exp(-time[t] * theta3_2[i])) # math
          #Y_trend[i,t,1:2] =  person_params_Intercept[i,] + person_params_Slope[i,]*time[t]
          Y_trend[i,t,1:2] = person_params[i,1:2]
          # 3. Generate subsequent deviations and observations
          # VAR process for deviations
          mean_dev = c(
            AR1[i] * Y_dev[i,t-1,1] + CR21[i] * Y_dev[i,t-1,2],# reading
            AR2[i] * Y_dev[i,t-1,2] + CR12[i] * Y_dev[i,t-1,1] # math
          )
          
          # Add innovations
          Y_dev[i,t,] = rmvnorm(1, mean=mean_dev, sigma=Sigma_innovation)
          
          # Combine trend and deviation
          #Y[i,t,] = Y_trend[i,t,] + Y_dev[i,t,]
          #Y[i,t,] = rmvnorm(1, mean=(Y_trend[i,t,] + Y_dev[i,t,]), sigma=Sigma_ME)
          Y[i,t,] = Y_trend[i,t,] + Y_dev[i,t,]
          
        } #Closing Burst1 loop
        
        #2. Generate first burst2 observations
        Y_trend[i,T1+1,1:2] = person_params[i,3:4] #Intercept1_Burst1, Intercept2_Burst1
        Y_dev[i,T1+1,] = rmvnorm(1, mean=c(0,0), sigma=Sigma_innovation)
        Y[i,T1+1,] = Y_dev[i,T1+1,]+Y_trend[i,T1+1,]
        
        for(t in (T1+2):T2){ # Burst 2
          Y_trend[i,t,1:2] = person_params[i,3:4]
          
          mean_dev = c(
            AR1[i] * Y_dev[i,t-1,1] + CR21[i] * Y_dev[i,t-1,2],# reading
            AR2[i] * Y_dev[i,t-1,2] + CR12[i] * Y_dev[i,t-1,1] # math
          )
          
          Y_dev[i,t,] = rmvnorm(1, mean=mean_dev, sigma=Sigma_innovation)
          
          Y[i,t,] = Y_trend[i,t,] + Y_dev[i,t,]
          
        } #Closing Burst2 loop
        
        
        #2. Generate first burst3 observations
        Y_trend[i,T2+1,1:2] = person_params[i,5:6] #Intercept1_Burst1, Intercept2_Burst1
        Y_dev[i,T2+1,] = rmvnorm(1, mean=c(0,0), sigma=Sigma_innovation)
        Y[i,T2+1,] = Y_dev[i,T2+1,]+Y_trend[i,T2+1,]
        
        for(t in (T2+2):T3){ # Burst 3
          Y_trend[i,t,1:2] = person_params[i,5:6]
          
          mean_dev = c(
            AR1[i] * Y_dev[i,t-1,1] + CR21[i] * Y_dev[i,t-1,2],# reading
            AR2[i] * Y_dev[i,t-1,2] + CR12[i] * Y_dev[i,t-1,1] # math
          )
          
          Y_dev[i,t,] = rmvnorm(1, mean=mean_dev, sigma=Sigma_innovation)
          
          Y[i,t,] = Y_trend[i,t,] + Y_dev[i,t,]
          
        } #Closing Burst3 loop
        
      } #Closing person loop
      
      ######################## Save Data ########################
      # Reshape to long format
      dat_long <- matrix(NA, nrow=nP*nT, ncol=4)  # [id, time, y1, y2]
      colnames(dat_long) <- c("id", "time", "y1", "y2")
      
      # Fill data
      for(i in 1:nP){
        for(t in 1:nT){
          row_idx = (i-1)*nT + t
          dat_long[row_idx, "id"] = i
          dat_long[row_idx, "time"] = time[t]
          dat_long[row_idx, "y1"] = Y[i,t,1]
          dat_long[row_idx, "y2"] = Y[i,t,2]
        }
      }
      
      dat_long = as.data.frame(dat_long)
      
      # Save data as CSV
      write.csv(dat_long, paste0(work_path, "/", dat_name, ".csv"), 
                row.names=FALSE)
      
      
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