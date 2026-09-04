# --------------------------------------------------------------------
# BurstiVAR public repository -- Study 2 (GoBurstiVAR).
# Original file: paper_simulation_archive/Study2_Gompertz/Generate_Summ_BurstIntercept_5Burst_TrendRecovery.R
# (originally at the archive root; staged beside the other code).
# Modification for this repository: the default result_path and out_path
# were absolute local paths and are now "./result" and
# "./TrendRecovery_B5"; run this script from this directory, or
# override the paths via the SIM_RESULT / SIM_OUT environment variables.
# The archived outputs it produced are in ../results/TrendRecovery_B5.
# All other content is identical to the version that produced the
# published results.
# --------------------------------------------------------------------

########################################################################
# Generate_Summ_BurstIntercept_5Burst_TrendRecovery.R   (2026-07-16)
#
# B = 5 adaptation of Study2_Gompertz_6Burst/Generate_Summ_BurstIntercept_6Burst.R
# for the advisor's Study 2 comment: report how well BurstiVAR's free burst
# intercepts recover the Gompertz-implied trend.
#
# The BurstiVAR free intercepts have well-defined estimands under the Gompertz
# data-generating process:
#   population mean of burst-b intercept  = E_i[ gompertz(b, theta_i) ]
#   between-person SD of burst-b intercept = SD_i[ gompertz(b, theta_i) ]
# computed by Monte Carlo over the SAME person-level distribution used in
# generation (independent normals; R_full = diag(10) in the generation script).
#
# READS  per-replication resulttables from the ORIGINAL result folder
#        (MLGVAR2601/mlGVARNoCorNoME_GompertzBurst/result), read-only.
# WRITES MCfile + MCfileSumm per condition to TrendRecovery_B5/ in this
#        archive folder (never into the original result folder).
# Also summarizes the VAR rows (Level2Mean/Sigma[11..14]) so the output can be
# cross-checked against the archived original MCfileSumm files.
########################################################################
rm(list=ls())

result_path <- Sys.getenv("SIM_RESULT",
  "./result")
out_path <- Sys.getenv("SIM_OUT",
  "./TrendRecovery_B5")
if(!dir.exists(result_path)) stop("Result directory not found!")
if(!dir.exists(out_path)) dir.create(out_path, recursive=TRUE)

Mean = function(x){ mean(x, na.rm=TRUE) }
SD   = function(x){ sd(x, na.rm=TRUE) }
relBias = function(x, truex){ if(is.na(truex) || abs(truex) < 1e-6) return(NA); mean((x - truex)/truex, na.rm=T) }
RMSE    = function(x, truex){ if(any(is.na(truex))) return(NA); sqrt(Mean((x - truex)^2)) }
find_col_index <- function(col_names, keywords){ for(key in keywords){ idx <- grep(key, col_names, ignore.case=TRUE); if(length(idx)>0) return(idx[1]) }; NA }
gompertz <- function(b, t1, t2, t3){ t1 * exp(-t2 * exp(-b * t3)) }

conds <- list(c(25,100), c(25,500), c(100,100))   # (nT, nP): Tb=5 x N100/N500, Tb=20 x N100
N_repl  <- 100
C_DataGen <- "mlGVARNoCorNoME_GompertzBurst"
C_Fit     <- "mlGVARNoCorNoME_BurstIntercept_5Burst"

# ---- True generating parameters (verified against the archived generation
#      script mlGVARNoCorNoME_GompertzBurst_DataGeneration.R lines 41-57) ----
nBurst <- 5
bvals  <- 1:5                              # integer burst index, as in generation
g_means <- c(10, 3, 1.0,  8, 2.5, 0.8)     # theta1_y1, theta2_y1, theta3_y1, theta1_y2, theta2_y2, theta3_y2
g_sds   <- c(2, 0.5, 0.2, 1.5, 0.4, 0.15)
mu_AR1=0.3; mu_AR2=0.2; mu_CR12=-0.15; mu_CR21=-0.1; sigma_AR=0.1; sigma_CR=0.1
true_sigma_val=1.0^2; true_cov_val=0.3*1.0^2

# ---- Monte Carlo estimands for the burst intercepts ----
set.seed(20260716)
M <- 2e6
draws <- cbind(
  rnorm(M, g_means[1], g_sds[1]), rnorm(M, g_means[2], g_sds[2]), rnorm(M, g_means[3], g_sds[3]),
  rnorm(M, g_means[4], g_sds[4]), rnorm(M, g_means[5], g_sds[5]), rnorm(M, g_means[6], g_sds[6])
)
true_int_mean <- matrix(NA, nBurst, 2)   # [burst, var]
true_int_sd   <- matrix(NA, nBurst, 2)
for(b in 1:nBurst){
  y1 <- gompertz(bvals[b], draws[,1], draws[,2], draws[,3])
  y2 <- gompertz(bvals[b], draws[,4], draws[,5], draws[,6])
  true_int_mean[b,1] <- mean(y1); true_int_sd[b,1] <- sd(y1)
  true_int_mean[b,2] <- mean(y2); true_int_sd[b,2] <- sd(y2)
}
cat("MC-derived TRUE burst-intercept means (rows=burst, cols=y1,y2):\n"); print(round(true_int_mean,4))
cat("MC-derived TRUE burst-intercept SDs:\n"); print(round(true_int_sd,4))
cat("Population-theta plug-in curve (reference, NOT the estimand):\n")
print(round(cbind(gompertz(bvals, g_means[1], g_means[2], g_means[3]),
                  gompertz(bvals, g_means[4], g_means[5], g_means[6])),4))

# ---- True_Values_Map / Param_Labels: interleaved b1y1,b1y2,...,b5y1,b5y2 (per JAGS model) ----
True_Values_Map <- numeric(0); Param_Labels <- character(0)
for(b in 1:nBurst){
  i_y1 <- (b-1)*2 + 1; i_y2 <- (b-1)*2 + 2
  True_Values_Map[paste0("Level2Mean[", i_y1, "]")] <- true_int_mean[b,1]
  True_Values_Map[paste0("Level2Mean[", i_y2, "]")] <- true_int_mean[b,2]
  Param_Labels[paste0("Level2Mean[", i_y1, "]")] <- paste0("Burst", b, "_y1")
  Param_Labels[paste0("Level2Mean[", i_y2, "]")] <- paste0("Burst", b, "_y2")
}
True_Values_Map["Level2Mean[11]"] <- mu_AR1;  Param_Labels["Level2Mean[11]"] <- "AR1"
True_Values_Map["Level2Mean[12]"] <- mu_AR2;  Param_Labels["Level2Mean[12]"] <- "AR2"
True_Values_Map["Level2Mean[13]"] <- mu_CR12; Param_Labels["Level2Mean[13]"] <- "CR_y1toy2"
True_Values_Map["Level2Mean[14]"] <- mu_CR21; Param_Labels["Level2Mean[14]"] <- "CR_y2toy1"
for(b in 1:nBurst){
  i_y1 <- (b-1)*2 + 1; i_y2 <- (b-1)*2 + 2
  True_Values_Map[paste0("Level2Sigma[", i_y1, "]")] <- true_int_sd[b,1]
  True_Values_Map[paste0("Level2Sigma[", i_y2, "]")] <- true_int_sd[b,2]
  Param_Labels[paste0("Level2Sigma[", i_y1, "]")] <- paste0("SD_Burst", b, "_y1")
  Param_Labels[paste0("Level2Sigma[", i_y2, "]")] <- paste0("SD_Burst", b, "_y2")
}
True_Values_Map["Level2Sigma[11]"] <- sigma_AR; Param_Labels["Level2Sigma[11]"] <- "SD_AR1"
True_Values_Map["Level2Sigma[12]"] <- sigma_AR; Param_Labels["Level2Sigma[12]"] <- "SD_AR2"
True_Values_Map["Level2Sigma[13]"] <- sigma_CR; Param_Labels["Level2Sigma[13]"] <- "SD_CR_y1toy2"
True_Values_Map["Level2Sigma[14]"] <- sigma_CR; Param_Labels["Level2Sigma[14]"] <- "SD_CR_y2toy1"
True_Values_Map["sigma_innovation[1,1]"] <- true_sigma_val; Param_Labels["sigma_innovation[1,1]"] <- "InnoVar_y1"
True_Values_Map["sigma_innovation[2,2]"] <- true_sigma_val; Param_Labels["sigma_innovation[2,2]"] <- "InnoVar_y2"
True_Values_Map["sigma_innovation[1,2]"] <- true_cov_val;   Param_Labels["sigma_innovation[1,2]"] <- "InnoCov"
True_Values_Map["sigma_innovation[2,1]"] <- true_cov_val;   Param_Labels["sigma_innovation[2,1]"] <- "InnoCov_sym"

for(cond in conds){
  nT <- cond[1]; nP <- cond[2]
  cat(paste0("\n========== Processing: nT=", nT, ", nP=", nP, " ==========\n"))
  target_params <- names(True_Values_Map)
  summ_cols_base <- c("Tru","Est","se","LL","UL","cFlag","pFlag")
  all_colnames <- c("r", as.vector(outer(summ_cols_base, target_params, paste, sep="_")))
  MCfile <- matrix(NA, nrow=N_repl, ncol=length(all_colnames)); colnames(MCfile) <- all_colnames

  col_map <- list()
  first_file <- list.files(result_path, pattern=paste0(C_Fit, "_resulttable_nT", nT, "_nP", nP, "_r1\\.csv$"), full.names=T)
  if(length(first_file)==0){ cat("No result files for this condition. Skipping.\n"); next }
  temp_df <- read.csv(first_file[1]); cnames <- colnames(temp_df)
  col_map$mean <- find_col_index(cnames, c("Mean","mean","Estimate","Est"))
  col_map$se   <- find_col_index(cnames, c("SD","PSD","se","StdDev"))
  col_map$ll   <- find_col_index(cnames, c("2.5","Low","LL","min"))
  col_map$ul   <- find_col_index(cnames, c("97.5","High","UL","max"))

  valid_reps <- 0
  for(r in 1:N_repl){
    full_path <- file.path(result_path, paste0(C_Fit, "_resulttable_nT", nT, "_nP", nP, "_r", r, ".csv"))
    if(file.exists(full_path)){
      res_table <- read.csv(full_path, row.names=1); rownames(res_table) <- trimws(rownames(res_table))
      MCfile[r,"r"] <- r
      for(par in target_params){ if(par %in% rownames(res_table)){
        tru_val <- True_Values_Map[par]
        get_val <- function(df,row,col_idx){ val <- df[row, col_idx-1]; if(length(val)==0||is.null(val)) return(NA); as.numeric(val) }
        est_val <- get_val(res_table,par,col_map$mean); se_val <- get_val(res_table,par,col_map$se)
        ll_val  <- get_val(res_table,par,col_map$ll);   ul_val <- get_val(res_table,par,col_map$ul)
        if(length(est_val)!=1) est_val<-NA; if(length(se_val)!=1) se_val<-NA
        if(length(ll_val)!=1) ll_val<-NA;   if(length(ul_val)!=1) ul_val<-NA
        cov_flag <- if(!is.na(tru_val) & !is.na(ll_val) & !is.na(ul_val)) ifelse(tru_val>=ll_val & tru_val<=ul_val,1,0) else NA
        power_flag <- NA
        if(!is.na(ll_val) & !is.na(ul_val)){
          is_var <- grepl("Level2Sigma",par) | par=="sigma_innovation[1,1]" | par=="sigma_innovation[2,2]"
          power_flag <- if(is_var) ifelse(ll_val>0,1,0) else ifelse(ll_val>0 | ul_val<0,1,0)
        }
        MCfile[r,paste0("Tru_",par)]<-tru_val; MCfile[r,paste0("Est_",par)]<-est_val
        MCfile[r,paste0("se_",par)]<-se_val;   MCfile[r,paste0("LL_",par)]<-ll_val
        MCfile[r,paste0("UL_",par)]<-ul_val;   MCfile[r,paste0("cFlag_",par)]<-cov_flag
        MCfile[r,paste0("pFlag_",par)]<-power_flag
      }}
      valid_reps <- valid_reps + 1
    }
  }
  cat(paste0("Found ", valid_reps, " valid replication files.\n"))
  if(valid_reps==0){ next }

  write.csv(MCfile, file.path(out_path, paste0(C_Fit, "_TrendRecovery_MCfile_nT", nT, "_nP", nP, ".csv")), row.names=FALSE)

  summ_stat_cols <- c("Label","TruePar","MeanPar_hat","RMSE","rBias","SD","Mean_SE_hat","RDSE",
                      "95%CI_LL","95%CI_UL","Coverage","Power","TypeIError","DetectRate","MissingPer")
  MCfile_summ <- matrix(NA, nrow=length(target_params), ncol=length(summ_stat_cols))
  rownames(MCfile_summ) <- target_params; colnames(MCfile_summ) <- summ_stat_cols
  MC_df <- as.data.frame(MCfile)
  for(par in target_params){
    est_vec <- as.numeric(MC_df[[paste0("Est_",par)]]); se_vec <- as.numeric(MC_df[[paste0("se_",par)]])
    flag_vec <- as.numeric(MC_df[[paste0("cFlag_",par)]]); pflag_vec <- as.numeric(MC_df[[paste0("pFlag_",par)]])
    if(all(is.na(est_vec))) next
    tpv <- True_Values_Map[par]
    MCfile_summ[par,"Label"]<-Param_Labels[par]; MCfile_summ[par,"TruePar"]<-tpv
    MCfile_summ[par,"MeanPar_hat"]<-Mean(est_vec); MCfile_summ[par,"SD"]<-SD(est_vec); MCfile_summ[par,"Mean_SE_hat"]<-Mean(se_vec)
    if(!is.na(tpv)){ MCfile_summ[par,"RMSE"]<-RMSE(est_vec,rep(tpv,length(est_vec))); MCfile_summ[par,"rBias"]<-relBias(Mean(est_vec),tpv) }
    sd_est <- SD(est_vec); if(!is.na(sd_est) && sd_est>0) MCfile_summ[par,"RDSE"]<-(Mean(se_vec)-sd_est)/sd_est
    MCfile_summ[par,"95%CI_LL"]<-quantile(est_vec,0.025,na.rm=T); MCfile_summ[par,"95%CI_UL"]<-quantile(est_vec,0.975,na.rm=T)
    if(!is.na(tpv)){ vf<-sum(!is.na(flag_vec)); if(vf>0) MCfile_summ[par,"Coverage"]<-sum(flag_vec,na.rm=T)/vf }
    vpf<-sum(!is.na(pflag_vec))
    if(vpf>0){ dr<-sum(pflag_vec,na.rm=T)/vpf; MCfile_summ[par,"DetectRate"]<-dr
      if(!is.na(tpv)){ if(abs(tpv)<1e-6) MCfile_summ[par,"TypeIError"]<-dr else MCfile_summ[par,"Power"]<-dr } }
    MCfile_summ[par,"MissingPer"]<-sum(is.na(est_vec))/N_repl
  }
  write.csv(as.data.frame(MCfile_summ), file.path(out_path, paste0(C_Fit, "_TrendRecovery_MCfileSumm_nT", nT, "_nP", nP, ".csv")))
  cat(paste0("Summary saved to TrendRecovery_B5/ for nT", nT, "_nP", nP, "\n"))
}
cat("\n========== BurstIntercept 5-burst trend-recovery summary complete ==========\n")
