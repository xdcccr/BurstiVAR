# ----------------------------------------------------------------------
# PROVENANCE (public BurstiVAR repository)
# Original path: paper_simulation_archive/Convergence_ESS_RHAT/Extract_ESS_RHAT_Summary.R
# Modified: this header was prepended, and the absolute path prefix
#   "D:/XXYDATAanalysis/IP_1b" in the path-configuration lines was replaced
#   with "." so the script runs from this directory (run from this
#   directory; the ./MLGVAR2601/... folders must hold your local
#   per-replication simulation output). All other content is identical to
#   the version that produced the published results.
# CAUTION: the "Study1_Tb3" block below points at a result folder that did
#   not exist, so Round 1 produced no Tb=3 output; and that block targets
#   the B=2 two-burst design, not the B=3 Study 1 Tb=3 cell. See README.md.
# ----------------------------------------------------------------------

########################################################################
# Extract_ESS_RHAT_Summary.R
# 
# Purpose: Extract ESS and RHAT from individual replication result tables
#          for Study 1 and Study 2, and produce summary statistics.
#
# Output per condition:
#   1. A CSV with per-replication, per-parameter ESS and RHAT
#   2. A summary CSV with min/median/mean/max ESS and max RHAT per parameter
#   3. Count of replications with any RHAT > 1.10 or any RHAT > 1.05
#   4. Count of replications with any ESS < 200 or any ESS < 100
########################################################################

rm(list = ls())

required_packages <- c("data.table")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# =====================================================================
# USER SETTINGS — modify paths as needed
# =====================================================================

# Study 1 base path
study1_base <- "./MLGVAR2601/mlGVARNoCorNoME_BurstIntercept_CREqual0"
study1_result <- paste0(study1_base, "/result")
study1_prefix <- "mlGVARNoCorNoME_BurstIntercept_CREqual0"

# Study 2 base path
study2_base <- "./MLGVAR2601/mlGVARNoCorNoME_GompertzBurst"
study2_result <- paste0(study2_base, "/result")
study2_prefix <- "mlGVARNoCorNoME_GompertzBurst"

# Output directory (will be created if not exists)
output_dir <- "./MLGVAR2601/ESS_RHAT_Audit"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Number of replications
N_repl <- 100

# =====================================================================
# Study 1 conditions
# Condition A: 2 bursts, Tb=3, T=6 (uses _2Burst_3T prefix)
# Condition B: 3 bursts, Tb=5, T=15
# Condition C: 3 bursts, Tb=20, T=60
# =====================================================================

# NOTE: Condition A (Tb=3) uses a DIFFERENT prefix and result directory.
# Adjust if needed. Based on the project files, Condition A files are:
#   mlGVARNoCorNoME_BurstIntercept_2Burst_3T_MCfileSumm_nT6_nP*.csv
# So its individual result tables should have a different prefix.
# If Condition A results are in a different folder, specify below.

study1_condA_base <- "./MLGVAR2601/mlGVARNoCorNoME_BurstIntercept_2Burst_3T"
study1_condA_result <- paste0(study1_condA_base, "/result")
study1_condA_prefix <- "mlGVARNoCorNoME_BurstIntercept_2Burst_3T"

study1_conditions <- list(
  # Condition A: Tb=3 (separate folder/prefix)
  list(name = "Study1_Tb3_N100",  result_path = study1_condA_result, 
       prefix = study1_condA_prefix, nT = 6, nP = 100),
  list(name = "Study1_Tb3_N500",  result_path = study1_condA_result, 
       prefix = study1_condA_prefix, nT = 6, nP = 500),
  # Condition B & C: same folder
  list(name = "Study1_Tb5_N100",  result_path = study1_result, 
       prefix = study1_prefix, nT = 15, nP = 100),
  list(name = "Study1_Tb5_N500",  result_path = study1_result, 
       prefix = study1_prefix, nT = 15, nP = 500),
  list(name = "Study1_Tb20_N100", result_path = study1_result, 
       prefix = study1_prefix, nT = 60, nP = 100),
  list(name = "Study1_Tb20_N500", result_path = study1_result, 
       prefix = study1_prefix, nT = 60, nP = 500)
)

# =====================================================================
# Study 2 conditions
# BurstiVAR fitted on GompertzBurst data — need to check what prefix
# the individual result files use. There may be TWO sets of results:
#   (a) BurstiVAR fitting (5Burst model on Gompertz data)
#   (b) GoBurstiVAR fitting (true model)
# Adjust prefixes accordingly.
# =====================================================================

# For the GoBurstiVAR TRUE MODEL fits:
study2_conditions_true <- list(
  list(name = "Study2_GoBurstiVAR_Tb5_N100",  result_path = study2_result, 
       prefix = "mlGVARNoCorNoME_GompertzBurst", nT = 25, nP = 100),
  list(name = "Study2_GoBurstiVAR_Tb5_N500",  result_path = study2_result, 
       prefix = "mlGVARNoCorNoME_GompertzBurst", nT = 25, nP = 500),
  list(name = "Study2_GoBurstiVAR_Tb20_N100", result_path = study2_result, 
       prefix = "mlGVARNoCorNoME_GompertzBurst", nT = 100, nP = 100)
)

# For the BurstiVAR (5Burst unrestricted) fits on Gompertz data:
# These likely have a different prefix — check your folder.
# Based on project CSV names: 
#   mlGVARNoCorNoME_BurstIntercept_5Burst_on_mlGVARNoCorNoME_GompertzBurst_MCfileSumm_*.csv
# So the individual result files probably use a similar naming convention.
study2_bursti_base <- "./MLGVAR2601/mlGVARNoCorNoME_BurstIntercept_5Burst_on_mlGVARNoCorNoME_GompertzBurst"
study2_bursti_result <- paste0(study2_bursti_base, "/result")
study2_bursti_prefix <- "mlGVARNoCorNoME_BurstIntercept_5Burst_on_mlGVARNoCorNoME_GompertzBurst"

study2_conditions_bursti <- list(
  list(name = "Study2_BurstiVAR_Tb5_N100",  result_path = study2_bursti_result, 
       prefix = study2_bursti_prefix, nT = 25, nP = 100),
  list(name = "Study2_BurstiVAR_Tb5_N500",  result_path = study2_bursti_result, 
       prefix = study2_bursti_prefix, nT = 25, nP = 500),
  list(name = "Study2_BurstiVAR_Tb20_N100", result_path = study2_bursti_result, 
       prefix = study2_bursti_prefix, nT = 100, nP = 100)
)

# Combine all conditions
all_conditions <- c(study1_conditions, study2_conditions_true, study2_conditions_bursti)

# =====================================================================
# CORE FUNCTION: Extract ESS and RHAT from one condition
# =====================================================================

extract_ess_rhat <- function(cond, N_repl = 100) {
  
  cat(sprintf("\n========== Processing: %s ==========\n", cond$name))
  cat(sprintf("  Path: %s\n", cond$result_path))
  cat(sprintf("  Prefix: %s, nT=%d, nP=%d\n", cond$prefix, cond$nT, cond$nP))
  
  all_records <- list()
  files_found <- 0
  files_missing <- 0
  
  for (r in 1:N_repl) {
    file_name <- paste0(cond$prefix, "_resulttable_nT", cond$nT, 
                        "_nP", cond$nP, "_r", r, ".csv")
    full_path <- file.path(cond$result_path, file_name)
    
    if (file.exists(full_path)) {
      files_found <- files_found + 1
      
      res <- tryCatch({
        read.csv(full_path, row.names = 1, stringsAsFactors = FALSE)
      }, error = function(e) {
        cat(sprintf("  WARNING: Error reading r=%d: %s\n", r, e$message))
        return(NULL)
      })
      
      if (is.null(res)) next
      
      # Find ESS and RHAT columns
      cnames <- colnames(res)
      ess_col <- grep("ESS", cnames, ignore.case = TRUE)
      rhat_col <- grep("RHAT", cnames, ignore.case = TRUE)
      
      if (length(ess_col) == 0 || length(rhat_col) == 0) {
        cat(sprintf("  WARNING: r=%d missing ESS/RHAT columns. Cols: %s\n", 
                    r, paste(cnames, collapse = ", ")))
        next
      }
      
      ess_col <- ess_col[1]
      rhat_col <- rhat_col[1]
      
      for (par in rownames(res)) {
        ess_val <- as.numeric(res[par, ess_col])
        rhat_val <- as.numeric(res[par, rhat_col])
        
        all_records[[length(all_records) + 1]] <- data.frame(
          Condition = cond$name,
          Replication = r,
          Parameter = par,
          ESS = ess_val,
          RHAT = rhat_val,
          stringsAsFactors = FALSE
        )
      }
    } else {
      files_missing <- files_missing + 1
    }
  }
  
  cat(sprintf("  Files found: %d, Files missing: %d\n", files_found, files_missing))
  
  if (length(all_records) == 0) {
    cat("  NO DATA EXTRACTED. Check file naming convention.\n")
    # Try to list what's actually in the directory
    if (dir.exists(cond$result_path)) {
      sample_files <- head(list.files(cond$result_path, pattern = "\\.csv$"), 5)
      cat("  Sample files in directory:\n")
      for (f in sample_files) cat(sprintf("    %s\n", f))
    } else {
      cat("  Directory does not exist!\n")
    }
    return(NULL)
  }
  
  df <- rbindlist(all_records)
  return(df)
}

# =====================================================================
# MAIN LOOP: Process all conditions
# =====================================================================

grand_summary <- list()

for (cond in all_conditions) {
  
  df <- extract_ess_rhat(cond, N_repl)
  
  if (is.null(df) || nrow(df) == 0) {
    cat(sprintf("  Skipping %s (no data)\n", cond$name))
    next
  }
  
  # --- Save raw per-replication data ---
  raw_file <- file.path(output_dir, paste0(cond$name, "_ESS_RHAT_raw.csv"))
  fwrite(df, raw_file)
  cat(sprintf("  Raw data saved: %s\n", raw_file))
  
  # --- Per-parameter summary ---
  param_summ <- df[, .(
    Min_ESS    = min(ESS, na.rm = TRUE),
    Q25_ESS    = quantile(ESS, 0.25, na.rm = TRUE),
    Median_ESS = median(ESS, na.rm = TRUE),
    Mean_ESS   = round(mean(ESS, na.rm = TRUE), 1),
    Max_ESS    = max(ESS, na.rm = TRUE),
    Min_RHAT   = round(min(RHAT, na.rm = TRUE), 4),
    Median_RHAT = round(median(RHAT, na.rm = TRUE), 4),
    Max_RHAT   = round(max(RHAT, na.rm = TRUE), 4),
    N_RHAT_gt_1.10 = sum(RHAT > 1.10, na.rm = TRUE),
    N_RHAT_gt_1.05 = sum(RHAT > 1.05, na.rm = TRUE),
    N_ESS_lt_200   = sum(ESS < 200, na.rm = TRUE),
    N_ESS_lt_100   = sum(ESS < 100, na.rm = TRUE),
    N_reps = .N
  ), by = Parameter]
  
  summ_file <- file.path(output_dir, paste0(cond$name, "_ESS_RHAT_summary.csv"))
  fwrite(param_summ, summ_file)
  cat(sprintf("  Summary saved: %s\n", summ_file))
  
  # --- Replication-level flags ---
  rep_flags <- df[, .(
    Max_RHAT = max(RHAT, na.rm = TRUE),
    Min_ESS  = min(ESS, na.rm = TRUE),
    Any_RHAT_gt_1.10 = any(RHAT > 1.10, na.rm = TRUE),
    Any_RHAT_gt_1.05 = any(RHAT > 1.05, na.rm = TRUE),
    Any_ESS_lt_200   = any(ESS < 200, na.rm = TRUE),
    Any_ESS_lt_100   = any(ESS < 100, na.rm = TRUE)
  ), by = Replication]
  
  n_total_reps <- nrow(rep_flags)
  
  overview <- data.frame(
    Condition = cond$name,
    N_reps = n_total_reps,
    N_nonconverged_Rhat1.10 = sum(rep_flags$Any_RHAT_gt_1.10),
    Pct_nonconverged_Rhat1.10 = round(100 * sum(rep_flags$Any_RHAT_gt_1.10) / n_total_reps, 1),
    N_nonconverged_Rhat1.05 = sum(rep_flags$Any_RHAT_gt_1.05),
    Pct_nonconverged_Rhat1.05 = round(100 * sum(rep_flags$Any_RHAT_gt_1.05) / n_total_reps, 1),
    N_anyESS_lt200 = sum(rep_flags$Any_ESS_lt_200),
    Pct_anyESS_lt200 = round(100 * sum(rep_flags$Any_ESS_lt_200) / n_total_reps, 1),
    N_anyESS_lt100 = sum(rep_flags$Any_ESS_lt_100),
    Pct_anyESS_lt100 = round(100 * sum(rep_flags$Any_ESS_lt_100) / n_total_reps, 1),
    Global_Min_ESS = min(rep_flags$Min_ESS),
    Global_Max_RHAT = round(max(rep_flags$Max_RHAT), 4),
    stringsAsFactors = FALSE
  )
  
  grand_summary[[length(grand_summary) + 1]] <- overview
  
  cat(sprintf("  === OVERVIEW ===\n"))
  cat(sprintf("  Reps with any RHAT > 1.10: %d / %d (%.1f%%)\n", 
              overview$N_nonconverged_Rhat1.10, n_total_reps, overview$Pct_nonconverged_Rhat1.10))
  cat(sprintf("  Reps with any RHAT > 1.05: %d / %d (%.1f%%)\n",
              overview$N_nonconverged_Rhat1.05, n_total_reps, overview$Pct_nonconverged_Rhat1.05))
  cat(sprintf("  Reps with any ESS < 200:   %d / %d (%.1f%%)\n",
              overview$N_anyESS_lt200, n_total_reps, overview$Pct_anyESS_lt200))
  cat(sprintf("  Reps with any ESS < 100:   %d / %d (%.1f%%)\n",
              overview$N_anyESS_lt100, n_total_reps, overview$Pct_anyESS_lt100))
  cat(sprintf("  Global Min ESS: %d | Global Max RHAT: %.4f\n",
              overview$Global_Min_ESS, overview$Global_Max_RHAT))
}

# =====================================================================
# SAVE GRAND SUMMARY
# =====================================================================

if (length(grand_summary) > 0) {
  grand_df <- rbindlist(grand_summary)
  grand_file <- file.path(output_dir, "GRAND_SUMMARY_ESS_RHAT.csv")
  fwrite(grand_df, grand_file)
  cat(sprintf("\n\n====== GRAND SUMMARY saved to: %s ======\n", grand_file))
  print(grand_df)
} else {
  cat("\nNo conditions were successfully processed.\n")
}

cat("\n\nDone!\n")
