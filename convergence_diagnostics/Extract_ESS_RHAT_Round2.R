# ----------------------------------------------------------------------
# PROVENANCE (public BurstiVAR repository)
# Original path: paper_simulation_archive/Convergence_ESS_RHAT/Extract_ESS_RHAT_Round2.R
# Modified: this header was prepended, and the absolute path prefix
#   "D:/XXYDATAanalysis/IP_1b" in the path-configuration lines was replaced
#   with "." so the script runs from this directory (run from this
#   directory; the ./MLGVAR2601/... folders must hold your local
#   per-replication simulation output). All other content is identical to
#   the version that produced the published results.
# CAUTION: the two "Study1_Tb3" conditions below read files with prefix
#   "..._2Burst_3T" (nT=6), i.e. the B=2 two-burst design of Supplement D,
#   NOT the B=3 Study 1 Tb=3 cell. Their output filenames are therefore
#   mislabeled; the repo copies of those outputs are renamed
#   SuppD_TwoBurst_*. See README.md.
# ----------------------------------------------------------------------

########################################################################
# Extract_ESS_RHAT_Round2.R
# 
# Extracts ESS and RHAT for the 5 missing conditions:
#   - Study 1: Tb=3 (nT=6), N=100 and N=500
#   - Study 2: BurstiVAR on Gompertz, Tb=5 N=100, Tb=5 N=500, Tb=20 N=100
########################################################################

rm(list = ls())

required_packages <- c("data.table")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# Output directory
output_dir <- "./MLGVAR2601/ESS_RHAT_Audit"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

N_repl <- 100

# =====================================================================
# Define conditions with CORRECT paths and prefixes
# =====================================================================

all_conditions <- list(
  # Study 1 Tb=3: same result folder as other Study 1, different file prefix
  list(
    name = "Study1_Tb3_N100",
    result_path = "./MLGVAR2601/mlGVARNoCorNoME_BurstIntercept_CREqual0/result",
    prefix = "mlGVARNoCorNoME_BurstIntercept_2Burst_3T",
    nT = 6, nP = 100
  ),
  list(
    name = "Study1_Tb3_N500",
    result_path = "./MLGVAR2601/mlGVARNoCorNoME_BurstIntercept_CREqual0/result",
    prefix = "mlGVARNoCorNoME_BurstIntercept_2Burst_3T",
    nT = 6, nP = 500
  ),
  # Study 2 BurstiVAR (5Burst) on Gompertz data: same folder as GoBurstiVAR results
  list(
    name = "Study2_BurstiVAR_Tb5_N100",
    result_path = "./MLGVAR2601/mlGVARNoCorNoME_GompertzBurst/result",
    prefix = "mlGVARNoCorNoME_BurstIntercept_5Burst",
    nT = 25, nP = 100
  ),
  list(
    name = "Study2_BurstiVAR_Tb5_N500",
    result_path = "./MLGVAR2601/mlGVARNoCorNoME_GompertzBurst/result",
    prefix = "mlGVARNoCorNoME_BurstIntercept_5Burst",
    nT = 25, nP = 500
  ),
  list(
    name = "Study2_BurstiVAR_Tb20_N100",
    result_path = "./MLGVAR2601/mlGVARNoCorNoME_GompertzBurst/result",
    prefix = "mlGVARNoCorNoME_BurstIntercept_5Burst",
    nT = 100, nP = 100
  )
)

# =====================================================================
# CORE FUNCTION
# =====================================================================

extract_ess_rhat <- function(cond, N_repl = 100) {
  
  cat(sprintf("\n========== Processing: %s ==========\n", cond$name))
  cat(sprintf("  Path: %s\n", cond$result_path))
  cat(sprintf("  File pattern: %s_resulttable_nT%d_nP%d_r*.csv\n", 
              cond$prefix, cond$nT, cond$nP))
  
  all_records <- list()
  files_found <- 0
  
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
      
      cnames <- colnames(res)
      ess_col <- grep("ESS", cnames, ignore.case = TRUE)
      rhat_col <- grep("RHAT", cnames, ignore.case = TRUE)
      
      if (length(ess_col) == 0 || length(rhat_col) == 0) {
        cat(sprintf("  WARNING: r=%d missing ESS/RHAT columns.\n", r))
        next
      }
      
      ess_col <- ess_col[1]
      rhat_col <- rhat_col[1]
      
      for (par in rownames(res)) {
        all_records[[length(all_records) + 1]] <- data.frame(
          Condition = cond$name,
          Replication = r,
          Parameter = par,
          ESS = as.numeric(res[par, ess_col]),
          RHAT = as.numeric(res[par, rhat_col]),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  cat(sprintf("  Files found: %d / %d\n", files_found, N_repl))
  
  if (length(all_records) == 0) {
    cat("  NO DATA. Listing sample files in directory:\n")
    if (dir.exists(cond$result_path)) {
      # Show files matching the prefix
      sample_files <- head(list.files(cond$result_path, 
                                       pattern = paste0(cond$prefix, ".*\\.csv$")), 5)
      if (length(sample_files) == 0) {
        sample_files <- head(list.files(cond$result_path, pattern = "\\.csv$"), 5)
      }
      for (f in sample_files) cat(sprintf("    %s\n", f))
    }
    return(NULL)
  }
  
  rbindlist(all_records)
}

# =====================================================================
# MAIN LOOP
# =====================================================================

grand_summary <- list()

for (cond in all_conditions) {
  
  df <- extract_ess_rhat(cond, N_repl)
  if (is.null(df) || nrow(df) == 0) next
  
  # Save raw
  fwrite(df, file.path(output_dir, paste0(cond$name, "_ESS_RHAT_raw.csv")))
  
  # Per-parameter summary
  param_summ <- df[, .(
    Min_ESS = min(ESS, na.rm = TRUE),
    Q25_ESS = quantile(ESS, 0.25, na.rm = TRUE),
    Median_ESS = median(ESS, na.rm = TRUE),
    Mean_ESS = round(mean(ESS, na.rm = TRUE), 1),
    Max_ESS = max(ESS, na.rm = TRUE),
    Min_RHAT = round(min(RHAT, na.rm = TRUE), 4),
    Median_RHAT = round(median(RHAT, na.rm = TRUE), 4),
    Max_RHAT = round(max(RHAT, na.rm = TRUE), 4),
    N_RHAT_gt_1.10 = sum(RHAT > 1.10, na.rm = TRUE),
    N_RHAT_gt_1.05 = sum(RHAT > 1.05, na.rm = TRUE),
    N_ESS_lt_200 = sum(ESS < 200, na.rm = TRUE),
    N_ESS_lt_100 = sum(ESS < 100, na.rm = TRUE),
    N_reps = .N
  ), by = Parameter]
  
  fwrite(param_summ, file.path(output_dir, paste0(cond$name, "_ESS_RHAT_summary.csv")))
  
  # Replication-level flags
  rep_flags <- df[, .(
    Max_RHAT = max(RHAT, na.rm = TRUE),
    Min_ESS = min(ESS, na.rm = TRUE),
    Any_RHAT_gt_1.10 = any(RHAT > 1.10, na.rm = TRUE),
    Any_RHAT_gt_1.05 = any(RHAT > 1.05, na.rm = TRUE),
    Any_ESS_lt_200 = any(ESS < 200, na.rm = TRUE),
    Any_ESS_lt_100 = any(ESS < 100, na.rm = TRUE)
  ), by = Replication]
  
  n_reps <- nrow(rep_flags)
  
  overview <- data.frame(
    Condition = cond$name,
    N_reps = n_reps,
    N_nonconverged_Rhat1.10 = sum(rep_flags$Any_RHAT_gt_1.10),
    Pct_nonconverged_Rhat1.10 = round(100 * sum(rep_flags$Any_RHAT_gt_1.10) / n_reps),
    N_nonconverged_Rhat1.05 = sum(rep_flags$Any_RHAT_gt_1.05),
    Pct_nonconverged_Rhat1.05 = round(100 * sum(rep_flags$Any_RHAT_gt_1.05) / n_reps),
    N_anyESS_lt200 = sum(rep_flags$Any_ESS_lt_200),
    Pct_anyESS_lt200 = round(100 * sum(rep_flags$Any_ESS_lt_200) / n_reps),
    N_anyESS_lt100 = sum(rep_flags$Any_ESS_lt_100),
    Pct_anyESS_lt100 = round(100 * sum(rep_flags$Any_ESS_lt_100) / n_reps),
    Global_Min_ESS = min(rep_flags$Min_ESS),
    Global_Max_RHAT = round(max(rep_flags$Max_RHAT), 4),
    stringsAsFactors = FALSE
  )
  
  grand_summary[[length(grand_summary) + 1]] <- overview
  
  cat(sprintf("  RHAT>1.10: %d/%d | RHAT>1.05: %d/%d | ESS<200: %d/%d | ESS<100: %d/%d\n",
              overview$N_nonconverged_Rhat1.10, n_reps,
              overview$N_nonconverged_Rhat1.05, n_reps,
              overview$N_anyESS_lt200, n_reps,
              overview$N_anyESS_lt100, n_reps))
  cat(sprintf("  Global Min ESS: %d | Global Max RHAT: %.4f\n",
              overview$Global_Min_ESS, overview$Global_Max_RHAT))
}

# Save grand summary for these 5 conditions
if (length(grand_summary) > 0) {
  grand_df <- rbindlist(grand_summary)
  grand_file <- file.path(output_dir, "GRAND_SUMMARY_ESS_RHAT_Round2.csv")
  fwrite(grand_df, grand_file)
  cat(sprintf("\n\nGrand summary saved: %s\n", grand_file))
  print(grand_df)
}

cat("\nDone!\n")
