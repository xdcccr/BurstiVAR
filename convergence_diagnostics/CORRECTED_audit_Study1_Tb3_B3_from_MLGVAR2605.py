#!/usr/bin/env python3
# ----------------------------------------------------------------------
# PROVENANCE (public BurstiVAR repository)
# Original path: paper_simulation_archive/Convergence_ESS_RHAT/
#   CORRECTED_audit_Study1_Tb3_B3_from_MLGVAR2605.py
# Modified: this header was inserted after the shebang, and the absolute
#   BASE path (D:\XXYDATAanalysis\IP_1b\MLGVAR2605\result) was replaced
#   with "./MLGVAR2605/result" so the script runs from this directory
#   (run from this directory; ./MLGVAR2605/result must hold your local
#   per-replication resulttables for the B=3, Tb=3 re-run). All other
#   content is identical to the version that produced the published
#   results.
# ----------------------------------------------------------------------
# ============================================================================
# CORRECTED convergence audit for Study 1, Tb = 3 (B = 3, T = 9).
#
# WHY THIS EXISTS:
#   The original R audit (MLGVAR2601/Extract_ESS_RHAT_*.R) computed its
#   "Study1_Tb3" rows from prefix "..._2Burst_3T", nT=6 -- i.e. the B=2, T=6
#   TWO-BURST MINIMUM DESIGN (Supplement D), NOT the B=3, T=9 Study 1 Tb=3
#   cell. So the archive's Study1_Tb3 ESS/RHAT summaries describe the wrong
#   condition. The real Study 1 Tb=3 per-replication resulttables live in
#   MLGVAR2605 (nT9). This script audits those, so the numbers can be compared
#   directly against the manuscript's Appendix B table (tab:convergence_summary).
#
# resulttable column layout (from summarizePost):
#   0:param  1:mean 2:PSD 3:PCI2.5 4:PCI97.5 5:HDI_Low 6:HDI_High
#   7:st-.05 8:(-.05,.05) 9:lt.05  10:ESS  11:RHAT
#
# Mirrors the per-condition "overview" that the R script produced.
# ============================================================================
import csv, os

BASE = "./MLGVAR2605/result"
PREFIX = "mlGVARNoCorNoME_BurstIntercept_3Burst_3T_resulttable_nT9_nP"
N_REPL = 100
ESS_COL, RHAT_COL = 10, 11

print("Corrected audit: Study 1 Tb=3 (B=3, T=9), source = MLGVAR2605 nT9 resulttables\n")
for nP in (100, 500):
    reps = 0
    n110 = n105 = ness200 = ness100 = 0
    gmin_ess = float("inf"); gmax_rhat = 0.0
    missing = []
    for r in range(1, N_REPL + 1):
        f = os.path.join(BASE, f"{PREFIX}{nP}_r{r}.csv")
        if not os.path.exists(f):
            missing.append(r); continue
        reps += 1
        any110 = any105 = anyess200 = anyess100 = False
        with open(f, newline="") as fh:
            rd = csv.reader(fh); next(rd, None)
            for row in rd:
                if len(row) <= RHAT_COL:
                    continue
                try:
                    ess = float(row[ESS_COL]); rhat = float(row[RHAT_COL])
                except ValueError:
                    continue
                gmax_rhat = max(gmax_rhat, rhat); gmin_ess = min(gmin_ess, ess)
                if rhat > 1.10: any110 = True
                if rhat > 1.05: any105 = True
                if ess < 200: anyess200 = True
                if ess < 100: anyess100 = True
        n110 += any110; n105 += any105; ness200 += anyess200; ness100 += anyess100
    print(f"=== Tb3 (B=3, T=9), N={nP}  |  reps read = {reps}, missing = {missing} ===")
    print(f"  % reps any RHAT>1.10 : {n110}/{reps}  ({round(100*n110/reps)}%)")
    print(f"  % reps any RHAT>1.05 : {n105}/{reps}  ({round(100*n105/reps)}%)")
    print(f"  % reps any ESS<200   : {ness200}/{reps}  ({round(100*ness200/reps)}%)")
    print(f"  % reps any ESS<100   : {ness100}/{reps}  ({round(100*ness100/reps)}%)")
    print(f"  Global Min ESS = {gmin_ess:.0f}   Global Max RHAT = {gmax_rhat:.4f}\n")
