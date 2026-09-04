# ----------------------------------------------------------------------
# PROVENANCE (public BurstiVAR repository)
# Original path: paper_simulation_archive/Convergence_ESS_RHAT/_verify_2burst_for_suppD.py
# Modified: this header was prepended; the absolute BASE path
#   (D:\XXYDATAanalysis\IP_1b\paper_simulation_archive\Convergence_ESS_RHAT)
#   was replaced with "." so the script runs from this directory; and the
#   input filename pattern Study1_Tb3_{tag}_ESS_RHAT_raw.csv was changed to
#   SuppD_TwoBurst_{tag}_ESS_RHAT_raw.csv to match the renamed repo copies
#   of those (originally mislabeled) files. All other content is identical
#   to the version that produced the published results.
# ----------------------------------------------------------------------

# Verify the Two-Burst (B=2, T=6) convergence numbers quoted in Supplement D.
# Source = the (mislabeled "Study1_Tb3") raw ESS/RHAT files, which the CORRECTED
# audit script confirms are actually the B=2 two-burst design.
# Recompute by the "ANY parameter within a rep crosses the threshold" rule.
import csv, os
BASE = "."

for nP, tag in ((100, "N100"), (500, "N500")):
    f = os.path.join(BASE, f"SuppD_TwoBurst_{tag}_ESS_RHAT_raw.csv")
    # per replication: did ANY param exceed each threshold?
    rep_rhat105 = {}   # rep -> bool any RHAT>1.05
    rep_rhat110 = {}
    rep_ess200 = {}    # rep -> bool any ESS<200
    rep_ess100 = {}
    global_min_ess = float("inf")
    global_max_rhat = 0.0
    params_per_rep = {}
    with open(f, newline="") as fh:
        rd = csv.DictReader(fh)
        for row in rd:
            rep = int(row["Replication"])
            ess = float(row["ESS"]); rhat = float(row["RHAT"])
            params_per_rep[rep] = params_per_rep.get(rep, 0) + 1
            rep_rhat105[rep] = rep_rhat105.get(rep, False) or (rhat > 1.05)
            rep_rhat110[rep] = rep_rhat110.get(rep, False) or (rhat > 1.10)
            rep_ess200[rep]  = rep_ess200.get(rep, False)  or (ess < 200)
            rep_ess100[rep]  = rep_ess100.get(rep, False)  or (ess < 100)
            global_min_ess = min(global_min_ess, ess)
            global_max_rhat = max(global_max_rhat, rhat)
    nrep = len(params_per_rep)
    pct = lambda d: 100.0 * sum(1 for v in d.values() if v) / nrep
    print(f"=== {tag}  (n_reps={nrep}, params/rep={set(params_per_rep.values())}) ===")
    print(f"  % reps any RHAT>1.10 : {pct(rep_rhat110):.0f}%")
    print(f"  % reps any RHAT>1.05 : {pct(rep_rhat105):.0f}%")
    print(f"  % reps any ESS<200   : {pct(rep_ess200):.0f}%")
    print(f"  % reps any ESS<100   : {pct(rep_ess100):.0f}%")
    print(f"  global Min ESS       : {global_min_ess:.0f}")
    print(f"  global Max RHAT      : {global_max_rhat:.4f}")
    print()
