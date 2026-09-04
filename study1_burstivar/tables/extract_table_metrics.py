# ----------------------------------------------------------------------------
# PUBLIC-REPO COPY. Original: IP_1b/paper_simulation_archive/Study1_BurstiVAR/extract_table_metrics.py
# Modified for this repository: the two input directories (R2605, RCRE) now
# point to ..\results. Run it from this tables directory (the script uses
# Windows path separators, as the original did). All other content is
# identical to the version that produced the published table values.
# ----------------------------------------------------------------------------
#!/usr/bin/env python3
# Extract per-cell metrics for the 11 dynamic-process parameters shown in
# Study 1 Table 2, to build (a) the modified Table 2 (add Mean est) and
# (b) the new appendix table (RMSE, posterior SD, empirical SD, RDSE).
# CSV cols: 0 name,1 TruePar,2 MeanPar_hat,3 RMSE,4 rBias,5 SD,
#           6 Mean_SE_hat,7 RDSE,8 LL,9 UL,10 Coverage,11 Power,12 TypeIError
import csv

R2605 = r"..\results"  # REPO EDIT: was r"D:\XXYDATAanalysis\IP_1b\MLGVAR2605\result"; run from the tables directory
RCRE  = r"..\results"  # REPO EDIT: was the CREqual0 result directory; run from the tables directory
cells = [
    ("Tb3_N100",  Rf"{R2605}\mlGVARNoCorNoME_BurstIntercept_3Burst_3T_MCfileSumm_nT9_nP100.csv"),
    ("Tb3_N500",  Rf"{R2605}\mlGVARNoCorNoME_BurstIntercept_3Burst_3T_MCfileSumm_nT9_nP500.csv"),
    ("Tb5_N100",  Rf"{RCRE}\mlGVARNoCorNoME_BurstIntercept_CREqual0_MCfileSumm_nT15_nP100.csv"),
    ("Tb5_N500",  Rf"{RCRE}\mlGVARNoCorNoME_BurstIntercept_CREqual0_MCfileSumm_nT15_nP500.csv"),
    ("Tb20_N100", Rf"{RCRE}\mlGVARNoCorNoME_BurstIntercept_CREqual0_MCfileSumm_nT60_nP100.csv"),
    ("Tb20_N500", Rf"{RCRE}\mlGVARNoCorNoME_BurstIntercept_CREqual0_MCfileSumm_nT60_nP500.csv"),
]
params = [
    ("Level2Mean[7]",  "mu_phi11 AR_y1"),
    ("Level2Mean[8]",  "mu_phi22 AR_y2"),
    ("Level2Mean[9]",  "mu_phi21 CR_y1->y2"),
    ("Level2Mean[10]", "mu_phi12 CR_y2->y1"),
    ("Level2Sigma[7]", "sigma_phi11"),
    ("Level2Sigma[8]", "sigma_phi22"),
    ("Level2Sigma[9]", "sigma_phi21"),
    ("Level2Sigma[10]","sigma_phi12"),
    ("sigma_innovation[1,1]", "sigma2_eps1"),
    ("sigma_innovation[2,2]", "sigma2_eps2"),
    ("sigma_innovation[1,2]", "sigma_eps12"),
]
data = {}
for cname, path in cells:
    d = {}
    with open(path, newline="") as fh:
        rd = csv.reader(fh); next(rd)
        for row in rd:
            if row and row[0]:
                d[row[0]] = row
    data[cname] = d

def f(x, n=3):
    try: return f"{float(x):.{n}f}"
    except Exception: return str(x)

for prow, plabel in params:
    true = data[cells[0][0]][prow][1]
    print(f"\n=== {plabel}   True={true} ===")
    print(f"{'cell':10}{'MeanEst':>9}{'RMSE':>9}{'postSD':>9}{'empSD':>9}{'RDSE':>9} | {'rBias':>8}{'Cov':>5}{'Pwr':>5}{'TypeI':>6}")
    for cname, _ in cells:
        v = data[cname][prow]
        print(f"{cname:10}{f(v[2]):>9}{f(v[3]):>9}{f(v[6]):>9}{f(v[5]):>9}{f(v[7]):>9} | {f(v[4]):>8}{f(v[10],2):>5}{f(v[11],2):>5}{f(v[12],2):>6}")
