# ----------------------------------------------------------------------------
# PUBLIC-REPO COPY. Original: IP_1b/paper_simulation_archive/Study1_BurstiVAR/gen_table_latex.py
# Modified for this repository: the two archive input directories (R3, R520) now point to ../results.
# All other content is identical to the version that produced the published
# table blocks.
# ----------------------------------------------------------------------------
#!/usr/bin/env python3
# Regenerate the two Study 1 tables after the 2026-07-16 restructure
# (advisor comment: RMSE is a point-estimate metric -> main table;
#  power/Type-I are inferential -> appendix table):
#   Table 2 (tab:study1_results):       Est Mean, rBias, RMSE
#   Appendix (tab:study1_supp_metrics): Cov, Pwr/T-I, Emp SD, RDSE
# Emits FULL table blocks (\begin{table} -> \end{table}) to
#   table2_restructured.tex  and  suppmetrics_restructured.tex
# in this folder, for review and range-replacement into overleaf/main.tex.
# Reads the ARCHIVE MCfileSumm copies (canonical per archive README).
# The pre-restructure generator is kept as gen_table_latex_v1_pre20260716.py.
# 2026-09-01 (coauthor IJAP comment 3): Pwr cells for the between-person SDs
# are now "---" (SDs are positive by construction, so CI-excludes-zero power
# is vacuous); both supp notes updated accordingly, and their "whose" wording
# synced to the "for which the" phrasing already applied in overleaf/main.tex.
# CSV cols: 0 name,1 True,2 MeanPar_hat,3 RMSE,4 rBias,5 SD,6 Mean_SE_hat,
#           7 RDSE,8 LL,9 UL,10 Cov,11 Power,12 TypeIError
import csv, os
HERE = os.path.dirname(os.path.abspath(__file__))
R3  = os.path.join(HERE, "..", "results")  # REPO EDIT: was os.path.join(HERE, "MCfiles_Tb3_from_MLGVAR2605")
R520 = os.path.join(HERE, "..", "results")  # REPO EDIT: was os.path.join(HERE, "MCfiles_Tb5_Tb20_from_CREqual0")
paths = {
 "Tb3_N100":  os.path.join(R3,  "mlGVARNoCorNoME_BurstIntercept_3Burst_3T_MCfileSumm_nT9_nP100.csv"),
 "Tb3_N500":  os.path.join(R3,  "mlGVARNoCorNoME_BurstIntercept_3Burst_3T_MCfileSumm_nT9_nP500.csv"),
 "Tb5_N100":  os.path.join(R520, "mlGVARNoCorNoME_BurstIntercept_CREqual0_MCfileSumm_nT15_nP100.csv"),
 "Tb5_N500":  os.path.join(R520, "mlGVARNoCorNoME_BurstIntercept_CREqual0_MCfileSumm_nT15_nP500.csv"),
 "Tb20_N100": os.path.join(R520, "mlGVARNoCorNoME_BurstIntercept_CREqual0_MCfileSumm_nT60_nP100.csv"),
 "Tb20_N500": os.path.join(R520, "mlGVARNoCorNoME_BurstIntercept_CREqual0_MCfileSumm_nT60_nP500.csv"),
}
D = {}
for k, p in paths.items():
    d = {}
    with open(p, newline="") as fh:
        rd = csv.reader(fh); next(rd)
        for r in rd:
            if r and r[0]: d[r[0]] = r
    D[k] = d

def lead(v, n):
    # APA-style leading-decimal formatting; keeps the sign even when the
    # rounded magnitude is .000 (matches the published table convention).
    try: x = float(v)
    except: return str(v)
    s = f"{abs(x):.{n}f}"
    if float(s) < 1: s = "." + s.split(".")[1]
    return f"$-{s}$" if x < 0 else s

groups = [
 ("\\quad \\textit{Fixed Effects}", [
   ("Level2Mean[7]",  r"$\mu_{\phi_{11}}$ (AR: $y_1$)",          "0.30"),
   ("Level2Mean[8]",  r"$\mu_{\phi_{22}}$ (AR: $y_2$)",          "0.20"),
   ("Level2Mean[9]",  r"$\mu_{\phi_{21}}$ (CR: $y_1 \to y_2$)",  "$-0.15$"),
   ("Level2Mean[10]", r"$\mu_{\phi_{12}}$ (CR: $y_2 \to y_1$)",  "0"),
 ]),
 ("\\quad \\textit{Random Effect SDs}", [
   ("Level2Sigma[7]", r"$\sigma_{\phi_{11}}$", "0.10"),
   ("Level2Sigma[8]", r"$\sigma_{\phi_{22}}$", "0.10"),
   ("Level2Sigma[9]", r"$\sigma_{\phi_{21}}$", "0.10"),
   ("Level2Sigma[10]",r"$\sigma_{\phi_{12}}$", "0.10"),
 ]),
 ("\\quad \\textit{Innovation Covariance}", [
   ("sigma_innovation[1,1]", r"$\sigma^2_{\varepsilon_1}$",  "1.0"),
   ("sigma_innovation[2,2]", r"$\sigma^2_{\varepsilon_2}$",  "1.0"),
   ("sigma_innovation[1,2]", r"$\sigma_{\varepsilon_{12}}$", "0.3"),
 ]),
]
TB = ["Tb3", "Tb5", "Tb20"]

def body(ncols, cellfn, grps):
    out = []
    for Nlab, Nsuf in [("$N = 100$", "N100"), ("$N = 500$", "N500")]:
        if Nsuf == "N500":
            out.append("\\midrule")
        out.append(f"\\multicolumn{{{ncols}}}{{l}}{{\\textit{{{Nlab}}}}} \\\\")
        out.append("\\addlinespace")
        for gi, (gtitle, plist) in enumerate(grps):
            out.append(f"\\multicolumn{{{ncols}}}{{l}}{{{gtitle}}} \\\\")
            for row, label, true in plist:
                parts = [cellfn(D[f"{tb}_{Nsuf}"][row], row, f"{tb}_{Nsuf}") for tb in TB]
                out.append(f"{label} & {true} & & " + " & ".join(parts) + r" \\")
            out.append("\\addlinespace")
    if out and out[-1] == "\\addlinespace":
        out.pop()
    return "\n".join(out)

# ---- Table 2: Est Mean, rBias, RMSE (point-estimate metrics) ----
def t2cell(v, row, cond):
    est = lead(v[2], 3)
    rb = "---" if row == "Level2Mean[10]" else lead(v[4], 3)   # rBias undefined for true value 0
    rmse = lead(v[3], 3)
    return f"{est} & {rb} & {rmse}"

TABLE2 = r"""\begin{table}[htbp]
\centering
\caption{Study 1: Monte Carlo Simulation Results for Dynamic Process Parameters}
\label{tab:study1_results}
\begin{threeparttable}
\small
\setlength{\tabcolsep}{3.5pt}
\begin{tabular}{llrrrrrrrrrrr}
\toprule
& & & \multicolumn{3}{c}{$T_b = 3$} & \multicolumn{3}{c}{$T_b = 5$} & \multicolumn{3}{c}{$T_b = 20$} \\
& & & \multicolumn{3}{c}{(Total $T = 9$)} & \multicolumn{3}{c}{(Total $T = 15$)} & \multicolumn{3}{c}{(Total $T = 60$)} \\
\cmidrule(lr){4-6} \cmidrule(lr){7-9} \cmidrule(lr){10-12}
Parameter & True & & \shortstack{Est\\Mean} & rBias & RMSE & \shortstack{Est\\Mean} & rBias & RMSE & \shortstack{Est\\Mean} & rBias & RMSE \\
\midrule
""" + body(12, t2cell, groups) + r"""
\bottomrule
\end{tabular}
\begin{tablenotes}[flushleft]
\small
\item \textit{Note.} Est Mean = mean posterior estimate across replications; rBias = relative bias (bias/true value), undefined (---) for the null cross-regressive effect $\mu_{\phi_{12}}$ (true value $0$); RMSE = root mean square error. Each cell summarizes 100 replications; $B = 3$ in all conditions. Coverage, power and Type-I error rates, and standard-error accuracy for the same parameters are reported in Table~\ref{tab:study1_supp_metrics}; per-condition $\hat{R}$ and ESS distributions appear in Appendix~\ref{sec:supplement_b}.
\end{tablenotes}
\end{threeparttable}
\end{table}"""

# ---- Appendix: Cov, Pwr/T-I, Emp SD, RDSE (interval and SE metrics) ----
def apcell(v, row, cond):
    cov = lead(v[10], 2)
    if row == "Level2Mean[10]":                 # null CR: Type-I error
        pw = lead(v[12], 2)
    elif row.startswith("sigma_innovation"):    # not reported; see table note
        pw = ""
    elif row.startswith("Level2Sigma"):         # positive by construction; no power (IJAP c3)
        pw = "---"
    else:
        pw = lead(v[11], 2)
    return f"{cov} & {pw} & {lead(v[5], 3)} & {lead(v[7], 2)}"

SUPP = r"""\begin{table}[htbp]
\centering
\caption{Study 1: Supplementary Recovery Metrics for Dynamic Process Parameters}
\label{tab:study1_supp_metrics}
\begin{threeparttable}
\footnotesize
\setlength{\tabcolsep}{3pt}
\begin{tabular}{llrrrrrrrrrrrrr}
\toprule
& & & \multicolumn{4}{c}{$T_b = 3$ ($T = 9$)} & \multicolumn{4}{c}{$T_b = 5$ ($T = 15$)} & \multicolumn{4}{c}{$T_b = 20$ ($T = 60$)} \\
\cmidrule(lr){4-7} \cmidrule(lr){8-11} \cmidrule(lr){12-15}
Parameter & True & & Cov & \shortstack{Pwr/\\T-I} & \shortstack{Emp\\SD} & RDSE & Cov & \shortstack{Pwr/\\T-I} & \shortstack{Emp\\SD} & RDSE & Cov & \shortstack{Pwr/\\T-I} & \shortstack{Emp\\SD} & RDSE \\
\midrule
""" + body(15, apcell, groups) + r"""
\bottomrule
\end{tabular}
\begin{tablenotes}[flushleft]
\footnotesize
\item \textit{Note.} Companion to Table~\ref{tab:study1_results}. Cov = coverage of 95\% credible intervals (nominal $.95$); Pwr/T-I = power for the nonzero fixed effects (proportion of replications for which the 95\% credible interval excluded zero), or the Type-I error rate for the null cross-regressive effect $\mu_{\phi_{12}}$ (true value $0$); power is not reported for the between-person SDs (---) or the innovation parameters, as variance components are positive by construction; Emp SD = empirical SD of the point estimates across replications (the estimator's Monte Carlo sampling variability); RDSE = relative difference in SEs, $(\overline{\mathrm{SE}}_{\text{post}} - \mathrm{Emp\ SD})/\mathrm{Emp\ SD}$, where $\overline{\mathrm{SE}}_{\text{post}}$ is the mean posterior SD; values near $0$ mean the posterior SD matches the estimator's true variability, positive values mean it overstates it. Each cell summarizes 100 replications; $B = 3$ in all conditions.
\end{tablenotes}
\end{threeparttable}
\end{table}"""

# =====================================================================
# 2026-07-16 (Step 3): trend parameters move into the MAIN TEXT as the new
# Table 3 (Est Mean, rBias, RMSE — Table 2 style), with their inferential
# metrics (Cov, Pwr/T-I, Emp SD, RDSE) in a second appendix table. These
# supersede tab:study1_intercept_sd (appendix section deleted) and the
# never-integrated Study-1 tables in overleaf/appendix_trend_parameters.tex.
# =====================================================================
trend_groups = [
 ("\\quad \\textit{Burst-Intercept Means}", [
   ("Level2Mean[1]", r"$\mu_{\alpha_{11}}$ ($y_1$, burst 1)", "0"),
   ("Level2Mean[2]", r"$\mu_{\alpha_{21}}$ ($y_2$, burst 1)", "1"),
   ("Level2Mean[3]", r"$\mu_{\alpha_{12}}$ ($y_1$, burst 2)", "1"),
   ("Level2Mean[4]", r"$\mu_{\alpha_{22}}$ ($y_2$, burst 2)", "1.5"),
   ("Level2Mean[5]", r"$\mu_{\alpha_{13}}$ ($y_1$, burst 3)", "2"),
   ("Level2Mean[6]", r"$\mu_{\alpha_{23}}$ ($y_2$, burst 3)", "2"),
 ]),
 ("\\quad \\textit{Burst-Intercept Between-Person SDs}", [
   ("Level2Sigma[1]", r"$\sigma_{\alpha_{11}}$ ($y_1$, burst 1)", "1.0"),
   ("Level2Sigma[2]", r"$\sigma_{\alpha_{21}}$ ($y_2$, burst 1)", "1.0"),
   ("Level2Sigma[3]", r"$\sigma_{\alpha_{12}}$ ($y_1$, burst 2)", "1.5"),
   ("Level2Sigma[4]", r"$\sigma_{\alpha_{22}}$ ($y_2$, burst 2)", "2.0"),
   ("Level2Sigma[5]", r"$\sigma_{\alpha_{13}}$ ($y_1$, burst 3)", "2.0"),
   ("Level2Sigma[6]", r"$\sigma_{\alpha_{23}}$ ($y_2$, burst 3)", "3.0"),
 ]),
]
# the Tb=3 breakdown cells (both N): bold their coverage in the supp table,
# matching the convention of the superseded tab:study1_intercept_sd
COLLAPSE = {("Tb3_N100", "Level2Sigma[4]"), ("Tb3_N100", "Level2Sigma[6]"),
            ("Tb3_N500", "Level2Sigma[4]"), ("Tb3_N500", "Level2Sigma[6]")}

def t3cell(v, row, cond):
    est = lead(v[2], 3)
    rb = "---" if row == "Level2Mean[1]" else lead(v[4], 3)   # rBias undefined for true value 0
    return f"{est} & {rb} & {lead(v[3], 3)}"

TABLE3 = r"""\begin{table}[htbp]
\centering
\caption{Study 1: Monte Carlo Simulation Results for Burst-Related (Trend) Parameters}
\label{tab:study1_trend}
\begin{threeparttable}
\small
\setlength{\tabcolsep}{3.5pt}
\begin{tabular}{llrrrrrrrrrrr}
\toprule
& & & \multicolumn{3}{c}{$T_b = 3$} & \multicolumn{3}{c}{$T_b = 5$} & \multicolumn{3}{c}{$T_b = 20$} \\
& & & \multicolumn{3}{c}{(Total $T = 9$)} & \multicolumn{3}{c}{(Total $T = 15$)} & \multicolumn{3}{c}{(Total $T = 60$)} \\
\cmidrule(lr){4-6} \cmidrule(lr){7-9} \cmidrule(lr){10-12}
Parameter & True & & \shortstack{Est\\Mean} & rBias & RMSE & \shortstack{Est\\Mean} & rBias & RMSE & \shortstack{Est\\Mean} & rBias & RMSE \\
\midrule
""" + body(12, t3cell, trend_groups) + r"""
\bottomrule
\end{tabular}
\begin{tablenotes}[flushleft]
\small
\item \textit{Note.} $\mu_{\alpha_{vb}}$ and $\sigma_{\alpha_{vb}}$ denote the population mean and the between-person SD of the burst-$b$ intercept of variable $v$. Est Mean = mean posterior estimate across replications; rBias = relative bias (bias/true value), undefined (---) for $\mu_{\alpha_{11}}$ (true value $0$); RMSE = root mean square error. Each cell summarizes 100 replications; $B = 3$ in all conditions. Coverage, power and Type-I error rates, and standard-error accuracy for the same parameters are reported in Table~\ref{tab:study1_trend_supp}.
\end{tablenotes}
\end{threeparttable}
\end{table}"""

def tscell(v, row, cond):
    cov = lead(v[10], 2)
    if (cond, row) in COLLAPSE:
        cov = f"\\textbf{{{cov}}}"
    if row == "Level2Mean[1]":
        pw = lead(v[12], 2)
    elif row.startswith("Level2Sigma"):         # positive by construction; no power (IJAP c3)
        pw = "---"
    else:
        pw = lead(v[11], 2)
    return f"{cov} & {pw} & {lead(v[5], 3)} & {lead(v[7], 2)}"

TRENDSUPP = r"""\begin{table}[htbp]
\centering
\caption{Study 1: Supplementary Recovery Metrics for Burst-Related (Trend) Parameters}
\label{tab:study1_trend_supp}
\begin{threeparttable}
\footnotesize
\setlength{\tabcolsep}{3pt}
\begin{tabular}{llrrrrrrrrrrrrr}
\toprule
& & & \multicolumn{4}{c}{$T_b = 3$ ($T = 9$)} & \multicolumn{4}{c}{$T_b = 5$ ($T = 15$)} & \multicolumn{4}{c}{$T_b = 20$ ($T = 60$)} \\
\cmidrule(lr){4-7} \cmidrule(lr){8-11} \cmidrule(lr){12-15}
Parameter & True & & Cov & \shortstack{Pwr/\\T-I} & \shortstack{Emp\\SD} & RDSE & Cov & \shortstack{Pwr/\\T-I} & \shortstack{Emp\\SD} & RDSE & Cov & \shortstack{Pwr/\\T-I} & \shortstack{Emp\\SD} & RDSE \\
\midrule
""" + body(15, tscell, trend_groups) + r"""
\bottomrule
\end{tabular}
\begin{tablenotes}[flushleft]
\footnotesize
\item \textit{Note.} Companion to Table~\ref{tab:study1_trend}. $\mu_{\alpha_{vb}}$ and $\sigma_{\alpha_{vb}}$ denote the population mean and the between-person SD of the burst-$b$ intercept of variable $v$. Cov = coverage of 95\% credible intervals (nominal $.95$); Pwr/T-I = power for the nonzero burst-intercept means (proportion of replications for which the 95\% credible interval excluded zero), or the Type-I error rate for $\mu_{\alpha_{11}}$ (true value $0$); power is not reported for the between-person SDs (---), which are positive by construction; Emp SD = empirical SD of the point estimates across replications (the estimator's Monte Carlo sampling variability); RDSE = relative difference in SEs, $(\overline{\mathrm{SE}}_{\text{post}} - \mathrm{Emp\ SD})/\mathrm{Emp\ SD}$, where $\overline{\mathrm{SE}}_{\text{post}}$ is the mean posterior SD. Bold coverage values mark the two largest $y_2$ intercept SDs ($\sigma_{\alpha_{22}} = 2$ and $\sigma_{\alpha_{23}} = 3$), which were underestimated by $25\%$--$35\%$ at $T_b = 3$ with coverage near zero; raising $N$ to $500$ did not help, and the problem was confined to $T_b = 3$. Each cell summarizes 100 replications; $B = 3$ in all conditions.
\end{tablenotes}
\end{threeparttable}
\end{table}"""

for fname, block in [("table2_restructured.tex", TABLE2),
                     ("suppmetrics_restructured.tex", SUPP),
                     ("table3_trend.tex", TABLE3),
                     ("trendsupp.tex", TRENDSUPP)]:
    with open(os.path.join(HERE, fname), "w", encoding="utf-8") as fh:
        fh.write(block + "\n")
    print(f"wrote {fname}  ({block.count(chr(10)) + 1} lines)")
