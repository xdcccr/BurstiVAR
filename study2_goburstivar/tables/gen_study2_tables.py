#!/usr/bin/env python3
# --------------------------------------------------------------------
# BurstiVAR public repository -- Study 2 (GoBurstiVAR).
# Original file: paper_simulation_archive/Study2_Gompertz/gen_study2_tables.py
# Modification for this repository: input CSV paths now point to
# ../results/MCfiles and ../results/TrendRecovery_B5; run this script
# from this directory. The six .tex table blocks are written here.
# All other content is identical to the version that produced the
# published results.
# --------------------------------------------------------------------

# Generate the six Study 2 table blocks for the 2026-07 advisor revision:
#   MAIN TEXT (point-estimate metrics, mirroring Study 1's Table 2/3 logic):
#     table4_var.tex       tab:study2_var        rBias + RMSE blocks, BurstiVAR vs GoBurstiVAR
#     table5_trend.tex     tab:study2_trend      NEW: BurstiVAR burst intercepts vs Gompertz-implied
#                                                truths (Est Mean, rBias, RMSE)
#     table6_gompertz.tex  tab:study2_gompertz   moved from Supplement A, restyled to Est/rBias/RMSE
#   APPENDIX (interval/SE metrics):
#     supp2_var.tex        tab:study2_supp_var       Coverage + RDSE blocks, both models
#     supp2_trend.tex      tab:study2_supp_trend     Coverage + RDSE, BurstiVAR intercepts
#     supp2_gompertz.tex   tab:study2_supp_gompertz  Coverage + RDSE, Gompertz parameters
# Sources: MCfiles/ (archived, verified) + TrendRecovery_B5/ (2026-07-16 run,
# pipeline cross-validated against MCfiles to 1e-9 on the shared VAR rows).
import csv, os
HERE = os.path.dirname(os.path.abspath(__file__))

CONDS = ["nT25_nP100", "nT25_nP500", "nT100_nP100"]
def load(path):
    d = {}
    with open(path, newline="") as fh:
        rd = csv.reader(fh); hdr = next(rd)
        for r in rd:
            if r and r[0]:
                d[r[0]] = dict(zip(hdr, r))
    return d
BI = {c: load(os.path.join(HERE, "..", "results", "MCfiles",
      f"mlGVARNoCorNoME_BurstIntercept_5Burst_on_mlGVARNoCorNoME_GompertzBurst_MCfileSumm_{c}.csv")) for c in CONDS}
GO = {c: load(os.path.join(HERE, "..", "results", "MCfiles",
      f"mlGVARNoCorNoME_GompertzBurst_TrueModel_MCfileSumm_{c}.csv")) for c in CONDS}
TR = {c: load(os.path.join(HERE, "..", "results", "TrendRecovery_B5",
      f"mlGVARNoCorNoME_BurstIntercept_5Burst_TrendRecovery_MCfileSumm_{c}.csv")) for c in CONDS}

def lead(v, n):
    # APA leading-decimal formatting, sign kept even when magnitude rounds to 0
    x = float(v)
    s = f"{abs(x):.{n}f}"
    if float(s) < 1: s = "." + s.split(".")[1]
    return f"$-{s}$" if x < 0 else s

# ---- VAR parameters: (label, true, BiVAR key, GoVAR key) ----
VAR_FE = [
 (r"$\mu_{\phi_{11}}$ (AR: $y_1$)",         "0.30",    "Level2Mean[11]", "Level2Mean[7]"),
 (r"$\mu_{\phi_{22}}$ (AR: $y_2$)",         "0.20",    "Level2Mean[12]", "Level2Mean[8]"),
 (r"$\mu_{\phi_{21}}$ (CR: $y_1 \to y_2$)", "$-0.15$", "Level2Mean[13]", "Level2Mean[9]"),
 (r"$\mu_{\phi_{12}}$ (CR: $y_2 \to y_1$)", "$-0.10$", "Level2Mean[14]", "Level2Mean[10]"),
]
VAR_RE = [
 (r"$\sigma_{\phi_{11}}$ (AR: $y_1$)",         "0.10", "Level2Sigma[11]", "Level2Sigma[7]"),
 (r"$\sigma_{\phi_{22}}$ (AR: $y_2$)",         "0.10", "Level2Sigma[12]", "Level2Sigma[8]"),
 (r"$\sigma_{\phi_{21}}$ (CR: $y_1 \to y_2$)", "0.10", "Level2Sigma[13]", "Level2Sigma[9]"),
 (r"$\sigma_{\phi_{12}}$ (CR: $y_2 \to y_1$)", "0.10", "Level2Sigma[14]", "Level2Sigma[10]"),
]
VAR_IN = [
 (r"$\sigma^2_{\varepsilon_1}$",  "1.0", "sigma_innovation[1,1]", "sigma_innovation[1,1]"),
 (r"$\sigma^2_{\varepsilon_2}$",  "1.0", "sigma_innovation[2,2]", "sigma_innovation[2,2]"),
 (r"$\sigma_{\varepsilon_{12}}$", "0.3", "sigma_innovation[1,2]", "sigma_innovation[1,2]"),
]

HEAD2 = r"""& & & \multicolumn{2}{c}{$T_b = 5$, $N = 100$} & & \multicolumn{2}{c}{$T_b = 5$, $N = 500$} & & \multicolumn{2}{c}{$T_b = 20$, $N = 100$} \\
& & & \multicolumn{2}{c}{($T = 25$)} & & \multicolumn{2}{c}{($T = 25$)} & & \multicolumn{2}{c}{($T = 100$)} \\
\cmidrule(lr){4-5} \cmidrule(lr){7-8} \cmidrule(lr){10-11}
& & & Bursti- & GoBursti- & & Bursti- & GoBursti- & & Bursti- & GoBursti- \\
Parameter & True & & VAR & VAR & & VAR & VAR & & VAR & VAR \\"""

def two_model_block(title, plist, col, dp, show_true, bold=None):
    out = [f"\\multicolumn{{11}}{{l}}{{\\textit{{{title}}}}} \\\\"]
    for label, true, kb, kg in plist:
        cells = []
        for c in CONDS:
            for src, key in ((BI, kb), (GO, kg)):
                v = lead(src[c][key][col], dp)
                if bold and col == "Coverage" and (c, key, src is GO) in bold:
                    v = f"\\textbf{{{v}}}"
                cells.append(v)
        t = true if show_true else ""
        out.append(f"{label} & {t} & & {cells[0]} & {cells[1]} & & {cells[2]} & {cells[3]} & & {cells[4]} & {cells[5]} \\\\")
    return "\n".join(out)

def two_model_table(caption, tablabel, spec, note, bold=None):
    body = []
    for i, (title, plist, col, dp, show_true) in enumerate(spec):
        if i and i % 2 == 0:
            body.append("\\midrule")
        body.append(two_model_block(title, plist, col, dp, show_true, bold))
        if i % 2 == 0:
            body.append("\\addlinespace")
        elif i < len(spec) - 1:
            body.append("\\addlinespace")
    return (f"\\begin{{table}}[htbp]\n\\centering\n\\caption{{{caption}}}\n\\label{{{tablabel}}}\n"
            "\\begin{threeparttable}\n\\small\n\\setlength{\\tabcolsep}{3.5pt}\n"
            "\\begin{tabular}{llrrrrrrrrr}\n\\toprule\n" + HEAD2 + "\n\\midrule\n"
            + "\n".join(body) + "\n\\bottomrule\n\\end{tabular}\n"
            "\\begin{tablenotes}[flushleft]\n\\small\n\\item \\textit{Note.} " + note +
            "\n\\end{tablenotes}\n\\end{threeparttable}\n\\end{table}")

# ---------------- Table 4: study2_var (rBias + RMSE) ----------------
NOTE4 = (r"$B = 5$ bursts in all conditions. rBias = relative bias (bias/true value); RMSE = root mean "
         r"square error. BurstiVAR freely estimates burst intercepts; GoBurstiVAR constrains intercepts "
         r"to follow a Gompertz function (correctly specified for these data). 100 replications per "
         r"condition. The $T_b = 5$ conditions are the primary comparisons; $T_b = 20$, $N = 100$ is a "
         r"supplementary condition. Coverage and standard-error accuracy for the same parameters are "
         r"reported in Table~\ref{tab:study2_supp_var}. Power was 1.00 for all nonzero VAR fixed effects "
         r"in both models across all conditions except at $T_b = 5$, $N = 100$, where BurstiVAR showed "
         r"reduced power on both CR effects: $\mu_{\phi_{12}}$ (CR: $y_2 \to y_1$; true "
         r"value = $-0.10$) had power of .84 versus .98 under GoBurstiVAR, and $\mu_{\phi_{21}}$ (CR: "
         r"$y_1 \to y_2$; true value = $-0.15$) had power of .98 versus 1.00 under GoBurstiVAR. Power "
         r"for all VAR random effect SDs was 1.00 in all conditions and both models.")
TABLE4 = two_model_table(
    "Study 2: VAR and Innovation Covariance Parameter Recovery by BurstiVAR Versus GoBurstiVAR",
    "tab:study2_var",
    [("Fixed Effects: rBias",        VAR_FE, "rBias", 3, True),
     ("Fixed Effects: RMSE",         VAR_FE, "RMSE",  3, False),
     ("Random Effect SDs: rBias",    VAR_RE, "rBias", 3, True),
     ("Random Effect SDs: RMSE",     VAR_RE, "RMSE",  3, False),
     ("Innovation Covariance: rBias", VAR_IN, "rBias", 3, True),
     ("Innovation Covariance: RMSE",  VAR_IN, "RMSE",  3, False)],
    NOTE4)

# ---------------- supp2_var (Coverage + RDSE) ----------------
BOLD_SUPP_VAR = set()
for c in CONDS:
    for _, _, kb, kg in VAR_FE + VAR_RE + VAR_IN:
        if float(BI[c][kb]["Coverage"]) < 0.90: BOLD_SUPP_VAR.add((c, kb, False))
        if float(GO[c][kg]["Coverage"]) < 0.90: BOLD_SUPP_VAR.add((c, kg, True))
NOTE_S4 = (r"Companion to Table~\ref{tab:study2_var}. Cov = coverage of 95\% credible intervals "
           r"(nominal $.95$); RDSE = relative difference in SEs, "
           r"$(\overline{\mathrm{SE}}_{\text{post}} - \mathrm{Emp\ SD})/\mathrm{Emp\ SD}$, where "
           r"$\overline{\mathrm{SE}}_{\text{post}}$ is the mean posterior SD and Emp SD is the empirical "
           r"SD of the point estimates across replications; values near $0$ mean the posterior SD matches "
           r"the estimator's true variability, positive values mean it overstates it. Bold coverage values "
           r"fall below $.90$. Each cell summarizes 100 replications; $B = 5$ in all conditions.")
SUPP2VAR = two_model_table(
    "Study 2: Supplementary Recovery Metrics for the VAR and Innovation Covariance Parameters",
    "tab:study2_supp_var",
    [("Fixed Effects: Coverage",        VAR_FE, "Coverage", 2, True),
     ("Fixed Effects: RDSE",            VAR_FE, "RDSE",     2, False),
     ("Random Effect SDs: Coverage",    VAR_RE, "Coverage", 2, True),
     ("Random Effect SDs: RDSE",        VAR_RE, "RDSE",     2, False),
     ("Innovation Covariance: Coverage", VAR_IN, "Coverage", 2, True),
     ("Innovation Covariance: RDSE",     VAR_IN, "RDSE",     2, False)],
    NOTE_S4, bold=BOLD_SUPP_VAR)

# ---------------- trend params (BurstiVAR intercepts) ----------------
def trend_rows():
    rows = []
    for b in range(1, 6):
        for v, vname in ((1, "y_1"), (2, "y_2")):
            key = f"Level2Mean[{(b-1)*2 + v}]"
            rows.append((rf"$\mu_{{\alpha_{{{v}{b}}}}}$ (${vname}$, burst {b})", key))
    sds = []
    for b in range(1, 6):
        for v, vname in ((1, "y_1"), (2, "y_2")):
            key = f"Level2Sigma[{(b-1)*2 + v}]"
            sds.append((rf"$\sigma_{{\alpha_{{{v}{b}}}}}$ (${vname}$, burst {b})", key))
    return rows, sds
TR_MEANS, TR_SDS = trend_rows()

def one_model_block(ncols, title, plist, src, cols_dps, with_true):
    out = [f"\\multicolumn{{{ncols}}}{{l}}{{\\quad \\textit{{{title}}}}} \\\\"]
    for label, key in plist:
        cells = []
        for c in CONDS:
            cells.append(" & ".join(lead(src[c][key][col], dp) for col, dp in cols_dps))
        prefix = label + (f" & {lead(src[CONDS[0]][key]['TruePar'], 2)} &" if with_true else "")
        out.append(prefix + " & " + " & & ".join(cells) + r" \\")
    return "\n".join(out)

# Table 5: NEW trend table (True, Est, rBias, RMSE x3)
T5_BODY = []
for Nlab, plist in (("Burst-Intercept Means", TR_MEANS), ("Burst-Intercept Between-Person SDs", TR_SDS)):
    T5_BODY.append(one_model_block(14, Nlab, plist, TR, [("MeanPar_hat", 3), ("rBias", 3), ("RMSE", 3)], True))
    T5_BODY.append("\\addlinespace")
T5_BODY = "\n".join(T5_BODY[:-1])
NOTE5 = (r"$\mu_{\alpha_{vb}}$ and $\sigma_{\alpha_{vb}}$ denote the population mean and the between-person "
         r"SD of the burst-$b$ intercept of variable $v$ under BurstiVAR, fitted to data generated from "
         r"GoBurstiVAR. True values are the Gompertz-implied population means and between-person SDs of "
         r"the burst intercepts, derived by Monte Carlo integration ($2 \times 10^6$ draws) over the "
         r"person-level Gompertz parameter distribution used in data generation. Est Mean = mean posterior "
         r"estimate across replications; rBias = relative bias (bias/true value); RMSE = root mean square "
         r"error. Each cell summarizes 100 replications; $B = 5$ in all conditions. Coverage and "
         r"standard-error accuracy for the same parameters are reported in Table~\ref{tab:study2_supp_trend}.")
TABLE5 = (r"""\begin{table}[htbp]
\centering
\caption{Study 2: BurstiVAR's Free Burst Intercepts Against the Gompertz-Implied Trend}
\label{tab:study2_trend}
\begin{threeparttable}
\small
\setlength{\tabcolsep}{3pt}
\begin{tabular}{llrrrrrrrrrrrr}
\toprule
& & & \multicolumn{3}{c}{$T_b = 5$, $N = 100$} & & \multicolumn{3}{c}{$T_b = 5$, $N = 500$} & & \multicolumn{3}{c}{$T_b = 20$, $N = 100$} \\
\cmidrule(lr){4-6} \cmidrule(lr){8-10} \cmidrule(lr){12-14}
Parameter & True & & \shortstack{Est\\Mean} & rBias & RMSE & & \shortstack{Est\\Mean} & rBias & RMSE & & \shortstack{Est\\Mean} & rBias & RMSE \\
\midrule
""" + T5_BODY + "\n\\bottomrule\n\\end{tabular}\n\\begin{tablenotes}[flushleft]\n\\small\n"
    r"\item \textit{Note.} " + NOTE5 + "\n\\end{tablenotes}\n\\end{threeparttable}\n\\end{table}")

# supp2_trend: Coverage + RDSE x3
S5_BODY = []
for Nlab, plist in (("Burst-Intercept Means", TR_MEANS), ("Burst-Intercept Between-Person SDs", TR_SDS)):
    S5_BODY.append(one_model_block(9, Nlab, plist, TR, [("Coverage", 2), ("RDSE", 2)], False))
    S5_BODY.append("\\addlinespace")
S5_BODY = "\n".join(S5_BODY[:-1])
NOTE_S5 = (r"Companion to Table~\ref{tab:study2_trend}. Cov = coverage of 95\% credible intervals against "
           r"the Gompertz-implied true values (nominal $.95$); RDSE = relative difference in SEs, "
           r"$(\overline{\mathrm{SE}}_{\text{post}} - \mathrm{Emp\ SD})/\mathrm{Emp\ SD}$. Each cell "
           r"summarizes 100 replications; $B = 5$ in all conditions.")
SUPP2TREND = (r"""\begin{table}[htbp]
\centering
\caption{Study 2: Supplementary Recovery Metrics for BurstiVAR's Free Burst Intercepts}
\label{tab:study2_supp_trend}
\begin{threeparttable}
\small
\setlength{\tabcolsep}{4pt}
\begin{tabular}{lrrrrrrrr}
\toprule
& \multicolumn{2}{c}{$T_b = 5$, $N = 100$} & & \multicolumn{2}{c}{$T_b = 5$, $N = 500$} & & \multicolumn{2}{c}{$T_b = 20$, $N = 100$} \\
\cmidrule(lr){2-3} \cmidrule(lr){5-6} \cmidrule(lr){8-9}
Parameter & Cov & RDSE & & Cov & RDSE & & Cov & RDSE \\
\midrule
""" + S5_BODY
    + "\n\\bottomrule\n\\end{tabular}\n\\begin{tablenotes}[flushleft]\n\\small\n"
    r"\item \textit{Note.} " + NOTE_S5 + "\n\\end{tablenotes}\n\\end{threeparttable}\n\\end{table}")

# ---------------- Gompertz params (GoBurstiVAR) ----------------
GOMP_FE = [
 (r"$\mu_{\theta_1}$ (Asymptote, $y_1$)",    "Level2Mean[1]"),
 (r"$\mu_{\theta_2}$ (Displacement, $y_1$)", "Level2Mean[2]"),
 (r"$\mu_{\theta_3}$ (Growth rate, $y_1$)",  "Level2Mean[3]"),
 (r"$\mu_{\theta_1}$ (Asymptote, $y_2$)",    "Level2Mean[4]"),
 (r"$\mu_{\theta_2}$ (Displacement, $y_2$)", "Level2Mean[5]"),
 (r"$\mu_{\theta_3}$ (Growth rate, $y_2$)",  "Level2Mean[6]"),
]
GOMP_RE = [
 (r"$\sigma_{\theta_1}$ ($y_1$)", "Level2Sigma[1]"),
 (r"$\sigma_{\theta_2}$ ($y_1$)", "Level2Sigma[2]"),
 (r"$\sigma_{\theta_3}$ ($y_1$)", "Level2Sigma[3]"),
 (r"$\sigma_{\theta_1}$ ($y_2$)", "Level2Sigma[4]"),
 (r"$\sigma_{\theta_2}$ ($y_2$)", "Level2Sigma[5]"),
 (r"$\sigma_{\theta_3}$ ($y_2$)", "Level2Sigma[6]"),
]
TRUE_FMT = {"Level2Mean[1]":"10","Level2Mean[2]":"3","Level2Mean[3]":"1.0",
            "Level2Mean[4]":"8","Level2Mean[5]":"2.5","Level2Mean[6]":"0.8",
            "Level2Sigma[1]":"2.0","Level2Sigma[2]":"0.5","Level2Sigma[3]":"0.2",
            "Level2Sigma[4]":"1.5","Level2Sigma[5]":"0.4","Level2Sigma[6]":"0.15"}

def gomp_block(ncols, title, plist, cols_dps, with_true, bold_cov=False):
    out = [f"\\multicolumn{{{ncols}}}{{l}}{{\\quad \\textit{{{title}}}}} \\\\"]
    for label, key in plist:
        cells = []
        for c in CONDS:
            vals = []
            for col, dp in cols_dps:
                v = lead(GO[c][key][col], dp)
                if bold_cov and col == "Coverage" and float(GO[c][key][col]) < 0.90:
                    v = f"\\textbf{{{v}}}"
                vals.append(v)
            cells.append(" & ".join(vals))
        prefix = label + (f" & {TRUE_FMT[key]} &" if with_true else "")
        out.append(prefix + " & " + " & & ".join(cells) + r" \\")
    return "\n".join(out)

T6_BODY = (gomp_block(14, "Gompertz Fixed Effects", GOMP_FE, [("MeanPar_hat",3),("rBias",3),("RMSE",3)], True)
           + "\n\\addlinespace\n"
           + gomp_block(14, "Gompertz Random Effect SDs", GOMP_RE, [("MeanPar_hat",3),("rBias",3),("RMSE",3)], True))
NOTE6 = (r"GoBurstiVAR constrains burst intercepts to a person-specific Gompertz function of burst number "
         r"(Equations~\ref{eq:gompertz}--\ref{eq:gompertz_re}); data were generated from this correctly "
         r"specified model. $\mu_{\theta}$ and $\sigma_{\theta}$ denote the population mean and the "
         r"between-person SD of each Gompertz parameter. Est Mean = mean posterior estimate across "
         r"replications; rBias = relative bias (bias/true value); RMSE = root mean square error. "
         r"$B = 5$ bursts in all conditions; 100 replications per condition. Coverage and standard-error "
         r"accuracy for the same parameters are reported in Table~\ref{tab:study2_supp_gompertz}.")
TABLE6 = (r"""\begin{table}[htbp]
\centering
\caption{Study 2: Gompertz Parameter Recovery by GoBurstiVAR}
\label{tab:study2_gompertz}
\begin{threeparttable}
\small
\setlength{\tabcolsep}{3pt}
\begin{tabular}{llrrrrrrrrrrrr}
\toprule
& & & \multicolumn{3}{c}{$T_b = 5$, $N = 100$} & & \multicolumn{3}{c}{$T_b = 5$, $N = 500$} & & \multicolumn{3}{c}{$T_b = 20$, $N = 100$} \\
\cmidrule(lr){4-6} \cmidrule(lr){8-10} \cmidrule(lr){12-14}
Parameter & True & & \shortstack{Est\\Mean} & rBias & RMSE & & \shortstack{Est\\Mean} & rBias & RMSE & & \shortstack{Est\\Mean} & rBias & RMSE \\
\midrule
""" + T6_BODY + "\n\\bottomrule\n\\end{tabular}\n\\begin{tablenotes}[flushleft]\n\\small\n"
    r"\item \textit{Note.} " + NOTE6 + "\n\\end{tablenotes}\n\\end{threeparttable}\n\\end{table}")

S6_BODY = (gomp_block(9, "Gompertz Fixed Effects", GOMP_FE, [("Coverage",2),("RDSE",2)], False, bold_cov=True)
           + "\n\\addlinespace\n"
           + gomp_block(9, "Gompertz Random Effect SDs", GOMP_RE, [("Coverage",2),("RDSE",2)], False, bold_cov=True))
NOTE_S6 = (r"Companion to Table~\ref{tab:study2_gompertz}. Cov = coverage of 95\% credible intervals "
           r"(nominal $.95$); RDSE = relative difference in SEs, "
           r"$(\overline{\mathrm{SE}}_{\text{post}} - \mathrm{Emp\ SD})/\mathrm{Emp\ SD}$. Bold coverage "
           r"values fall below $.90$. Each cell summarizes 100 replications; $B = 5$ in all conditions.")
SUPP2GOMP = (r"""\begin{table}[htbp]
\centering
\caption{Study 2: Supplementary Recovery Metrics for the Gompertz Parameters}
\label{tab:study2_supp_gompertz}
\begin{threeparttable}
\small
\setlength{\tabcolsep}{4pt}
\begin{tabular}{lrrrrrrrr}
\toprule
& \multicolumn{2}{c}{$T_b = 5$, $N = 100$} & & \multicolumn{2}{c}{$T_b = 5$, $N = 500$} & & \multicolumn{2}{c}{$T_b = 20$, $N = 100$} \\
\cmidrule(lr){2-3} \cmidrule(lr){5-6} \cmidrule(lr){8-9}
Parameter & Cov & RDSE & & Cov & RDSE & & Cov & RDSE \\
\midrule
""" + S6_BODY + "\n\\bottomrule\n\\end{tabular}\n\\begin{tablenotes}[flushleft]\n\\small\n"
    r"\item \textit{Note.} " + NOTE_S6 + "\n\\end{tablenotes}\n\\end{threeparttable}\n\\end{table}")

for fname, block in [("table4_var.tex", TABLE4), ("table5_trend.tex", TABLE5),
                     ("table6_gompertz.tex", TABLE6), ("supp2_var.tex", SUPP2VAR),
                     ("supp2_trend.tex", SUPP2TREND), ("supp2_gompertz.tex", SUPP2GOMP)]:
    with open(os.path.join(HERE, fname), "w", encoding="utf-8") as fh:
        fh.write(block + "\n")
    print(f"wrote {fname}  ({block.count(chr(10)) + 1} lines)")
