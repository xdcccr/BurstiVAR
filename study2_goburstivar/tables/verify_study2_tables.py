#!/usr/bin/env python3
# --------------------------------------------------------------------
# BurstiVAR public repository -- Study 2 (GoBurstiVAR).
# Original file: paper_simulation_archive/Study2_Gompertz/verify_study2_tables.py
# Modification for this repository: (1) input CSV paths now point to
# ../results/MCfiles and ../results/TrendRecovery_B5; run this script
# from this directory, after gen_study2_tables.py. (2) MAIN now points
# to a main.tex beside this script (the manuscript is not distributed),
# and the open(MAIN) call was moved below the early exit so the disabled
# old-table regression legs (RUN_OLD_REGRESSION = False) no longer
# require the manuscript file to be present.
# All other content is identical to the version that produced the
# published results.
# --------------------------------------------------------------------

# Verify the six generated Study 2 table blocks cell-by-cell against the source
# CSVs, and regression-check against the OLD in-manuscript tables:
#   - old tab:study2_var rBias/Coverage cells == new table4/supp2_var cells
#   - old appendix tab:study2_gompertz rBias/Cov == new table6/supp2_gompertz
# Exits nonzero on any mismatch.
import csv, os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__))
MAIN = os.path.join(HERE, "main.tex")  # manuscript .tex (not distributed; needed only if RUN_OLD_REGRESSION is re-enabled)

CONDS = ["nT25_nP100", "nT25_nP500", "nT100_nP100"]
def load(path):
    d = {}
    with open(path, newline="") as fh:
        rd = csv.reader(fh); hdr = next(rd)
        for r in rd:
            if r and r[0]: d[r[0]] = dict(zip(hdr, r))
    return d
BI = {c: load(os.path.join(HERE, "..", "results", "MCfiles", f"mlGVARNoCorNoME_BurstIntercept_5Burst_on_mlGVARNoCorNoME_GompertzBurst_MCfileSumm_{c}.csv")) for c in CONDS}
GO = {c: load(os.path.join(HERE, "..", "results", "MCfiles", f"mlGVARNoCorNoME_GompertzBurst_TrueModel_MCfileSumm_{c}.csv")) for c in CONDS}
TR = {c: load(os.path.join(HERE, "..", "results", "TrendRecovery_B5", f"mlGVARNoCorNoME_BurstIntercept_5Burst_TrendRecovery_MCfileSumm_{c}.csv")) for c in CONDS}

def lead(v, n):
    x = float(v)
    s = f"{abs(x):.{n}f}"
    if float(s) < 1: s = "." + s.split(".")[1]
    return f"$-{s}$" if x < 0 else s
def unbold(s):
    m = re.fullmatch(r"\\textbf\{(.+)\}", s.strip())
    return m.group(1) if m else s.strip()

errs = 0
def chk(where, got, want):
    global errs
    if got != want:
        errs += 1; print(f"  MISMATCH {where}: table={got!r} expected={want!r}")

def datarows(text):
    """yield (label, [cells]) for rows 'label & c & c ... \\\\' skipping multicolumn/headers"""
    for ln in text.splitlines():
        ln = ln.strip()
        if not ln.endswith(r"\\") or "multicolumn" in ln or ln.startswith(("%", r"\cmidrule", r"\midrule")):
            continue
        parts = [p.strip() for p in ln[:-2].split("&")]
        if len(parts) >= 3 and ("$" in parts[0] or parts[0]):
            yield parts[0], parts[1:]

# maps: label -> (BiVAR key, GoVAR key) for two-model tables
LBL2 = {}
for lbl, key_b, key_g in [
    (r"$\mu_{\phi_{11}}$ (AR: $y_1$)", "Level2Mean[11]", "Level2Mean[7]"),
    (r"$\mu_{\phi_{22}}$ (AR: $y_2$)", "Level2Mean[12]", "Level2Mean[8]"),
    (r"$\mu_{\phi_{21}}$ (CR: $y_1 \to y_2$)", "Level2Mean[13]", "Level2Mean[9]"),
    (r"$\mu_{\phi_{12}}$ (CR: $y_2 \to y_1$)", "Level2Mean[14]", "Level2Mean[10]"),
    (r"$\sigma_{\phi_{11}}$ (AR: $y_1$)", "Level2Sigma[11]", "Level2Sigma[7]"),
    (r"$\sigma_{\phi_{22}}$ (AR: $y_2$)", "Level2Sigma[12]", "Level2Sigma[8]"),
    (r"$\sigma_{\phi_{21}}$ (CR: $y_1 \to y_2$)", "Level2Sigma[13]", "Level2Sigma[9]"),
    (r"$\sigma_{\phi_{12}}$ (CR: $y_2 \to y_1$)", "Level2Sigma[14]", "Level2Sigma[10]"),
    (r"$\sigma^2_{\varepsilon_1}$", "sigma_innovation[1,1]", "sigma_innovation[1,1]"),
    (r"$\sigma^2_{\varepsilon_2}$", "sigma_innovation[2,2]", "sigma_innovation[2,2]"),
    (r"$\sigma_{\varepsilon_{12}}$", "sigma_innovation[1,2]", "sigma_innovation[1,2]")]:
    LBL2[lbl] = (key_b, key_g)

def verify_two_model(fname, blocks_cols):
    """blocks_cols: list of (block_title_regex, csv_col, dp) in order of appearance"""
    text = open(os.path.join(HERE, fname), encoding="utf-8").read()
    # split into blocks by the multicolumn titles
    segs = re.split(r"\\multicolumn\{11\}\{l\}\{\\textit\{([^}]*)\}\} \\\\", text)
    n = 0
    for i in range(1, len(segs), 2):
        title, seg = segs[i], segs[i+1]
        col = dp = None
        for trex, c, d in blocks_cols:
            if re.fullmatch(trex, title): col, dp = c, d
        assert col, f"unmatched block title {title!r}"
        for lbl, cells in datarows(seg):
            kb, kg = LBL2[lbl]
            data = cells[1:]  # drop True col
            vals = [x for x in data if x != ""]
            assert len(vals) == 6, (lbl, cells)
            for j, c in enumerate(CONDS):
                chk(f"{fname} [{title}] {lbl} {c} BiVAR", unbold(vals[2*j]), lead(BI[c][kb][col], dp)); n += 1
                chk(f"{fname} [{title}] {lbl} {c} GoVAR", unbold(vals[2*j+1]), lead(GO[c][kg][col], dp)); n += 1
    print(f"{fname}: {n} cells checked")

verify_two_model("table4_var.tex", [(r".*rBias", "rBias", 3), (r".*RMSE", "RMSE", 3)])
verify_two_model("supp2_var.tex", [(r".*Coverage", "Coverage", 2), (r".*RDSE", "RDSE", 2)])

# one-model tables: trend (TR) and gompertz (GO)
def key_from_alpha(lbl):
    m = re.match(r"\$(\\mu|\\sigma)_\{\\alpha_\{(\d)(\d)\}\}\$", lbl)
    if not m: return None
    v, b = int(m.group(2)), int(m.group(3))
    idx = (b-1)*2 + v
    return ("Level2Mean" if m.group(1) == r"\mu" else "Level2Sigma") + f"[{idx}]"
GKEYS = {r"$\mu_{\theta_1}$ (Asymptote, $y_1$)": "Level2Mean[1]",
         r"$\mu_{\theta_2}$ (Displacement, $y_1$)": "Level2Mean[2]",
         r"$\mu_{\theta_3}$ (Growth rate, $y_1$)": "Level2Mean[3]",
         r"$\mu_{\theta_1}$ (Asymptote, $y_2$)": "Level2Mean[4]",
         r"$\mu_{\theta_2}$ (Displacement, $y_2$)": "Level2Mean[5]",
         r"$\mu_{\theta_3}$ (Growth rate, $y_2$)": "Level2Mean[6]",
         r"$\sigma_{\theta_1}$ ($y_1$)": "Level2Sigma[1]",
         r"$\sigma_{\theta_2}$ ($y_1$)": "Level2Sigma[2]",
         r"$\sigma_{\theta_3}$ ($y_1$)": "Level2Sigma[3]",
         r"$\sigma_{\theta_1}$ ($y_2$)": "Level2Sigma[4]",
         r"$\sigma_{\theta_2}$ ($y_2$)": "Level2Sigma[5]",
         r"$\sigma_{\theta_3}$ ($y_2$)": "Level2Sigma[6]"}

def verify_one_model(fname, src, keyfn, cols_dps, has_true, true_dp=2):
    text = open(os.path.join(HERE, fname), encoding="utf-8").read()
    n = 0
    for lbl, cells in datarows(text):
        key = keyfn(lbl)
        if key is None: continue
        vals = [x for x in cells if x != ""]
        if has_true:
            tv = vals[0]; vals = vals[1:]
            if fname.startswith("table5"):
                chk(f"{fname} {lbl} True", tv, lead(src[CONDS[0]][key]["TruePar"], true_dp)); n += 1
        per = len(cols_dps)
        assert len(vals) == 3*per, (lbl, cells)
        for j, c in enumerate(CONDS):
            for k, (col, dp) in enumerate(cols_dps):
                chk(f"{fname} {lbl} {c} {col}", unbold(vals[per*j+k]), lead(src[c][key][col], dp)); n += 1
    print(f"{fname}: {n} cells checked")

verify_one_model("table5_trend.tex", TR, key_from_alpha, [("MeanPar_hat",3),("rBias",3),("RMSE",3)], True)
verify_one_model("supp2_trend.tex", TR, key_from_alpha, [("Coverage",2),("RDSE",2)], False)
verify_one_model("table6_gompertz.tex", GO, lambda l: GKEYS.get(l), [("MeanPar_hat",3),("rBias",3),("RMSE",3)], True)
verify_one_model("supp2_gompertz.tex", GO, lambda l: GKEYS.get(l), [("Coverage",2),("RDSE",2)], False)

# ---- regression vs OLD manuscript tables ----
# OBSOLETE as of the 2026-07-16/17 installation: the old rBias+Cov tables were
# replaced/deleted in main.tex, so these legs would now re-parse the NEW blocks
# with the OLD layout. They passed against the pre-installation manuscript
# (203/204 identical; 1 known $-.000$ sign-display cell). Kept for the record.
RUN_OLD_REGRESSION = False
if not RUN_OLD_REGRESSION:
    print(f"\n(old-table regression legs skipped — tables already replaced)\nTOTAL mismatches = {errs}")
    sys.exit(1 if errs else 0)
main = open(MAIN, encoding="utf-8").read()
def grab(label):
    ci = main.index(rf"\label{{{label}}}")
    bi = main.rindex(r"\begin{table}", 0, ci)
    ei = main.index(r"\end{table}", ci)
    return main[bi:ei]
n3 = 0
old = grab("tab:study2_var")
segs = re.split(r"\\multicolumn\{11\}\{l\}\{\\textit\{([^}]*)\}\} \\\\", old)
for i in range(1, len(segs), 2):
    title, seg = segs[i], segs[i+1]
    col, dp = ("rBias", 3) if "rBias" in title else ("Coverage", 2)
    for lbl, cells in datarows(seg):
        if lbl not in LBL2: continue
        kb, kg = LBL2[lbl]
        vals = [x for x in cells[1:] if x != ""]
        for j, c in enumerate(CONDS):
            chk(f"OLD study2_var [{title}] {lbl} {c} BiVAR", unbold(vals[2*j]), lead(BI[c][kb][col], dp)); n3 += 1
            chk(f"OLD study2_var [{title}] {lbl} {c} GoVAR", unbold(vals[2*j+1]), lead(GO[c][kg][col], dp)); n3 += 1
print(f"OLD tab:study2_var regression: {n3} cells checked")

n4 = 0
oldg = grab("tab:study2_gompertz")
for lbl, cells in datarows(oldg):
    key = GKEYS.get(lbl)
    if key is None: continue
    vals = [x for x in cells if x != ""]  # True, then (rBias, Cov) x3
    assert len(vals) == 7, (lbl, cells)
    for j, c in enumerate(CONDS):
        chk(f"OLD gompertz {lbl} {c} rBias", unbold(vals[1+2*j]), lead(GO[c][key]["rBias"], 3)); n4 += 1
        chk(f"OLD gompertz {lbl} {c} Cov", unbold(vals[2+2*j]), lead(GO[c][key]["Coverage"], 2)); n4 += 1
print(f"OLD tab:study2_gompertz regression: {n4} cells checked")

print(f"\nTOTAL mismatches = {errs}")
sys.exit(1 if errs else 0)
