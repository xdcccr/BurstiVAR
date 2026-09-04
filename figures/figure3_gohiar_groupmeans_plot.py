# Provenance: MOL/figure3_gohiar_groupmeans_plot.py in the authors'
# project directory (D:/XXYDATAanalysis/IP_1b).
# Modified for this public repository: the two absolute output paths
# (OUT_PDF, which pointed into the Overleaf folder, and OUT_PNG, which
# pointed to a local scratch folder) were replaced with relative
# filenames so the script writes into the directory it is run from.
# The input path was not changed: it already resolves relative to this
# script, to BurstiVAR_results_GroupCovFull/
# GoHIAR_BurstiVAR_4Burst_GroupCovFull_GroupMeans.csv, an aggregate
# posterior summary (group-by-phase means and 95% CIs; no person-level
# data) written by GoHIAR_BurstiVAR_4Burst_GroupCovFull_Modelfit.R.
# Run it from this directory. All other content is identical to the
# version that produced the published results.

r"""
Figure 3 (manuscript: fig:gohiar_phase_means) — phase-level posterior means
by group for the mobile positive-psychology PREVENTION study application.

Data source (final analysis, CovFull specification):
    BurstiVAR_results_GroupCovFull/GoHIAR_BurstiVAR_4Burst_GroupCovFull_GroupMeans.csv
    (written by GoHIAR_BurstiVAR_4Burst_GroupCovFull_Modelfit.R; posterior means
    and 95% CIs of the reconstructed group-by-phase means, mu_ctrl / mu_intv.
    Verified against manuscript in-text values: P1 Meaning 75.64 [72.85, 78.45],
    P1 Relationship 75.17 [72.45, 77.89].)

Output:
    D:\XXYDATAanalysis\IP_1b\overleaf\figure3_gohiar_groupmeans.pdf

History:
    2026-07 — legend "Intervention" -> "PPI+Med" (matches the group label used
    in the manuscript and in Heshmati et al., 2025). After checking both source
    papers (Heshmati 2025; Li 2025 GoHiAR), the manuscript keeps INTERVENTION
    terminology throughout (matching the sources), so the x tick stays
    "P2: Intv". Plot style replicates the original PDF (matplotlib default
    DejaVu Sans; blue solid circles = Control, red dashed squares = PPI+Med;
    95% CI error bars; dotted reference line at 75).

Run:  py -3.12 figure3_gohiar_groupmeans_plot.py   (from MOL/)
"""
import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
CSV = HERE / 'BurstiVAR_results_GroupCovFull' / 'GoHIAR_BurstiVAR_4Burst_GroupCovFull_GroupMeans.csv'
OUT_PDF = Path('figure3_gohiar_groupmeans.pdf')
OUT_PNG = Path('fig3_new.png')

plt.rcParams.update({
    'font.size': 12,
    'axes.titlesize': 15,
    'axes.labelsize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.linewidth': 1.0,
})

# --- read the posterior group means ---
rows = list(csv.DictReader(open(CSV, encoding='utf-8')))
def series(variable, group):
    r = [x for x in rows if x['Variable'] == variable and x['Group'] == group]
    r.sort(key=lambda x: int(x['Burst']))
    mean = [float(x['Mean']) for x in r]
    lo = [float(x['LL95']) for x in r]
    hi = [float(x['UL95']) for x in r]
    return mean, lo, hi

GROUPS = [
    # (CSV group name, display label, color, linestyle, marker)
    ('Control',      'Control', '#1F77B4', '-',  'o'),
    ('Intervention', 'PPI+Med', '#D62728', '--', 's'),
]
PHASES = ['P1:\nPre', 'P2:\nIntv', 'P3:\nPost-I', 'P4:\nPost-II']
DODGE = {'Control': -0.07, 'Intervention': +0.07}

fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.4), sharey=True)

for ax, (var, panel) in zip(axes, [('Meaning', '(A) Meaning'),
                                   ('Relationship', '(B) Relationship')]):
    ax.axhline(75, color='0.75', linestyle=':', linewidth=1.0, zorder=1)
    for gname, glabel, color, ls, marker in GROUPS:
        mean, lo, hi = series(var, gname)
        x = [i + 1 + DODGE[gname] for i in range(4)]
        yerr = [[m - l for m, l in zip(mean, lo)],
                [h - m for m, h in zip(mean, hi)]]
        ax.errorbar(x, mean, yerr=yerr, color=color, linestyle=ls,
                    marker=marker, markersize=8, linewidth=2.2,
                    elinewidth=1.8, capsize=4, capthick=1.8,
                    label=glabel, zorder=3)
    ax.set_title(panel, fontweight='bold', loc='left')
    ax.set_xlim(0.6, 4.4)
    ax.set_xticks([1, 2, 3, 4])
    ax.set_xticklabels(PHASES)
    ax.set_ylim(66, 84)
    ax.set_yticks(range(66, 85, 2))

axes[0].set_ylabel('Posterior Mean (mPERMA, 0–100)')
axes[0].legend(loc='lower left', frameon=True, edgecolor='#cccccc',
               framealpha=1.0)

plt.tight_layout()
fig.savefig(OUT_PDF, facecolor='white', bbox_inches='tight')
fig.savefig(OUT_PNG, dpi=200, facecolor='white', bbox_inches='tight')
print('saved:', OUT_PDF)
