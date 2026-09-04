# Provenance: paper/figure1_v2_dotted_3person.py in the authors' project
# directory (D:/XXYDATAanalysis/IP_1b).
# Modified for this public repository: the two absolute output paths
# (OUT_PDF, which pointed into the Overleaf folder, and OUT_PNG, which
# pointed to a local scratch folder) were replaced with relative
# filenames so the script writes into the directory it is run from.
# Run it from this directory. All other content is identical to the
# version that produced the published results.

"""
Figure 1 v2 (2026-07): Panel A + three-person Panel B, replacing
figure1_horizontal_final_short.pdf in the Overleaf repo.

Panel A  — same data as the original short version (T=15, seed=42,
           sigma_eps=0.8) but with the M3-talk accessibility fix:
           black SOLID (beta=0.42, phi=0.2) vs red DOTTED (beta=0.40,
           phi=0.5), so the two series are not distinguished by color
           alone. Gray point halos kept from the talk version.
           (Source of the fix: _m3_fig/make_panelA_dotted.py)

Panel B  — composition follows the M3 talk's "BurstiVAR Idea" final
           frame (imgshorttalk/idea_p{1,2,3}_full.pdf): three per-person
           subpanels, each with a DISTINCT form-free trend type:
             Person 1: decelerating increase        (phi_i = 0.5)
             Person 2: damped oscillation           (phi_i = 0.3)
             Person 3: reversed-Gompertz decline    (phi_i = 0.6)
           Visual scheme kept from the ORIGINAL paper Panel B:
           person-specific colors (#E69F00 / #0072B2 / #009E73),
           data line + dots per burst, thinner person-colored dashed
           step lines for the burst intercepts alpha_ib, gray burst
           separators, numeric axes.
           (Sources: paper/figure1_horizontal.py styling +
            _m3_fig/make_idea_subfigs.py data/means/seeds)
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 11,
    'axes.linewidth': 1,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'legend.fontsize': 9,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'axes.spines.top': False,
    'axes.spines.right': False,
})

OUT_PDF = 'figure1_horizontal_final_short.pdf'
OUT_PNG = 'fig1_v2_preview.png'

# =============================================================================
# Panel A data — identical to the original generator (seed=42)
# =============================================================================
def simulate_ar1_with_trend(T, beta, phi, sigma_eps=1.0, seed=None):
    if seed is not None:
        np.random.seed(seed)
    t = np.arange(1, T + 1)
    trend = beta * t
    eps = np.random.normal(0, sigma_eps, T)
    eta = np.zeros(T)
    eta[0] = eps[0]
    for i in range(1, T):
        eta[i] = phi * eta[i - 1] + eps[i]
    return t, trend + eta

T_A = 15
t1, y1 = simulate_ar1_with_trend(T_A, 0.42, 0.2, 0.8, seed=42)   # black solid
t2, y2 = simulate_ar1_with_trend(T_A, 0.40, 0.5, 0.8, seed=42)   # red dotted

# =============================================================================
# Panel B data — identical to _m3_fig/make_idea_subfigs.py
# =============================================================================
B, Tb = 5, 10
people = [
    dict(means=[2.0, 6.0, 8.2, 9.2, 9.6], phi=0.50, seed=11, color='#E69F00',
         label=r'Person 1 ($\phi_i = 0.5$)'),
    dict(means=[8.5, 2.0, 6.5, 4.0, 5.2], phi=0.30, seed=22, color='#0072B2',
         label=r'Person 2 ($\phi_i = 0.3$)'),
    dict(means=[9.6, 9.0, 5.0, 1.3, 0.6], phi=0.60, seed=33, color='#009E73',
         label=r'Person 3 ($\phi_i = 0.6$)'),
]

def gen(means, phi, seed, sig=0.55):
    rng = np.random.RandomState(seed)
    t_all, y_all = [], []
    for b in range(B):
        m = means[b]
        t0 = b * Tb + 1
        yb = np.zeros(Tb)
        yb[0] = m + rng.normal(0, sig)
        for k in range(1, Tb):
            yb[k] = m + phi * (yb[k - 1] - m) + rng.normal(0, sig)
        t_all += list(np.arange(t0, t0 + Tb))
        y_all += list(yb)
    return np.array(t_all), np.array(y_all)

# =============================================================================
# Figure: one row — A (wider) + three B subpanels sharing the y axis
# =============================================================================
fig = plt.figure(figsize=(14, 4.1))
outer = fig.add_gridspec(1, 2, width_ratios=[1.55, 3.15],
                         wspace=0.14, left=0.05, right=0.995,
                         top=0.82, bottom=0.145)
gsB = outer[0, 1].subgridspec(1, 3, wspace=0.0)

# ----- Panel A -----
axA = fig.add_subplot(outer[0, 0])
axA.scatter(t1, y1, color='gray', s=42, alpha=0.35, edgecolors='none', zorder=1)
axA.scatter(t2, y2, color='gray', s=42, alpha=0.35, edgecolors='none', zorder=1)
axA.plot(t1, y1, color='#1a1a1a', linewidth=2.2, linestyle='-',
         label=r'$\beta = 0.42$, $\phi = 0.2$', zorder=3)
axA.plot(t2, y2, color='#D62728', linewidth=2.2, linestyle=':',
         label=r'$\beta = 0.40$, $\phi = 0.5$', zorder=2)
axA.set_xlabel('Time')
axA.set_ylabel(r'$y_t$')
axA.set_xlim(0, 16)
axA.legend(loc='upper left', frameon=True, fancybox=False,
           edgecolor='#888888', framealpha=0.95)

# ----- Panel B: three per-person subpanels, directly adjoining -----
axesB = []
for j, p in enumerate(people):
    ax = fig.add_subplot(gsB[0, j], sharey=axesB[0] if axesB else None)
    axesB.append(ax)
    t_p, y_p = gen(p['means'], p['phi'], p['seed'])
    for b in range(B):
        s_, e_ = b * Tb, (b + 1) * Tb
        # data line + scatter, per burst (paper Panel B scheme)
        ax.plot(t_p[s_:e_], y_p[s_:e_], color=p['color'], linewidth=1.2,
                alpha=0.7, zorder=3)
        ax.scatter(t_p[s_:e_], y_p[s_:e_], color=p['color'], s=15, alpha=0.5,
                   edgecolors='none', zorder=4)
        # person-level alpha_ib — thinner dashed step line, person color
        ax.hlines(p['means'][b], b * Tb + 1, (b + 1) * Tb,
                  colors=p['color'], linewidth=1.0, linestyle='--',
                  alpha=0.85, zorder=5)
    for sep in range(1, B):
        ax.axvline(sep * Tb + 0.5, color='gray', linestyle='-',
                   linewidth=0.5, alpha=0.3)
    ax.set_xlim(0, B * Tb + 1)
    # interior ticks only: endpoint labels would collide at the panel seams
    ax.set_xticks([10, 20, 30, 40])
    ax.set_title(p['label'], fontsize=11, fontweight='normal', pad=4)

axesB[1].set_xlabel('Time')          # one shared label under the middle panel
axesB[0].set_ylim(-0.4, 10.8)
axesB[0].set_yticks([0, 5, 10])
axesB[0].set_ylabel(r'$y_t$')
for ax in axesB[1:]:
    ax.tick_params(labelleft=False, left=False)

# ----- Panel headers, aligned on one top line -----
posA = axA.get_position()
posB1 = axesB[0].get_position()
fig.text(posA.x0, 0.92, '(A) Parameter Entanglement',
         fontsize=13, fontweight='bold', ha='left', va='center')
fig.text(posB1.x0, 0.92, '(B) Multilevel BurstiVAR',
         fontsize=13, fontweight='bold', ha='left', va='center')

fig.savefig(OUT_PDF, facecolor='white')
fig.savefig(OUT_PNG, dpi=170, facecolor='white')
print('saved:', OUT_PDF)
