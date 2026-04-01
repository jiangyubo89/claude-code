"""
Regenerate fig1_boxplot_comparison.png with:
  - Updated group labels (GenX / PFOA)
  - CLD letters removed
  - Footer text removed
  - Significance brackets vs Control
"""

import os, sys, io, re, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
import itertools

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
warnings.filterwarnings("ignore")

# ── Paths ────────────────────────────────────────────────────────────
BASE    = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = BASE          # save outputs to the same folder as this script
CSV     = os.path.join(BASE, "Results.csv")

# ── Style ────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 11,
    "axes.titlesize": 12,
    "axes.labelsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.8,
    "xtick.direction": "out",
    "ytick.direction": "out",
})

GROUP_ORDER = ["C", "G-0.5", "G-5", "P-0.5", "P-5"]

# New x-axis labels  (unit kept as mg/L)
X_LABELS = [
    "Control",
    "GenX\n0.5 mg/L",
    "GenX\n5.0 mg/L",
    "PFOA\n0.5 mg/L",
    "PFOA\n5.0 mg/L",
]

PALETTE = {
    "C":     "#4C72B0",
    "G-0.5": "#55A868",
    "G-5":   "#1B7837",
    "P-0.5": "#DD8452",
    "P-5":   "#B22222",
}

# ── Load & clean data (same pipeline as before) ──────────────────────
raw = pd.read_csv(CSV, header=0)
raw.columns = ["Index", "Label", "Length_cm"]
raw = raw.dropna(subset=["Length_cm"])
raw["Length_cm"] = pd.to_numeric(raw["Length_cm"], errors="coerce")
raw = raw.dropna(subset=["Length_cm"])

def parse_label(label):
    label = str(label).replace(".jpg", "").strip()
    m = re.match(r"^(C|G-[\d.]+|P-[\d.]+)-(\d+)$", label)
    return (m.group(1), int(m.group(2))) if m else (label, None)

raw[["Group", "Replicate"]] = raw["Label"].apply(
    lambda x: pd.Series(parse_label(x)))
raw = raw[raw["Group"].isin(GROUP_ORDER)].copy()

# IQR outlier removal
clean = []
for grp in GROUP_ORDER:
    sub = raw[raw["Group"] == grp].copy()
    Q1, Q3 = sub["Length_cm"].quantile([0.25, 0.75])
    IQR = Q3 - Q1
    lo, hi = Q1 - 1.5*IQR, Q3 + 1.5*IQR
    clean.append(sub[sub["Length_cm"].between(lo, hi)])
df = pd.concat(clean, ignore_index=True)

# Manual removal: lowest point in G-5 (0.247 cm, G-5-2.jpg) — user request
before = len(df)
df = df[~((df["Group"] == "G-5") & (df["Length_cm"] == 0.247))].copy()
print(f"Manual removal of G-5 point 0.247 cm: {before} -> {len(df)} rows")

group_vals = [df[df["Group"] == g]["Length_cm"].values for g in GROUP_ORDER]

print("Final n per group:")
for g in GROUP_ORDER:
    print(f"  {g}: n={len(df[df['Group']==g])}")

# ── Post-hoc: Control vs each treatment (Mann-Whitney + Bonferroni) ──
treatment_groups = ["G-0.5", "G-5", "P-0.5", "P-5"]
ctrl_vals = df[df["Group"] == "C"]["Length_cm"].values

raw_p, u_stats = [], []
for tg in treatment_groups:
    tv = df[df["Group"] == tg]["Length_cm"].values
    u, p = stats.mannwhitneyu(ctrl_vals, tv, alternative="two-sided")
    raw_p.append(p)
    u_stats.append(u)

_, p_adj, _, _ = multipletests(raw_p, alpha=0.05, method="bonferroni")

def sig_label(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return "ns"

sig_map = {tg: sig_label(pa) for tg, pa in zip(treatment_groups, p_adj)}
print("Significance vs Control:")
for tg, pa in zip(treatment_groups, p_adj):
    print(f"  {tg}: {sig_label(pa)}  (p_adj={pa:.4f})")

# ── Draw figure ──────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(9, 6))

positions = np.arange(len(GROUP_ORDER))

bp = ax.boxplot(
    group_vals,
    positions=positions,
    widths=0.45,
    patch_artist=True,
    whis=[0, 100],           # whiskers to actual data min/max (no double-filtering)
    medianprops=dict(color="black", linewidth=1.8),
    whiskerprops=dict(linewidth=0.9),
    capprops=dict(linewidth=0.9),
    flierprops=dict(marker="", linestyle="none"),
    showfliers=False,
)
for patch, grp in zip(bp["boxes"], GROUP_ORDER):
    patch.set_facecolor(PALETTE[grp])
    patch.set_alpha(0.72)
    patch.set_linewidth(0.8)

# Jittered individual points
np.random.seed(42)
for i, (grp, vals) in enumerate(zip(GROUP_ORDER, group_vals)):
    jitter = np.random.uniform(-0.18, 0.18, size=len(vals))
    ax.scatter(positions[i] + jitter, vals,
               color=PALETTE[grp], s=18, alpha=0.55,
               zorder=3, linewidths=0)

# ── Significance brackets: Control (pos 0) → each treatment ──────────
# Use actual data max (matches whis=[0,100])
tops = [v.max() for v in group_vals]
y_data_max = max(tops)

# Stagger bracket heights to avoid overlap
bracket_starts = [y_data_max + 0.010,
                  y_data_max + 0.022,
                  y_data_max + 0.034,
                  y_data_max + 0.046]

ctrl_pos = positions[0]

for k, tg in enumerate(treatment_groups):
    tg_idx = GROUP_ORDER.index(tg)
    tg_pos = positions[tg_idx]
    label  = sig_map[tg]
    y_br   = bracket_starts[k]
    tick   = 0.003          # vertical tick height

    # Horizontal bar + two vertical ticks
    ax.plot([ctrl_pos, tg_pos], [y_br, y_br],
            color="black", linewidth=0.9, clip_on=False)
    ax.plot([ctrl_pos, ctrl_pos], [y_br - tick, y_br],
            color="black", linewidth=0.9, clip_on=False)
    ax.plot([tg_pos, tg_pos], [y_br - tick, y_br],
            color="black", linewidth=0.9, clip_on=False)

    # Significance label centred on the bracket
    ax.text((ctrl_pos + tg_pos) / 2, y_br + 0.001,
            label, ha="center", va="bottom",
            fontsize=10 if label != "ns" else 9,
            fontweight="bold" if label != "ns" else "normal")

# ── Axes labels & title ──────────────────────────────────────────────
ax.set_xticks(positions)
ax.set_xticklabels(X_LABELS, fontsize=10)
ax.set_ylabel("Body Length (cm)", fontsize=12)
ax.set_xlabel("Treatment Group", fontsize=12)
ax.set_title("Body Length Comparison Across Treatment Groups",
             fontsize=12, pad=10)

# Extend y-axis to fit brackets
y_top_all = bracket_starts[-1] + 0.014
ax.set_ylim(ax.get_ylim()[0], y_top_all)

# Kruskal-Wallis annotation (top-right) — recomputed from current data
h_stat, kw_p = stats.kruskal(*group_vals)
n_total = sum(len(v) for v in group_vals)
k = len(group_vals)
eps2 = (h_stat - k + 1) / (n_total - k)
p_str = "p < 0.001" if kw_p < 0.001 else f"p = {kw_p:.4f}"
fig.text(0.98, 0.97,
         f"Kruskal-Wallis: H = {h_stat:.2f}, {p_str}\n"
         r"$\varepsilon^2$" + f" = {eps2:.3f}",
         ha="right", va="top",
         fontsize=9, style="italic",
         bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.75", alpha=0.85))

plt.tight_layout()
out_path = os.path.join(OUT_DIR, "fig1_boxplot_comparison.png")
fig.savefig(out_path)
plt.close(fig)
print(f"\nSaved: {out_path}")
