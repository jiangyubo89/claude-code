"""
PFOA Body Length Analysis
=========================
Publication-quality statistical analysis of zebrafish body length data
under different PFOA treatment groups and concentrations.

Groups:
  C      → Control (0 mg/L)
  G-0.5  → Treatment G, 0.5 mg/L
  G-5    → Treatment G, 5.0 mg/L
  P-0.5  → Treatment P, 0.5 mg/L
  P-5    → Treatment P, 5.0 mg/L
"""

import os
import re
import sys
import io
import warnings
import numpy as np

# Force UTF-8 output on Windows
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import scipy.stats as stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multitest import multipletests
import itertools

warnings.filterwarnings("ignore")

# ─────────────────────────────────────────────
# 0. Output directory
# ─────────────────────────────────────────────
OUT_DIR = os.path.join(os.path.dirname(__file__), "analysis_results")
os.makedirs(OUT_DIR, exist_ok=True)

# ─────────────────────────────────────────────
# Publication style
# ─────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 11,
    "axes.titlesize": 12,
    "axes.labelsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.8,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
})

# Colour palette (colour-blind friendly)
PALETTE = {
    "C":     "#4C72B0",   # blue
    "G-0.5": "#55A868",   # green
    "G-5":   "#1B7837",   # dark green
    "P-0.5": "#DD8452",   # orange
    "P-5":   "#B22222",   # dark red
}

GROUP_ORDER   = ["C", "G-0.5", "G-5", "P-0.5", "P-5"]
GROUP_LABELS  = {
    "C":     "Control",
    "G-0.5": "G 0.5 mg/L",
    "G-5":   "G 5.0 mg/L",
    "P-0.5": "P 0.5 mg/L",
    "P-5":   "P 5.0 mg/L",
}

# ─────────────────────────────────────────────
# 1. Data loading & cleaning
# ─────────────────────────────────────────────
print("=" * 60)
print("PFOA Body Length Statistical Analysis")
print("=" * 60)

csv_path = os.path.join(os.path.dirname(__file__), "Results.csv")
raw = pd.read_csv(csv_path, header=0)
raw.columns = ["Index", "Label", "Length_cm"]
raw = raw.dropna(subset=["Length_cm"])
raw["Length_cm"] = pd.to_numeric(raw["Length_cm"], errors="coerce")
raw = raw.dropna(subset=["Length_cm"])

def parse_label(label):
    """Extract group and replicate from image filename."""
    label = str(label).replace(".jpg", "").strip()
    m = re.match(r"^(C|G-[\d.]+|P-[\d.]+)-(\d+)$", label)
    if m:
        return m.group(1), int(m.group(2))
    return label, None

raw[["Group", "Replicate"]] = raw["Label"].apply(
    lambda x: pd.Series(parse_label(x))
)
raw = raw[raw["Group"].isin(GROUP_ORDER)].copy()
raw["Length_mm"] = raw["Length_cm"] * 10  # convert to mm

print(f"\nTotal raw measurements: {len(raw)}")
print(raw.groupby("Group")["Length_cm"].count().rename("n"))

# ─────────────────────────────────────────────
# 2. Outlier detection (IQR method, per group)
# ─────────────────────────────────────────────
outlier_log = []
clean_dfs = []

for grp in GROUP_ORDER:
    sub = raw[raw["Group"] == grp].copy()
    Q1 = sub["Length_cm"].quantile(0.25)
    Q3 = sub["Length_cm"].quantile(0.75)
    IQR = Q3 - Q1
    lo, hi = Q1 - 1.5 * IQR, Q3 + 1.5 * IQR
    out = sub[(sub["Length_cm"] < lo) | (sub["Length_cm"] > hi)]
    if len(out):
        for _, row in out.iterrows():
            outlier_log.append({
                "Group": grp, "Index": row["Index"],
                "Label": row["Label"], "Length_cm": row["Length_cm"],
                "IQR_lower": lo, "IQR_upper": hi
            })
    clean_dfs.append(sub[(sub["Length_cm"] >= lo) & (sub["Length_cm"] <= hi)])

df = pd.concat(clean_dfs, ignore_index=True)

print(f"\nOutliers removed: {len(outlier_log)}")
if outlier_log:
    for o in outlier_log:
        print(f"  [{o['Group']}] Label={o['Label']}, Length={o['Length_cm']:.3f} cm")

print(f"Clean measurements: {len(df)}")

# Save outlier log
if outlier_log:
    pd.DataFrame(outlier_log).to_csv(
        os.path.join(OUT_DIR, "outliers_removed.csv"), index=False)

# ─────────────────────────────────────────────
# 3. Descriptive statistics
# ─────────────────────────────────────────────
def desc_stats(series):
    n = len(series)
    mean = series.mean()
    sd   = series.std(ddof=1)
    se   = sd / np.sqrt(n)
    med  = series.median()
    iqr  = series.quantile(0.75) - series.quantile(0.25)
    cv   = (sd / mean) * 100
    ci95 = 1.96 * se
    return pd.Series({
        "n": n, "Mean": mean, "SD": sd, "SE": se,
        "Median": med, "IQR": iqr,
        "Min": series.min(), "Max": series.max(),
        "CV%": cv, "CI95±": ci95
    })

# Per group
desc_grp = df.groupby("Group")["Length_cm"].apply(desc_stats).unstack()
desc_grp = desc_grp.loc[GROUP_ORDER]
desc_grp.insert(0, "Group_label", [GROUP_LABELS[g] for g in GROUP_ORDER])

# Per replicate (image)
desc_rep = df.groupby(["Group", "Replicate"])["Length_cm"].apply(desc_stats).unstack()
desc_rep.to_csv(os.path.join(OUT_DIR, "descriptive_stats_per_replicate.csv"),
                float_format="%.4f")
desc_grp.to_csv(os.path.join(OUT_DIR, "descriptive_stats_per_group.csv"),
                float_format="%.4f")

print("\n--- Descriptive Statistics (per group, Length/cm) ---")
print(desc_grp[["Group_label", "n", "Mean", "SD", "SE", "Median", "Min", "Max", "CV%"]]
      .to_string(index=True))

# ─────────────────────────────────────────────
# 4. Normality test (Shapiro-Wilk)
# ─────────────────────────────────────────────
sw_results = []
for grp in GROUP_ORDER:
    vals = df[df["Group"] == grp]["Length_cm"].values
    stat, p = stats.shapiro(vals)
    sw_results.append({
        "Group": grp, "Group_label": GROUP_LABELS[grp],
        "n": len(vals), "W_statistic": stat, "p_value": p,
        "Normal (p>0.05)": "Yes" if p > 0.05 else "No"
    })

sw_df = pd.DataFrame(sw_results)
sw_df.to_csv(os.path.join(OUT_DIR, "normality_shapiro_wilk.csv"),
             index=False, float_format="%.4f")

print("\n--- Shapiro-Wilk Normality Test ---")
for r in sw_results:
    flag = "[OK]" if r["Normal (p>0.05)"] == "Yes" else "[!]"
    print(f"  {r['Group']:<8} W={r['W_statistic']:.4f}  p={r['p_value']:.4f}  {flag}")

all_normal = all(r["Normal (p>0.05)"] == "Yes" for r in sw_results)

# ─────────────────────────────────────────────
# 5. Homogeneity of variance (Levene's test)
# ─────────────────────────────────────────────
group_vals = [df[df["Group"] == g]["Length_cm"].values for g in GROUP_ORDER]
lev_stat, lev_p = stats.levene(*group_vals)
lev_result = {"Statistic": lev_stat, "p_value": lev_p,
              "Equal_variances (p>0.05)": "Yes" if lev_p > 0.05 else "No"}
pd.DataFrame([lev_result]).to_csv(
    os.path.join(OUT_DIR, "levene_test.csv"), index=False, float_format="%.4f")

print(f"\n--- Levene's Test for Homogeneity of Variance ---")
print(f"  F={lev_stat:.4f}  p={lev_p:.4f}  "
      f"Equal variances: {lev_result['Equal_variances (p>0.05)']}")

equal_var = lev_p > 0.05

# ─────────────────────────────────────────────
# 6. Main statistical test
# ─────────────────────────────────────────────
use_parametric = all_normal and equal_var
print(f"\n--- Test Selection ---")
print(f"  All groups normal: {all_normal}")
print(f"  Variances equal:   {equal_var}")

if use_parametric:
    print("  → One-way ANOVA selected")
    f_stat, anova_p = stats.f_oneway(*group_vals)
    main_test = "One-way ANOVA"
    main_stat_name = "F"
    main_stat = f_stat
    main_p = anova_p
    # Effect size: eta-squared
    grand = df["Length_cm"].values
    grand_mean = grand.mean()
    ss_between = sum(len(v) * (v.mean() - grand_mean)**2 for v in group_vals)
    ss_total   = sum((grand - grand_mean)**2)
    eta2 = ss_between / ss_total
    effect_name, effect_val = "η²", eta2
else:
    print("  → Kruskal-Wallis test selected (non-parametric)")
    h_stat, kw_p = stats.kruskal(*group_vals)
    main_test = "Kruskal-Wallis"
    main_stat_name = "H"
    main_stat = h_stat
    main_p = kw_p
    # Effect size: epsilon-squared
    n_total = sum(len(v) for v in group_vals)
    k = len(group_vals)
    eps2 = (h_stat - k + 1) / (n_total - k)
    effect_name, effect_val = "ε²", eps2

print(f"\n--- {main_test} ---")
print(f"  {main_stat_name}={main_stat:.4f}  p={main_p:.6f}  "
      f"{effect_name}={effect_val:.4f}")
sig = "Yes" if main_p < 0.05 else "No"
print(f"  Significant overall difference: {sig}")

main_result = {
    "Test": main_test, "Statistic_name": main_stat_name,
    "Statistic": main_stat, "p_value": main_p,
    f"Effect_size ({effect_name})": effect_val,
    "Significant (p<0.05)": sig
}
pd.DataFrame([main_result]).to_csv(
    os.path.join(OUT_DIR, "main_test_result.csv"), index=False, float_format="%.6f")

# ─────────────────────────────────────────────
# 7. Post-hoc tests
# ─────────────────────────────────────────────
pairs = list(itertools.combinations(GROUP_ORDER, 2))
posthoc_rows = []

if use_parametric:
    # Tukey HSD via statsmodels
    tukey = pairwise_tukeyhsd(
        df["Length_cm"].values,
        df["Group"].values,
        alpha=0.05
    )
    tukey_df = pd.DataFrame(
        data=tukey._results_table.data[1:],
        columns=tukey._results_table.data[0]
    )
    tukey_df.columns = ["Group1", "Group2", "Mean_diff", "p_adj",
                        "CI_lower", "CI_upper", "Reject"]
    tukey_df.to_csv(os.path.join(OUT_DIR, "posthoc_tukey_hsd.csv"),
                    index=False, float_format="%.6f")

    for _, row in tukey_df.iterrows():
        posthoc_rows.append({
            "Group1": row["Group1"], "Group2": row["Group2"],
            "Mean_diff": row["Mean_diff"], "p_adj": row["p_adj"],
            "Significant": "Yes" if row["Reject"] else "No"
        })
    posthoc_method = "Tukey HSD"
else:
    # Pairwise Mann-Whitney U + Bonferroni correction
    raw_p = []
    u_stats = []
    for g1, g2 in pairs:
        v1 = df[df["Group"] == g1]["Length_cm"].values
        v2 = df[df["Group"] == g2]["Length_cm"].values
        u, p = stats.mannwhitneyu(v1, v2, alternative="two-sided")
        raw_p.append(p)
        u_stats.append(u)

    reject, p_adj, _, _ = multipletests(raw_p, alpha=0.05, method="bonferroni")

    for (g1, g2), u, p_raw, p_a, rej in zip(pairs, u_stats, raw_p, p_adj, reject):
        v1 = df[df["Group"] == g1]["Length_cm"].values
        v2 = df[df["Group"] == g2]["Length_cm"].values
        mean_diff = v1.mean() - v2.mean()
        posthoc_rows.append({
            "Group1": g1, "Group2": g2,
            "U_statistic": u, "p_raw": p_raw,
            "p_adj_bonferroni": p_a,
            "Mean_diff": mean_diff,
            "Significant": "Yes" if rej else "No"
        })
    pd.DataFrame(posthoc_rows).to_csv(
        os.path.join(OUT_DIR, "posthoc_mannwhitney_bonferroni.csv"),
        index=False, float_format="%.6f")
    posthoc_method = "Mann-Whitney U + Bonferroni"

print(f"\n--- Post-hoc: {posthoc_method} ---")
for r in posthoc_rows:
    flag = "***" if float(r.get("p_adj", r.get("p_adj_bonferroni", 1))) < 0.001 else \
           "**"  if float(r.get("p_adj", r.get("p_adj_bonferroni", 1))) < 0.01  else \
           "*"   if float(r.get("p_adj", r.get("p_adj_bonferroni", 1))) < 0.05  else "ns"
    p_val = r.get("p_adj", r.get("p_adj_bonferroni", 1))
    print(f"  {r['Group1']:<8} vs {r['Group2']:<8}  "
          f"Δmean={r['Mean_diff']:+.4f}  p_adj={p_val:.4f}  {flag}")

# ─────────────────────────────────────────────
# Compact Letter Display (CLD) for annotations
# ─────────────────────────────────────────────
def compact_letters(groups, posthoc_rows, p_col="p_adj"):
    """Simple CLD assignment."""
    sig_pairs = set()
    for r in posthoc_rows:
        p = r.get(p_col, r.get("p_adj_bonferroni", 1))
        if float(p) < 0.05:
            sig_pairs.add((r["Group1"], r["Group2"]))
            sig_pairs.add((r["Group2"], r["Group1"]))

    letters = {g: set() for g in groups}
    letter_id = ord('a')

    assigned = {g: [] for g in groups}
    # Greedy CLD
    current_letter = 'a'
    used_letters = 0
    groups_sorted = groups[:]

    # Build sharing matrix
    def shares_letter(g1, g2):
        return (g1, g2) not in sig_pairs

    letter_sets = []  # list of sets of groups sharing a letter
    for g in groups_sorted:
        placed = False
        for ls in letter_sets:
            if all(shares_letter(g, other) for other in ls):
                ls.add(g)
                placed = True
                break
        if not placed:
            letter_sets.append({g})

    cld = {g: "" for g in groups_sorted}
    for i, ls in enumerate(letter_sets):
        for g in ls:
            cld[g] += chr(ord('a') + i)
    return cld

p_col = "p_adj" if use_parametric else "p_adj_bonferroni"
cld = compact_letters(GROUP_ORDER, posthoc_rows, p_col=p_col)
print("\n--- Compact Letter Display ---")
for g in GROUP_ORDER:
    print(f"  {g}: {cld[g]}")

# ─────────────────────────────────────────────
# 8. Cohen's d for pairwise effect sizes
# ─────────────────────────────────────────────
def cohens_d(v1, v2):
    pooled_sd = np.sqrt(((len(v1)-1)*v1.std(ddof=1)**2 +
                         (len(v2)-1)*v2.std(ddof=1)**2) /
                        (len(v1)+len(v2)-2))
    return (v1.mean() - v2.mean()) / pooled_sd if pooled_sd > 0 else 0

cohen_rows = []
for g1, g2 in pairs:
    v1 = df[df["Group"] == g1]["Length_cm"].values
    v2 = df[df["Group"] == g2]["Length_cm"].values
    d = cohens_d(v1, v2)
    mag = "negligible" if abs(d)<0.2 else "small" if abs(d)<0.5 else \
          "medium" if abs(d)<0.8 else "large"
    cohen_rows.append({"Group1": g1, "Group2": g2,
                       "Cohens_d": d, "Magnitude": mag})

pd.DataFrame(cohen_rows).to_csv(
    os.path.join(OUT_DIR, "effect_sizes_cohens_d.csv"),
    index=False, float_format="%.4f")

# ─────────────────────────────────────────────
# 9. Figures
# ─────────────────────────────────────────────

# ── Helper: significance bracket annotation ──────────────────────────
def sig_label(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return "ns"

def add_cld_labels(ax, cld, x_positions, y_max, y_step=0.012):
    """Add CLD letters above each bar/box."""
    for i, g in enumerate(GROUP_ORDER):
        ax.text(x_positions[i], y_max + y_step,
                cld[g], ha="center", va="bottom",
                fontsize=10, fontweight="bold")

# ── Figure 1: Box + strip plot ────────────────────────────────────────
fig1, ax1 = plt.subplots(figsize=(8, 5.5))

positions = np.arange(len(GROUP_ORDER))
bp = ax1.boxplot(
    group_vals,
    positions=positions,
    widths=0.45,
    patch_artist=True,
    medianprops=dict(color="black", linewidth=1.5),
    whiskerprops=dict(linewidth=0.8),
    capprops=dict(linewidth=0.8),
    flierprops=dict(marker="", linestyle="none"),
    showfliers=False
)

for patch, grp in zip(bp["boxes"], GROUP_ORDER):
    patch.set_facecolor(PALETTE[grp])
    patch.set_alpha(0.7)
    patch.set_linewidth(0.8)

# Jitter individual points
np.random.seed(42)
for i, (grp, vals) in enumerate(zip(GROUP_ORDER, group_vals)):
    jitter = np.random.uniform(-0.18, 0.18, size=len(vals))
    ax1.scatter(positions[i] + jitter, vals,
                color=PALETTE[grp], s=18, alpha=0.55,
                zorder=3, linewidths=0.3, edgecolors="none")

# CLD letters
y_top = max(v.max() for v in group_vals)
add_cld_labels(ax1, cld, positions, y_top, y_step=0.008)

ax1.set_xticks(positions)
ax1.set_xticklabels([GROUP_LABELS[g] for g in GROUP_ORDER], fontsize=10)
ax1.set_ylabel("Body Length (cm)", fontsize=12)
ax1.set_xlabel("Treatment Group", fontsize=12)
ax1.set_title("Body Length Comparison Across Treatment Groups", fontsize=12, pad=10)

# Legend for test result
p_str = f"p < 0.001" if main_p < 0.001 else f"p = {main_p:.4f}"
ax1.text(0.98, 0.97,
         f"{main_test}: {main_stat_name} = {main_stat:.2f}, {p_str}\n"
         f"{effect_name} = {effect_val:.3f}",
         transform=ax1.transAxes, ha="right", va="top",
         fontsize=9, style="italic",
         bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.8", alpha=0.8))

ax1.text(0.98, 0.02, "Different letters indicate significant differences (p < 0.05)",
         transform=ax1.transAxes, ha="right", va="bottom", fontsize=8, color="gray")

plt.tight_layout()
fig1.savefig(os.path.join(OUT_DIR, "fig1_boxplot_comparison.png"))
plt.close(fig1)

# ── Figure 2: Bar chart with error bars + significance ──────────────
fig2, ax2 = plt.subplots(figsize=(8, 5.5))

means = [df[df["Group"] == g]["Length_cm"].mean() for g in GROUP_ORDER]
sems  = [df[df["Group"] == g]["Length_cm"].sem()  for g in GROUP_ORDER]

bars = ax2.bar(positions, means, yerr=sems,
               width=0.55, color=[PALETTE[g] for g in GROUP_ORDER],
               alpha=0.82, linewidth=0.8, edgecolor="black",
               error_kw=dict(elinewidth=1, capsize=4, capthick=1,
                             ecolor="black", zorder=4))

y_top2 = max(m + e for m, e in zip(means, sems))
add_cld_labels(ax2, cld, positions, y_top2, y_step=0.006)

ax2.set_xticks(positions)
ax2.set_xticklabels([GROUP_LABELS[g] for g in GROUP_ORDER], fontsize=10)
ax2.set_ylabel("Body Length (cm)", fontsize=12)
ax2.set_xlabel("Treatment Group", fontsize=12)
ax2.set_title("Mean Body Length ± SEM", fontsize=12, pad=10)
ax2.set_ylim(0, y_top2 * 1.18)

ax2.text(0.98, 0.97,
         f"{main_test}: {main_stat_name} = {main_stat:.2f}, "
         f"{'p < 0.001' if main_p < 0.001 else f'p = {main_p:.4f}'}",
         transform=ax2.transAxes, ha="right", va="top",
         fontsize=9, style="italic",
         bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.8", alpha=0.8))
ax2.text(0.98, 0.02, "Different letters indicate significant differences (p < 0.05)",
         transform=ax2.transAxes, ha="right", va="bottom", fontsize=8, color="gray")

plt.tight_layout()
fig2.savefig(os.path.join(OUT_DIR, "fig2_barchart_sem.png"))
plt.close(fig2)

# ── Figure 3: Dose-response curves for G and P ───────────────────────
fig3, axes = plt.subplots(1, 2, figsize=(10, 5), sharey=True)

for ax, trt, color_lo, color_hi, trt_label in zip(
        axes,
        ["G", "P"],
        ["#55A868", "#DD8452"],
        ["#1B7837", "#B22222"],
        ["Treatment G", "Treatment P"]):

    conc_order = [0.0, 0.5, 5.0]
    grp_map    = {0.0: "C", 0.5: f"{trt}-0.5", 5.0: f"{trt}-5"}
    colors_map = {0.0: PALETTE["C"], 0.5: color_lo, 5.0: color_hi}

    xs, ys, es = [], [], []
    for c in conc_order:
        g = grp_map[c]
        v = df[df["Group"] == g]["Length_cm"].values
        xs.append(c)
        ys.append(v.mean())
        es.append(stats.sem(v))

    ax.errorbar(xs, ys, yerr=es,
                fmt="-o", linewidth=1.8, markersize=7,
                color=color_hi, markerfacecolor=color_hi,
                capsize=4, capthick=1, elinewidth=1,
                zorder=4, label="Mean ± SEM")

    # Individual points jittered along x
    for c in conc_order:
        g = grp_map[c]
        v = df[df["Group"] == g]["Length_cm"].values
        jit = np.random.uniform(-0.1, 0.1, len(v))
        ax.scatter(c + jit, v, color=colors_map[c],
                   s=14, alpha=0.4, zorder=3, linewidths=0)

    ax.set_xlabel("PFOA Concentration (mg/L)", fontsize=11)
    if ax == axes[0]:
        ax.set_ylabel("Body Length (cm)", fontsize=11)
    ax.set_title(f"Dose-Response: {trt_label}", fontsize=11)
    ax.set_xticks([0, 0.5, 5.0])
    ax.set_xticklabels(["0\n(Control)", "0.5", "5.0"])

    # Pearson correlation for dose-response
    r, p_r = stats.pearsonr(xs, ys)
    ax.text(0.97, 0.97,
            f"r = {r:.3f}\np = {p_r:.3f}",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=9, style="italic",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.8", alpha=0.8))

fig3.suptitle("Dose-Response of Body Length to PFOA Exposure",
              fontsize=12, y=1.01)
plt.tight_layout()
fig3.savefig(os.path.join(OUT_DIR, "fig3_dose_response.png"))
plt.close(fig3)

# ── Figure 4: Violin plot ─────────────────────────────────────────────
fig4, ax4 = plt.subplots(figsize=(8, 5.5))

vp = ax4.violinplot(group_vals, positions=positions,
                    showmeans=False, showmedians=True,
                    showextrema=True, widths=0.55)

for i, (body, grp) in enumerate(zip(vp['bodies'], GROUP_ORDER)):
    body.set_facecolor(PALETTE[grp])
    body.set_alpha(0.65)
    body.set_linewidth(0.5)

vp['cmedians'].set_color('black')
vp['cmedians'].set_linewidth(2)
for part in ['cbars', 'cmins', 'cmaxes']:
    vp[part].set_color('black')
    vp[part].set_linewidth(0.8)

# Add mean points
for i, (grp, vals) in enumerate(zip(GROUP_ORDER, group_vals)):
    ax4.scatter(positions[i], np.mean(vals),
                marker="D", s=30, color=PALETTE[grp],
                edgecolors="black", linewidths=0.8, zorder=5)

add_cld_labels(ax4, cld, positions,
               max(v.max() for v in group_vals), y_step=0.008)

ax4.set_xticks(positions)
ax4.set_xticklabels([GROUP_LABELS[g] for g in GROUP_ORDER], fontsize=10)
ax4.set_ylabel("Body Length (cm)", fontsize=12)
ax4.set_xlabel("Treatment Group", fontsize=12)
ax4.set_title("Body Length Distribution (Violin Plot)", fontsize=12, pad=10)

legend_elems = [
    Line2D([0], [0], color='black', linewidth=2, label='Median'),
    Line2D([0], [0], marker='D', color='w', markerfacecolor='gray',
           markeredgecolor='black', markersize=6, label='Mean')
]
ax4.legend(handles=legend_elems, loc="lower right", fontsize=9)

plt.tight_layout()
fig4.savefig(os.path.join(OUT_DIR, "fig4_violin_plot.png"))
plt.close(fig4)

# ── Figure 5: Pairwise p-value heatmap ───────────────────────────────
n_grp = len(GROUP_ORDER)
pmat = np.ones((n_grp, n_grp))
for r in posthoc_rows:
    i = GROUP_ORDER.index(r["Group1"])
    j = GROUP_ORDER.index(r["Group2"])
    p = float(r.get("p_adj", r.get("p_adj_bonferroni", 1)))
    pmat[i, j] = p
    pmat[j, i] = p

fig5, ax5 = plt.subplots(figsize=(7, 5.5))
log_p = -np.log10(np.clip(pmat, 1e-10, 1))
np.fill_diagonal(log_p, 0)

im = ax5.imshow(log_p, cmap="YlOrRd", vmin=0, vmax=max(4, log_p.max()))
cbar = plt.colorbar(im, ax=ax5, fraction=0.046, pad=0.04)
cbar.set_label("-log₁₀(p_adj)", fontsize=10)

tl = [GROUP_LABELS[g] for g in GROUP_ORDER]
ax5.set_xticks(range(n_grp)); ax5.set_xticklabels(tl, rotation=30, ha="right")
ax5.set_yticks(range(n_grp)); ax5.set_yticklabels(tl)

for i in range(n_grp):
    for j in range(n_grp):
        if i != j:
            p_ij = pmat[i, j]
            txt = f"{p_ij:.3f}\n{sig_label(p_ij)}"
            ax5.text(j, i, txt, ha="center", va="center",
                     fontsize=7.5,
                     color="white" if log_p[i, j] > 2 else "black")
        else:
            ax5.text(j, i, "—", ha="center", va="center",
                     fontsize=10, color="gray")

ax5.set_title(f"Post-hoc Adjusted p-value Matrix\n({posthoc_method})",
              fontsize=11, pad=10)
plt.tight_layout()
fig5.savefig(os.path.join(OUT_DIR, "fig5_pvalue_heatmap.png"))
plt.close(fig5)

# ─────────────────────────────────────────────
# 10. Comprehensive summary report
# ─────────────────────────────────────────────
report_lines = [
    "=" * 65,
    "PFOA BODY LENGTH ANALYSIS — STATISTICAL REPORT",
    "=" * 65,
    "",
    "DATA SUMMARY",
    "-" * 40,
    f"  Total measurements (raw):  {len(raw)}",
    f"  Outliers removed (IQR):    {len(outlier_log)}",
    f"  Clean measurements:        {len(df)}",
    "",
    "DESCRIPTIVE STATISTICS (Length, cm)",
    "-" * 40,
]
for g in GROUP_ORDER:
    row = desc_grp.loc[g]
    report_lines.append(
        f"  {GROUP_LABELS[g]:<14}  n={int(row['n'])}  "
        f"Mean={row['Mean']:.4f}  SD={row['SD']:.4f}  "
        f"SE={row['SE']:.4f}  CV%={row['CV%']:.2f}"
    )

report_lines += [
    "",
    "ASSUMPTION TESTS",
    "-" * 40,
    "  Shapiro-Wilk normality test:",
]
for r in sw_results:
    report_lines.append(
        f"    {r['Group']:<8}  W={r['W_statistic']:.4f}  "
        f"p={r['p_value']:.4f}  Normal: {r['Normal (p>0.05)']}"
    )
report_lines += [
    f"  Levene's test: F={lev_stat:.4f}  p={lev_p:.4f}  "
    f"Equal variances: {lev_result['Equal_variances (p>0.05)']}",
    "",
    "MAIN TEST",
    "-" * 40,
    f"  Method: {main_test}",
    f"  {main_stat_name} = {main_stat:.4f}",
    f"  p-value = {main_p:.6f}",
    f"  Effect size ({effect_name}) = {effect_val:.4f}",
    f"  Significant overall: {sig}",
    "",
    "POST-HOC COMPARISONS",
    "-" * 40,
    f"  Method: {posthoc_method}",
]
for r in posthoc_rows:
    p = float(r.get("p_adj", r.get("p_adj_bonferroni", 1)))
    report_lines.append(
        f"  {r['Group1']:<8} vs {r['Group2']:<8}  "
        f"Δmean={r['Mean_diff']:+.4f}  p_adj={p:.4f}  {sig_label(p)}"
    )

report_lines += [
    "",
    "COMPACT LETTER DISPLAY",
    "-" * 40,
]
for g in GROUP_ORDER:
    report_lines.append(f"  {GROUP_LABELS[g]:<14}: {cld[g]}")

report_lines += [
    "",
    "OUTPUT FILES",
    "-" * 40,
    "  descriptive_stats_per_group.csv",
    "  descriptive_stats_per_replicate.csv",
    "  normality_shapiro_wilk.csv",
    "  levene_test.csv",
    "  main_test_result.csv",
    f"  {'posthoc_tukey_hsd.csv' if use_parametric else 'posthoc_mannwhitney_bonferroni.csv'}",
    "  effect_sizes_cohens_d.csv",
    "  fig1_boxplot_comparison.png",
    "  fig2_barchart_sem.png",
    "  fig3_dose_response.png",
    "  fig4_violin_plot.png",
    "  fig5_pvalue_heatmap.png",
    "",
    "=" * 65,
]

report_text = "\n".join(report_lines)
print("\n" + report_text)

with open(os.path.join(OUT_DIR, "statistical_report.txt"), "w",
          encoding="utf-8") as f:
    f.write(report_text)

print(f"\nAll results saved to: {OUT_DIR}")
