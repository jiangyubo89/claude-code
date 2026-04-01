"""
Generate a publication-ready boxplot for body length from 体长.csv.

Input (wide CSV): D:/桌面/claude/idtrackerai/收样数据/体长.csv
Output:
  - fig1_body_length_boxplot.png
  - this script (regen_fig1.py)
"""

import io
import os
import sys
import warnings

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scipy.stats as stats

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
warnings.filterwarnings("ignore")

OUT_DIR = os.path.dirname(os.path.abspath(__file__))
CSV_PATH = os.path.join(OUT_DIR, "体长.csv")

plt.rcParams.update(
    {
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
    }
)

GROUP_ORDER = ["C", "50", "250", "1250"]
X_LABELS = ["Control", "50 mg/kg", "250 mg/kg", "1250 mg/kg"]
PALETTE = {"C": "#4C72B0", "50": "#55A868", "250": "#DD8452", "1250": "#B22222"}


def sig_label(p: float) -> str:
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def read_csv_fallback(path: str) -> pd.DataFrame:
    last_err: Exception | None = None
    for enc in ("utf-8-sig", "utf-8", "gbk", "cp936"):
        try:
            return pd.read_csv(path, header=0, encoding=enc)
        except Exception as e:  # noqa: BLE001
            last_err = e
    raise RuntimeError(f"Failed to read CSV with common encodings: {path}. Last error: {last_err}")


raw = read_csv_fallback(CSV_PATH)

# Expected wide format: first column is label, remaining columns are groups
first_col = raw.columns[0]
wide = raw.drop(columns=[first_col], errors="ignore").copy()
wide = wide.applymap(lambda x: str(x).strip() if isinstance(x, str) else x)

df_long = wide.melt(var_name="Group", value_name="Length_cm")
df_long["Length_cm"] = pd.to_numeric(df_long["Length_cm"], errors="coerce")
df_long = df_long.dropna(subset=["Length_cm"])
df_long["Group"] = df_long["Group"].astype(str).str.strip()
df_long = df_long[df_long["Group"].isin(GROUP_ORDER)].copy()

group_vals = [df_long[df_long["Group"] == g]["Length_cm"].values for g in GROUP_ORDER]

print("Final n per group:")
for g in GROUP_ORDER:
    print(f"  {g}: n={len(df_long[df_long['Group']==g])}")

# Single planned comparison: Control vs 1250 (t-test)
ctrl_vals = df_long[df_long["Group"] == "C"]["Length_cm"].values
tg = "1250"
tv_1250 = df_long[df_long["Group"] == tg]["Length_cm"].values
t_res = stats.ttest_ind(ctrl_vals, tv_1250, equal_var=False, nan_policy="omit")
p_1250 = float(t_res.pvalue) if np.isfinite(t_res.pvalue) else 1.0
sig_map = {tg: sig_label(p_1250)}

print("Planned comparison (t-test): Control vs 1250")
print(f"  1250: {sig_label(p_1250)} (p={p_1250:.6f})")

fig, ax = plt.subplots(figsize=(7.6, 5.2))
positions = np.arange(len(GROUP_ORDER))

bp = ax.boxplot(
    group_vals,
    positions=positions,
    widths=0.45,
    patch_artist=True,
    whis=[0, 100],
    medianprops=dict(color="black", linewidth=1.8),
    whiskerprops=dict(linewidth=0.9),
    capprops=dict(linewidth=0.9),
    flierprops=dict(marker="", linestyle="none"),
    showfliers=False,
)
for patch, grp in zip(bp["boxes"], GROUP_ORDER):
    patch.set_facecolor(PALETTE.get(grp, "#999999"))
    patch.set_alpha(0.72)
    patch.set_linewidth(0.8)

# Jittered individual points
np.random.seed(42)
for i, (grp, vals) in enumerate(zip(GROUP_ORDER, group_vals)):
    jitter = np.random.uniform(-0.18, 0.18, size=len(vals))
    ax.scatter(
        positions[i] + jitter,
        vals,
        color=PALETTE.get(grp, "#999999"),
        s=18,
        alpha=0.55,
        zorder=3,
        linewidths=0,
    )

# Significance bracket: Control vs 1250 only
tops = [v.max() for v in group_vals if len(v)]
y_data_max = max(tops) if tops else 0.0
ctrl_pos = positions[0]

if tg in GROUP_ORDER:
    tg_idx = GROUP_ORDER.index(tg)
    tg_pos = positions[tg_idx]
    label = sig_map.get(tg, "ns")
    y_br = y_data_max + 0.08
    tick = 0.02

    ax.plot([ctrl_pos, tg_pos], [y_br, y_br], color="black", linewidth=0.9, clip_on=False)
    ax.plot([ctrl_pos, ctrl_pos], [y_br - tick, y_br], color="black", linewidth=0.9, clip_on=False)
    ax.plot([tg_pos, tg_pos], [y_br - tick, y_br], color="black", linewidth=0.9, clip_on=False)
    ax.text(
        (ctrl_pos + tg_pos) / 2,
        y_br + 0.01,
        label,
        ha="center",
        va="bottom",
        fontsize=10 if label != "ns" else 9,
        fontweight="bold" if label != "ns" else "normal",
    )

ax.set_xticks(positions)
ax.set_xticklabels(X_LABELS)
ax.set_ylabel("Body Length (cm)")
ax.set_xlabel("Group")
ax.set_title("Body Length Comparison", pad=26)

# Annotation aligned with title
p_str = "p < 0.001" if p_1250 < 0.001 else f"p = {p_1250:.4f}"
ax.text(
    0.5,
    1.02,
    f"t-test: {p_str}",
    transform=ax.transAxes,
    ha="center",
    va="bottom",
    fontsize=9,
    style="italic",
)

y_top_all = y_data_max + 0.18
ax.set_ylim(ax.get_ylim()[0], y_top_all)

plt.tight_layout()
out_path = os.path.join(OUT_DIR, "fig1_body_length_boxplot.png")
fig.savefig(out_path)
plt.close(fig)
print(f"Saved: {out_path}")
