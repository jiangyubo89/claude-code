"""
Generate main panel figure and bubble plots (v2)
Synced with dual-tier DEP classification from proteomics_analysis.py
"""
import pandas as pd, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import os

try:
    from adjustText import adjust_text
    HAS_ADJUST = True
except ImportError:
    HAS_ADJUST = False

plt.rcParams['font.family'] = 'Arial'

OUT = r"D:\桌面\proteomics_output"
os.makedirs(OUT, exist_ok=True)

# ── Load data and reproduce classification ───────────────────────────────────
INPUT_CSV = os.path.join(OUT, "2026.3.12 pr-ptx3 60h 蛋白质谱3vs3.csv")
df = pd.read_csv(INPUT_CSV)
df.columns = df.columns.str.strip()
num_cols = ['CTRL-1', 'CTRL-2', 'CTRL-3', 'PTX3-1', 'PTX3-2', 'PTX3-3',
            'log2FC', 'P.value', 'FC', 'mean c', 'mean ptx3']
for col in num_cols:
    df[col] = pd.to_numeric(df[col], errors='coerce')
df['-log10P'] = -np.log10(df['P.value'].clip(lower=1e-15))

ctrl_cols = ['CTRL-1', 'CTRL-2', 'CTRL-3']
ptx3_cols = ['PTX3-1', 'PTX3-2', 'PTX3-3']

# BH FDR
from statsmodels.stats.multitest import multipletests
valid_mask = df['P.value'].notna()
df['FDR'] = np.nan
_, fdr_vals, _, _ = multipletests(df.loc[valid_mask, 'P.value'], method='fdr_bh')
df.loc[valid_mask, 'FDR'] = fdr_vals

# On/off detection
ONOFF_THRESH = 0.01
FC_CUTOFF = 1.0
P_STRINGENT = 0.01
P_EXPLORATORY = 0.05

ctrl_max = df[ctrl_cols].max(axis=1)
ptx3_max = df[ptx3_cols].max(axis=1)
df['is_onoff'] = (ctrl_max < ONOFF_THRESH) | (ptx3_max < ONOFF_THRESH)

def classify_protein(row):
    if pd.isna(row['P.value']) or pd.isna(row['log2FC']):
        return 'NS'
    if row['is_onoff']:
        if row['P.value'] < P_EXPLORATORY:
            return 'On/Off'
        return 'NS'
    if abs(row['log2FC']) >= FC_CUTOFF:
        if row['P.value'] < P_STRINGENT:
            return 'Stringent_Up' if row['log2FC'] > 0 else 'Stringent_Down'
        elif row['P.value'] < P_EXPLORATORY:
            return 'Exploratory_Up' if row['log2FC'] > 0 else 'Exploratory_Down'
    return 'NS'

df['Class'] = df.apply(classify_protein, axis=1)

stringent_up   = df[df['Class'] == 'Stringent_Up']
stringent_down = df[df['Class'] == 'Stringent_Down']
exploratory_up   = df[df['Class'] == 'Exploratory_Up']
exploratory_down = df[df['Class'] == 'Exploratory_Down']
all_up   = df[df['Class'].isin(['Stringent_Up', 'Exploratory_Up'])]
all_down = df[df['Class'].isin(['Stringent_Down', 'Exploratory_Down'])]
onoff_df = df[df['Class'] == 'On/Off']

# ── Load enrichment CSVs ─────────────────────────────────────────────────────
def load_enr(key):
    f = f'{OUT}/enrichment_{key}.csv'
    return pd.read_csv(f) if os.path.exists(f) else pd.DataFrame()

kegg_up = load_enr('KEGG_up')
kegg_dn = load_enr('KEGG_down')
gobp_up = load_enr('GOBP_up')
gobp_dn = load_enr('GOBP_down')
gocc_up = load_enr('GOCC_up')
gomf_up = load_enr('GOMF_up')
wiki_up = load_enr('WikiPath_up')


# ══════════════════════════════════════════════════════════════════════════════
# MAIN 5-PANEL FIGURE
# ══════════════════════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(18, 14))
gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.40, wspace=0.35)

ax_vol  = fig.add_subplot(gs[0, 0])
ax_heat = fig.add_subplot(gs[0:2, 1])
ax_kegg = fig.add_subplot(gs[0, 2])
ax_go   = fig.add_subplot(gs[1, 0])
ax_kegg2 = fig.add_subplot(gs[1, 2])

# ── A. Volcano ───────────────────────────────────────────────────────────────
ns = df[df['Class'] == 'NS']
ax_vol.scatter(ns['log2FC'], ns['-log10P'], s=5, c='#CCCCCC',
               alpha=0.4, linewidths=0, rasterized=True)
ax_vol.scatter(exploratory_up['log2FC'], exploratory_up['-log10P'],
               s=12, c='#E8927C', alpha=0.7, linewidths=0)
ax_vol.scatter(exploratory_down['log2FC'], exploratory_down['-log10P'],
               s=12, c='#7FAED4', alpha=0.7, linewidths=0)
ax_vol.scatter(stringent_up['log2FC'], stringent_up['-log10P'],
               s=18, c='#D73027', alpha=0.9, linewidths=0,
               label=f'Up ({len(all_up)})')
ax_vol.scatter(stringent_down['log2FC'], stringent_down['-log10P'],
               s=18, c='#4575B4', alpha=0.9, linewidths=0,
               label=f'Down ({len(all_down)})')

# On/off triangles
onoff_up_v = onoff_df[onoff_df['log2FC'] > 0]
onoff_dn_v = onoff_df[onoff_df['log2FC'] < 0]
ax_vol.scatter(onoff_up_v['log2FC'], onoff_up_v['-log10P'],
               s=14, c='#D73027', alpha=0.5, marker='^', linewidths=0.3, edgecolors='k')
ax_vol.scatter(onoff_dn_v['log2FC'], onoff_dn_v['-log10P'],
               s=14, c='#4575B4', alpha=0.5, marker='^', linewidths=0.3, edgecolors='k')

ax_vol.axhline(-np.log10(0.05), ls='--', lw=0.8, c='#555555')
ax_vol.axvline(0, ls='-', lw=0.5, c='#888888')
ax_vol.axvline(FC_CUTOFF, ls=':', lw=0.7, c='#999999')
ax_vol.axvline(-FC_CUTOFF, ls=':', lw=0.7, c='#999999')

# Label top 5+5 (quantitative) + PTX3 (experimental target)
top_label = pd.concat([all_up.nlargest(5, '-log10P'), all_down.nlargest(5, '-log10P')])
ptx3_row = onoff_df[onoff_df['Gene'] == 'PTX3']
top_label = pd.concat([top_label, ptx3_row])
texts = []
for _, row in top_label.iterrows():
    col = '#A52A2A' if row['log2FC'] > 0 else '#1A3A6B'
    gene_name = row['Gene']
    if gene_name == 'PTX3':
        t = ax_vol.annotate('PTX3 *', (row['log2FC'], row['-log10P']),
                            fontsize=7.5, color='#8B0000', fontweight='bold',
                            xytext=(3, 2), textcoords='offset points',
                            bbox=dict(boxstyle='round,pad=0.1', fc='#FFFFCC', ec='#8B0000', alpha=0.8, lw=0.5))
    else:
        t = ax_vol.annotate(gene_name, (row['log2FC'], row['-log10P']),
                            fontsize=7, color=col, fontweight='bold',
                            xytext=(3, 2), textcoords='offset points')
    texts.append(t)
if HAS_ADJUST and texts:
    adjust_text(texts, ax=ax_vol, arrowprops=dict(arrowstyle='-', color='gray', lw=0.4))

ax_vol.set_xlabel('log2(Fold Change)', fontsize=10)
ax_vol.set_ylabel('-log10(P value)', fontsize=10)
ax_vol.set_title('A  Volcano Plot', fontsize=11, fontweight='bold', loc='left')
ax_vol.legend(fontsize=8, frameon=True, loc='upper left')
ax_vol.spines[['top', 'right']].set_visible(False)
ax_vol.tick_params(labelsize=9)

# ── B. Heatmap ───────────────────────────────────────────────────────────────
n_up_h = min(10, len(all_up))
n_dn_h = min(10, len(all_down))
top_u = all_up.nlargest(n_up_h, '-log10P')
top_d = all_down.nlargest(n_dn_h, '-log10P')
heat_df = pd.concat([top_u, top_d]).set_index('Gene')
mat = heat_df[ctrl_cols + ptx3_cols].astype(float)
mat_z = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1) + 1e-9, axis=0)

sns.heatmap(mat_z, cmap='RdBu_r', center=0, vmin=-2.5, vmax=2.5,
            linewidths=0.3, linecolor='white',
            cbar_kws={'label': 'Z-score', 'shrink': 0.5, 'pad': 0.02},
            xticklabels=['C1', 'C2', 'C3', 'P1', 'P2', 'P3'],
            yticklabels=True, ax=ax_heat)
ax_heat.set_title(f'B  Top {len(heat_df)} DEPs Heatmap (Z-score)',
                  fontsize=11, fontweight='bold', loc='left')
ax_heat.tick_params(axis='y', labelsize=8)
ax_heat.tick_params(axis='x', labelsize=9, rotation=0)
up_genes_set = set(top_u['Gene'])
for lbl in ax_heat.get_yticklabels():
    lbl.set_color('#B22222' if lbl.get_text() in up_genes_set else '#2B4D9E')

# ── C/D/E. Enrichment bar (using -log10 Adj.P) ──────────────────────────────
def plot_enr_ax(ax, edf, title, top_n=10):
    if edf is None or edf.empty:
        ax.text(0.5, 0.5, 'No significant terms\n(n_genes >= 3 filter)',
                ha='center', va='center', transform=ax.transAxes, fontsize=10)
        ax.set_title(title, fontsize=10, fontweight='bold', loc='left')
        return
    # Filter n_genes >= 3 if column exists
    if 'n_genes' in edf.columns:
        edf = edf[edf['n_genes'] >= 3]
    elif 'Genes' in edf.columns:
        edf = edf.copy()
        edf['n_genes'] = edf['Genes'].str.split(';').str.len()
        edf = edf[edf['n_genes'] >= 3]

    if edf.empty:
        ax.text(0.5, 0.5, 'No terms with n >= 3 genes',
                ha='center', va='center', transform=ax.transAxes, fontsize=10)
        ax.set_title(title, fontsize=10, fontweight='bold', loc='left')
        return

    dplot = edf.head(top_n).copy()
    dplot['Term_s'] = dplot['Term'].str.split(' Homo').str[0].str[:48]
    dplot['-log10AdjP'] = -np.log10(dplot['Adj.P'].clip(lower=1e-10))
    dplot = dplot.sort_values('-log10AdjP')  # bottom = least significant

    norm = plt.Normalize(dplot['-log10AdjP'].min(), dplot['-log10AdjP'].max())
    colors_bar = plt.cm.YlOrRd(norm(dplot['-log10AdjP'].values))

    bars = ax.barh(range(len(dplot)), dplot['-log10AdjP'],
                   color=colors_bar, height=0.65, edgecolor='none')

    # Fade insignificant
    for bar, adjp in zip(bars, dplot['Adj.P']):
        if adjp >= 0.05:
            bar.set_alpha(0.4)

    # Significance line
    ax.axvline(-np.log10(0.05), ls='--', lw=0.8, c='#333333')

    ax.set_yticks(range(len(dplot)))
    ax.set_yticklabels(dplot['Term_s'], fontsize=8.5)
    ax.set_xlabel('-log10(Adj.P)', fontsize=9)
    ax.set_title(title, fontsize=10, fontweight='bold', loc='left')
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(labelsize=8)

    n_genes_col = 'n_genes' if 'n_genes' in dplot.columns else None
    for i, (_, row) in enumerate(dplot.iterrows()):
        n = row.get('n_genes', len(str(row.get('Genes', '')).split(';')))
        ax.text(0.2, i, f"n={n}", va='center', fontsize=6.5,
                color='white', fontweight='bold')

plot_enr_ax(ax_kegg, kegg_up, 'C  KEGG Pathways – Up')
plot_enr_ax(ax_go, gobp_up, 'D  GO Biological Process – Up', top_n=12)
plot_enr_ax(ax_kegg2, kegg_dn, 'E  KEGG Pathways – Down')

plt.suptitle('PTX3 Treatment vs Control (60 h): Quantitative Proteomics (n=3+3)\n'
             f'DEPs: {len(all_up)} up + {len(all_down)} down '
             f'(p<0.05, |log2FC|>={FC_CUTOFF}, excl On/Off)',
             fontsize=13, fontweight='bold', y=1.02)
plt.savefig(f'{OUT}/00_MAIN_PANEL.pdf', dpi=300, bbox_inches='tight')
plt.savefig(f'{OUT}/00_MAIN_PANEL.png', dpi=300, bbox_inches='tight')
plt.close()
print("Main panel saved")


# ══════════════════════════════════════════════════════════════════════════════
# BUBBLE PLOT: GO BP up & down  (using -log10 Adj.P)
# ══════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(16, 9))

def bubble_enr(ax, edf, title, top_n=15):
    if edf is None or edf.empty:
        ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
        ax.set_title(title, fontsize=11, fontweight='bold')
        return

    # Filter n_genes >= 3
    dplot = edf.copy()
    if 'n_genes' not in dplot.columns:
        dplot['n_genes'] = dplot['Genes'].str.split(';').str.len()
    dplot = dplot[dplot['n_genes'] >= 3]

    if dplot.empty:
        ax.text(0.5, 0.5, 'No terms with n >= 3', ha='center', va='center',
                transform=ax.transAxes, fontsize=10)
        ax.set_title(title, fontsize=11, fontweight='bold')
        return

    dplot = dplot.head(top_n).copy()
    dplot['Term_s'] = dplot['Term'].str.split(' Homo').str[0].str[:55]
    dplot['-log10AdjP'] = -np.log10(dplot['Adj.P'].clip(lower=1e-10))
    dplot = dplot.sort_values('-log10AdjP')  # bottom = least significant

    sc = ax.scatter(dplot['-log10AdjP'], range(len(dplot)),
                    s=dplot['n_genes'] * 25,
                    c=dplot['-log10AdjP'],
                    cmap='YlOrRd',
                    vmin=0, vmax=max(4, dplot['-log10AdjP'].max()),
                    alpha=0.9, edgecolors='gray', lw=0.5)
    ax.set_yticks(range(len(dplot)))
    ax.set_yticklabels(dplot['Term_s'], fontsize=9)
    ax.set_xlabel('-log10(Adjusted P value)', fontsize=12)
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.spines[['top', 'right']].set_visible(False)

    # Significance line
    ax.axvline(-np.log10(0.05), ls='--', lw=0.8, c='#333333', label='Adj.P = 0.05')
    ax.legend(fontsize=8, loc='lower right')

    cb = plt.colorbar(sc, ax=ax, label='-log10(Adj.P)', shrink=0.6)
    cb.ax.tick_params(labelsize=9)

    # Size legend
    for size_n in [3, 5, 10]:
        if size_n <= dplot['n_genes'].max():
            ax.scatter([], [], s=size_n * 25, c='gray', alpha=0.7, label=f'{size_n} genes')
    ax.legend(title='Gene count', fontsize=8, title_fontsize=8, loc='lower right')

bubble_enr(axes[0], gobp_up,
           'GO Biological Process\nUp-regulated (PTX3 vs Control, 60 h)')
bubble_enr(axes[1], gobp_dn,
           'GO Biological Process\nDown-regulated (PTX3 vs Control, 60 h)')

plt.tight_layout()
plt.savefig(f'{OUT}/13_GOBP_bubble.pdf', dpi=300, bbox_inches='tight')
plt.savefig(f'{OUT}/13_GOBP_bubble.png', dpi=300, bbox_inches='tight')
plt.close()
print("Bubble plot saved")


# ══════════════════════════════════════════════════════════════════════════════
# KEGG combined dotplot (up + down side by side)  using -log10(Adj.P)
# ══════════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(11, 9))

def prep_kegg(edf, direction, top_n=10):
    if edf is None or edf.empty:
        return pd.DataFrame()
    d = edf.copy()
    if 'n_genes' not in d.columns:
        d['n_genes'] = d['Genes'].str.split(';').str.len()
    d = d[d['n_genes'] >= 3].head(top_n).copy()
    d['direction'] = direction
    d['-log10AdjP'] = -np.log10(d['Adj.P'].clip(lower=1e-10))
    sign = 1 if direction == 'Up' else -1
    d['signed_logP'] = sign * d['-log10AdjP']
    return d

kegg_u10 = prep_kegg(kegg_up, 'Up')
kegg_d10 = prep_kegg(kegg_dn, 'Down')

if not kegg_u10.empty or not kegg_d10.empty:
    comb = pd.concat([kegg_u10, kegg_d10])
    comb['Term_s'] = comb['Term'].str.split(' Homo').str[0].str[:50]
    comb = comb.sort_values('signed_logP')

    colors_bar = ['#D73027' if d == 'Up' else '#4575B4' for d in comb['direction']]

    bars = ax.barh(range(len(comb)), comb['signed_logP'],
                   color=colors_bar, height=0.65, edgecolor='none')

    # Fade insignificant
    for bar, adjp in zip(bars, comb['Adj.P']):
        if adjp >= 0.05:
            bar.set_alpha(0.4)

    ax.set_yticks(range(len(comb)))
    ax.set_yticklabels(comb['Term_s'], fontsize=9)
    ax.axvline(0, color='black', lw=0.8)

    # Significance lines
    sig_line = -np.log10(0.05)
    ax.axvline(sig_line, ls='--', lw=0.7, c='#999999')
    ax.axvline(-sig_line, ls='--', lw=0.7, c='#999999')

    ax.set_xlabel('-log10(Adj.P)  (positive=Up, negative=Down)', fontsize=11)
    ax.set_title('KEGG Pathway Enrichment\nPTX3 vs Control (60 h)',
                 fontsize=12, fontweight='bold')
    ax.spines[['top', 'right']].set_visible(False)

    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='#D73027', label='Up-regulated'),
                       Patch(facecolor='#4575B4', label='Down-regulated')]
    ax.legend(handles=legend_elements, fontsize=10, loc='lower right')
else:
    ax.text(0.5, 0.5, 'No KEGG terms with n >= 3 genes',
            ha='center', va='center', transform=ax.transAxes, fontsize=12)

plt.tight_layout()
plt.savefig(f'{OUT}/14_KEGG_combined.pdf', dpi=300, bbox_inches='tight')
plt.savefig(f'{OUT}/14_KEGG_combined.png', dpi=300, bbox_inches='tight')
plt.close()
print("KEGG combined saved")


# ══════════════════════════════════════════════════════════════════════════════
# WRITE COMPREHENSIVE REPORT (data-driven, no hardcoded numbers)
# ══════════════════════════════════════════════════════════════════════════════
def fmt_enr(edf, top_n=8):
    if edf is None or edf.empty:
        return "  (No terms with n >= 3 genes)\n"
    lines = []
    for _, row in edf.head(top_n).iterrows():
        term = str(row['Term']).split(' Homo')[0][:55]
        n = row.get('n_genes', len(str(row.get('Genes', '')).split(';')))
        sig_mark = "**" if row['Adj.P'] < 0.05 else "  "
        lines.append(f"  {sig_mark}{row['Rank']:>2}. {term:<55} adj.p={row['Adj.P']:.2e} n={n}")
    return "\n".join(lines)

# Sample QC stats
log_mat = df[ctrl_cols + ptx3_cols].clip(lower=0.001).apply(np.log2)
corr = log_mat.corr()
ctrl_intra = [corr.loc[a, b] for i, a in enumerate(ctrl_cols)
              for j, b in enumerate(ctrl_cols) if i < j]
ptx3_intra = [corr.loc[a, b] for i, a in enumerate(ptx3_cols)
              for j, b in enumerate(ptx3_cols) if i < j]

report = f"""
{"="*72}
PROTEOMICS ANALYSIS REPORT — PTX3 vs CONTROL (60 h)  [v2]
Quantitative Label-Free Proteomics | 3+3 Biological Replicates
Generated: 2026-03-14
{"="*72}

1. EXPERIMENTAL OVERVIEW
{"─"*60}
  Comparison : PTX3-treated cells vs untreated Control
  Duration   : 60 hours
  Design     : 3 biological replicates per group (CTRL-1/2/3, PTX3-1/2/3)
  Platform   : Label-free quantification (LFQ) mass spectrometry
  Database   : UniProt/Swiss-Prot (Human)

2. DATA SUMMARY & DEP CLASSIFICATION
{"─"*60}
  Proteins quantified : {len(df):,}
  log2FC range        : {df['log2FC'].min():.2f} to {df['log2FC'].max():.2f}

  Classification criteria:
    FC cutoff      : |log2FC| >= {FC_CUTOFF}
    On/Off thresh  : max replicate intensity < {ONOFF_THRESH}
    FDR method     : Benjamini-Hochberg

  STRINGENT DEPs (p < {P_STRINGENT}, |log2FC| >= {FC_CUTOFF}, excl on/off):
    Up-regulated   : {len(stringent_up)}
    Down-regulated : {len(stringent_down)}
  EXPLORATORY DEPs (p < {P_EXPLORATORY}, |log2FC| >= {FC_CUTOFF}, excl on/off):
    Up-regulated   : {len(exploratory_up)}
    Down-regulated : {len(exploratory_down)}
  TOTAL quantitative DEPs: {len(all_up) + len(all_down)} ({len(all_up)} up, {len(all_down)} down)
  On/Off proteins (p < 0.05): {len(onoff_df)}

  NOTE: BH FDR is very conservative with n=3 per group.
  Only {(df['FDR'] < 0.05).sum()} proteins pass FDR < 0.05 (all on/off artifacts).

3. QUALITY CONTROL
{"─"*60}
  Intra-group Pearson r (Ctrl)  : {np.mean(ctrl_intra):.4f} +/- {np.std(ctrl_intra):.4f}
  Intra-group Pearson r (PTX3)  : {np.mean(ptx3_intra):.4f} +/- {np.std(ptx3_intra):.4f}

4. TOP STRINGENT DEPs
{"─"*60}
  [UP-REGULATED — top by log2FC, p < {P_STRINGENT}]
  {"Gene":<14} {"log2FC":>8}  {"P value":>10}  {"FDR":>10}  {"FC":>8}
  {"-"*62}
"""

for _, row in stringent_up.nlargest(15, 'log2FC').iterrows():
    fdr_s = f"{row['FDR']:.4f}" if pd.notna(row['FDR']) else "NA"
    report += f"  {row['Gene']:<14} {row['log2FC']:>8.3f}  {row['P.value']:>10.4f}  {fdr_s:>10}  {row['FC']:>8.2f}\n"

report += f"""
  [DOWN-REGULATED — top by |log2FC|, p < {P_STRINGENT}]
  {"Gene":<14} {"log2FC":>8}  {"P value":>10}  {"FDR":>10}  {"FC":>8}
  {"-"*62}
"""

for _, row in stringent_down.nsmallest(15, 'log2FC').iterrows():
    fdr_s = f"{row['FDR']:.4f}" if pd.notna(row['FDR']) else "NA"
    report += f"  {row['Gene']:<14} {row['log2FC']:>8.3f}  {row['P.value']:>10.4f}  {fdr_s:>10}  {row['FC']:>8.2f}\n"

report += f"""
5. ON/OFF PROTEINS
{"─"*60}
  Proteins detected in one condition but not the other (max replicate
  intensity < {ONOFF_THRESH} in at least one group). These show extreme log2FC
  values that are artifacts of near-zero denominators.
  Total: {len(onoff_df)} (Gained: {len(onoff_df[onoff_df['log2FC']>0])}, Lost: {len(onoff_df[onoff_df['log2FC']<0])})
  These are EXCLUDED from pathway enrichment analysis.

6. PATHWAY ENRICHMENT ANALYSIS (Enrichr)
{"─"*60}
  Gene lists: {len(all_up)} up-regulated, {len(all_down)} down-regulated
  Minimum gene overlap: n >= 3
  ** = Adj.P < 0.05

  [KEGG – Up-regulated]
{fmt_enr(kegg_up)}

  [KEGG – Down-regulated]
{fmt_enr(kegg_dn)}

  [GO Biological Process – Up-regulated]
{fmt_enr(gobp_up)}

  [GO Biological Process – Down-regulated]
{fmt_enr(gobp_dn)}

  [GO Molecular Function – Up-regulated]
{fmt_enr(gomf_up)}

  [GO Cellular Component – Up-regulated]
{fmt_enr(gocc_up)}

  [WikiPathways – Up-regulated]
{fmt_enr(wiki_up)}

7. STATISTICAL CONSIDERATIONS
{"─"*60}
  - P values are unadjusted two-sample t-tests (from original data)
  - BH FDR computed but too conservative for n=3 (see Section 2)
  - Dual-tier approach balances discovery vs false positive control
  - On/off proteins excluded from enrichment to avoid inflated FC bias
  - Enrichment terms filtered: minimum 3 overlapping genes
  - Enrichment x-axis uses -log10(Adj.P) — Enrichr BH-corrected
  - Adj.P < 0.05 threshold line shown on all enrichment plots

8. RECOMMENDED FOLLOW-UP
{"─"*60}
  1. Validate top hits by Western blot or targeted mass spectrometry
  2. Assess apoptosis: Annexin V/PI staining + caspase activity assay
  3. Innate immune signaling: NF-kB luciferase reporter
  4. Mitochondrial function: Seahorse XF respirometry
  5. Cell cycle: BrdU/EdU incorporation or FACS-based analysis
  6. PPI network: STRING/Cytoscape analysis of DEPs

{"="*72}
OUTPUT FILES (in {OUT}):
  00_MAIN_PANEL.pdf/png       - 5-panel summary figure
  01_volcano_plot.pdf/png     - Volcano plot (dual-tier + on/off)
  02_heatmap_top_DEPs.pdf/png - Heatmap of top quantitative DEPs
  03_bar_chart_top_DEPs.pdf/png - Top up/down bar charts
  04_sample_correlation.pdf/png - Sample correlation heatmap
  05_distribution_plots.pdf/png - P-value & FC distributions
  06_scatter_ctrl_vs_ptx3.pdf/png - Intensity scatter plot
  07-12_enrichment.pdf/png    - Pathway enrichment bar charts
  13_GOBP_bubble.pdf/png      - GO BP bubble plot
  14_KEGG_combined.pdf/png    - KEGG combined bar
  DEPs_table.csv              - Quantitative DEPs ({len(all_up)+len(all_down)} proteins)
  OnOff_proteins_table.csv    - On/off proteins ({len(onoff_df)} proteins)
  full_protein_table.csv      - All {len(df):,} proteins
  enrichment_*.csv            - Pathway enrichment tables
  ANALYSIS_REPORT.txt         - This report
{"="*72}
"""

with open(f'{OUT}/ANALYSIS_REPORT.txt', 'w', encoding='utf-8') as f:
    f.write(report)
print("Report saved!")
print("=== ALL DONE ===")
