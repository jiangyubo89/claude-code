"""
PTX3 Proteomics Analysis - Publication Quality (v2)
3 Control vs 3 PTX3-treated (60h) - Human Cells
Fixes: BH FDR correction, on/off artifact handling, dual-tier DEP classification,
       enrichment quality filters, improved figure labeling
"""

import csv, math, json, os, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import seaborn as sns
from scipy import stats
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests
import urllib.request, urllib.parse
warnings.filterwarnings('ignore')

try:
    from adjustText import adjust_text
    HAS_ADJUST = True
except ImportError:
    HAS_ADJUST = False
    print("WARNING: adjustText not installed; labels may overlap")

# ── Global plot settings ─────────────────────────────────────────────────────
plt.rcParams['font.family'] = 'Arial'

# ── Output directory ─────────────────────────────────────────────────────────
OUT = r"D:\桌面\proteomics_output"
os.makedirs(OUT, exist_ok=True)

# ── Load data ────────────────────────────────────────────────────────────────
INPUT_CSV = os.path.join(OUT, "2026.3.12 pr-ptx3 60h 蛋白质谱3vs3.csv")
df = pd.read_csv(INPUT_CSV)
df.columns = df.columns.str.strip()

ctrl_cols = ['CTRL-1', 'CTRL-2', 'CTRL-3']
ptx3_cols = ['PTX3-1', 'PTX3-2', 'PTX3-3']

# Numeric conversion
for col in ctrl_cols + ptx3_cols + ['mean c', 'mean ptx3', 'FC', 'log2FC', 'P.value']:
    df[col] = pd.to_numeric(df[col], errors='coerce')

# ══════════════════════════════════════════════════════════════════════════════
# DEP CLASSIFICATION  (replaces original Sig-based approach)
# ══════════════════════════════════════════════════════════════════════════════
# 1. Compute BH-corrected FDR
valid_mask = df['P.value'].notna()
df['FDR'] = np.nan
_, fdr_vals, _, _ = multipletests(df.loc[valid_mask, 'P.value'], method='fdr_bh')
df.loc[valid_mask, 'FDR'] = fdr_vals

# 2. Mark on/off proteins (any group max replicate < 0.01)
ctrl_max = df[ctrl_cols].max(axis=1)
ptx3_max = df[ptx3_cols].max(axis=1)
ONOFF_THRESH = 0.01
df['is_onoff'] = (ctrl_max < ONOFF_THRESH) | (ptx3_max < ONOFF_THRESH)

# 3. -log10P
df['-log10P'] = -np.log10(df['P.value'].clip(lower=1e-15))

# 4. Dual-tier classification
#    Stringent: p < 0.01 AND |log2FC| >= 1.0 AND not on/off  (primary DEPs)
#    Exploratory: p < 0.05 AND |log2FC| >= 1.0 AND not on/off (for enrichment)
#    On/Off: separately flagged

FC_CUTOFF = 1.0
P_STRINGENT = 0.01
P_EXPLORATORY = 0.05

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

# Convenience subsets
stringent_up   = df[df['Class'] == 'Stringent_Up']
stringent_down = df[df['Class'] == 'Stringent_Down']
exploratory_up   = df[df['Class'] == 'Exploratory_Up']
exploratory_down = df[df['Class'] == 'Exploratory_Down']
onoff_df = df[df['Class'] == 'On/Off']

# Combined sets for plotting and enrichment
all_up   = df[df['Class'].isin(['Stringent_Up', 'Exploratory_Up'])]
all_down = df[df['Class'].isin(['Stringent_Down', 'Exploratory_Down'])]
all_sig  = pd.concat([all_up, all_down])

# For enrichment: use exploratory tier (broader gene list)
enrich_up   = pd.concat([stringent_up, exploratory_up])
enrich_down = pd.concat([stringent_down, exploratory_down])

print(f"Total proteins         : {len(df)}")
print(f"Stringent Up (p<0.01)  : {len(stringent_up)}")
print(f"Stringent Down (p<0.01): {len(stringent_down)}")
print(f"Exploratory Up (p<0.05): {len(exploratory_up)}")
print(f"Exploratory Down       : {len(exploratory_down)}")
print(f"On/Off (p<0.05)        : {len(onoff_df)}")
print(f"Total DEPs (excl on/off): {len(all_sig)}")

# ══════════════════════════════════════════════════════════════════════════════
# 1. VOLCANO PLOT
# ══════════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(8, 7))

# NS
ns = df[df['Class'] == 'NS']
ax.scatter(ns['log2FC'], ns['-log10P'], s=8, c='#BBBBBB', alpha=0.5,
           linewidths=0, rasterized=True, label='NS')

# Exploratory (lighter colors)
ax.scatter(exploratory_up['log2FC'], exploratory_up['-log10P'],
           s=16, c='#E8927C', alpha=0.7, linewidths=0, label=f'Exploratory Up ({len(exploratory_up)})')
ax.scatter(exploratory_down['log2FC'], exploratory_down['-log10P'],
           s=16, c='#7FAED4', alpha=0.7, linewidths=0, label=f'Exploratory Down ({len(exploratory_down)})')

# Stringent (saturated colors)
ax.scatter(stringent_up['log2FC'], stringent_up['-log10P'],
           s=22, c='#D73027', alpha=0.9, linewidths=0, label=f'Stringent Up ({len(stringent_up)})')
ax.scatter(stringent_down['log2FC'], stringent_down['-log10P'],
           s=22, c='#4575B4', alpha=0.9, linewidths=0, label=f'Stringent Down ({len(stringent_down)})')

# On/Off (triangle markers)
onoff_up = onoff_df[onoff_df['log2FC'] > 0]
onoff_dn = onoff_df[onoff_df['log2FC'] < 0]
ax.scatter(onoff_up['log2FC'], onoff_up['-log10P'],
           s=20, c='#D73027', alpha=0.6, marker='^', linewidths=0.3, edgecolors='k',
           label=f'On/Off Up ({len(onoff_up)})')
ax.scatter(onoff_dn['log2FC'], onoff_dn['-log10P'],
           s=20, c='#4575B4', alpha=0.6, marker='^', linewidths=0.3, edgecolors='k',
           label=f'On/Off Down ({len(onoff_dn)})')

# Reference lines
ax.axhline(-np.log10(0.05), ls='--', lw=0.8, c='#444444', label='p = 0.05')
ax.axvline(0, ls='-', lw=0.5, c='#888888')
ax.axvline(FC_CUTOFF, ls=':', lw=0.8, c='#999999', label=f'|log2FC| = {FC_CUTOFF}')
ax.axvline(-FC_CUTOFF, ls=':', lw=0.8, c='#999999')
# FDR line — compute -log10(max p that passes FDR<0.05)
fdr_05_pvals = df.loc[df['FDR'] < 0.05, 'P.value']
if len(fdr_05_pvals) > 0:
    fdr_line = -np.log10(fdr_05_pvals.max())
    ax.axhline(fdr_line, ls='-.', lw=0.8, c='#E67E22', label='FDR = 0.05')

# Label top 5 up + top 5 down (quantitative DEPs)
label_up   = all_up.nlargest(5, '-log10P')
label_down = all_down.nlargest(5, '-log10P')
label_all  = pd.concat([label_up, label_down])

# Also label key On/Off proteins: PTX3 (experimental target) + top 2 on/off each direction
ptx3_row = onoff_df[onoff_df['Gene'] == 'PTX3']
onoff_label_up = onoff_up[onoff_up['Gene'] != 'PTX3'].nlargest(2, '-log10P')
onoff_label_dn = onoff_dn.nlargest(2, '-log10P')
label_onoff = pd.concat([ptx3_row, onoff_label_up, onoff_label_dn])
label_all = pd.concat([label_all, label_onoff])

texts = []
for _, row in label_all.iterrows():
    col = '#B22222' if row['log2FC'] > 0 else '#2B4D9E'
    gene_name = row['Gene']
    # Highlight PTX3 specially
    if gene_name == 'PTX3':
        t = ax.annotate('PTX3 *', (row['log2FC'], row['-log10P']),
                        fontsize=9, color='#8B0000', fontweight='bold',
                        xytext=(4, 3), textcoords='offset points',
                        ha='left', va='bottom',
                        bbox=dict(boxstyle='round,pad=0.15', fc='#FFFFCC', ec='#8B0000', alpha=0.8))
    else:
        t = ax.annotate(gene_name, (row['log2FC'], row['-log10P']),
                        fontsize=8, color=col, fontweight='bold',
                        xytext=(4, 3), textcoords='offset points',
                        ha='left', va='bottom')
    texts.append(t)

if HAS_ADJUST and texts:
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))

ax.set_xlabel('log2(Fold Change)', fontsize=13)
ax.set_ylabel('-log10(P value)', fontsize=13)
ax.set_title('PTX3 Treatment vs Control (60 h)\nVolcano Plot', fontsize=14, fontweight='bold')
ax.legend(frameon=True, fontsize=7.5, loc='upper left', ncol=1)
ax.tick_params(labelsize=11)
ax.spines[['top', 'right']].set_visible(False)

ax.text(0.98, 0.98,
        f'n = {len(df)} proteins\n'
        f'Stringent: {len(stringent_up)}↑ {len(stringent_down)}↓\n'
        f'Exploratory: {len(exploratory_up)}↑ {len(exploratory_down)}↓\n'
        f'On/Off: {len(onoff_df)}',
        transform=ax.transAxes, fontsize=8, va='top', ha='right',
        bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='gray', alpha=0.8))

plt.tight_layout()
plt.savefig(f'{OUT}/01_volcano_plot.pdf', dpi=300, bbox_inches='tight')
plt.savefig(f'{OUT}/01_volcano_plot.png', dpi=300, bbox_inches='tight')
plt.close()
print("✓ Volcano plot saved")


# ══════════════════════════════════════════════════════════════════════════════
# 2. HEATMAP – top 20 DEPs (quantitative only)
# ══════════════════════════════════════════════════════════════════════════════
n_heat_per = min(10, max(len(all_up), 1))
top_up_heat   = all_up.nlargest(n_heat_per, '-log10P')
n_heat_per_dn = min(10, max(len(all_down), 1))
top_down_heat = all_down.nlargest(n_heat_per_dn, '-log10P')
heat_df = pd.concat([top_up_heat, top_down_heat])

if len(heat_df) > 0:
    heat_df = heat_df.set_index('Gene')
    mat = heat_df[ctrl_cols + ptx3_cols].copy().astype(float)
    mat_z = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1) + 1e-9, axis=0)

    fig_h = max(5, len(heat_df) * 0.4)
    fig, ax = plt.subplots(figsize=(9, fig_h))
    sns.heatmap(mat_z, cmap='RdBu_r', center=0, vmin=-2.5, vmax=2.5,
                linewidths=0.3, linecolor='white',
                cbar_kws={'label': 'Z-score', 'shrink': 0.5},
                xticklabels=['Ctrl-1', 'Ctrl-2', 'Ctrl-3', 'PTX3-1', 'PTX3-2', 'PTX3-3'],
                ax=ax)
    ax.set_title(f'Top {len(heat_df)} Differentially Expressed Proteins\n(Z-score normalized, quantitative DEPs only)',
                 fontsize=13, fontweight='bold', pad=12)
    ax.set_ylabel('')
    ax.tick_params(axis='y', labelsize=9)
    ax.tick_params(axis='x', labelsize=10, rotation=30)

    up_genes_set = set(top_up_heat['Gene'])
    for label in ax.get_yticklabels():
        label.set_color('#B22222' if label.get_text() in up_genes_set else '#2B4D9E')

    plt.tight_layout()
    plt.savefig(f'{OUT}/02_heatmap_top_DEPs.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{OUT}/02_heatmap_top_DEPs.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("✓ Heatmap saved")
else:
    print("  No quantitative DEPs for heatmap")


# ══════════════════════════════════════════════════════════════════════════════
# 3. BAR CHART – top up + top down (with FDR values)
# ══════════════════════════════════════════════════════════════════════════════
n_bar = min(20, len(all_up))
n_bar_dn = min(20, len(all_down))
top_bar_u = all_up.nlargest(n_bar, 'log2FC').sort_values('log2FC')
top_bar_d = all_down.nsmallest(n_bar_dn, 'log2FC').sort_values('log2FC', ascending=False)

fig, axes = plt.subplots(1, 2, figsize=(14, max(5, max(n_bar, n_bar_dn) * 0.35)))

axes[0].barh(top_bar_u['Gene'], top_bar_u['log2FC'], color='#D73027', edgecolor='none', height=0.7)
axes[0].set_xlabel('log2(Fold Change)', fontsize=12)
axes[0].set_title(f'Top {n_bar} Up-regulated', fontsize=12, fontweight='bold', color='#D73027')
axes[0].axvline(0, color='black', lw=0.8)
axes[0].spines[['top', 'right']].set_visible(False)
axes[0].tick_params(labelsize=9)
for i, (_, row) in enumerate(top_bar_u.iterrows()):
    fdr_str = f"FDR={row['FDR']:.3f}" if pd.notna(row['FDR']) and row['FDR'] < 1 else f"p={row['P.value']:.3f}"
    axes[0].text(row['log2FC'] + 0.1, i, fdr_str, va='center', fontsize=8, color='#555555')

axes[1].barh(top_bar_d['Gene'], top_bar_d['log2FC'], color='#4575B4', edgecolor='none', height=0.7)
axes[1].set_xlabel('log2(Fold Change)', fontsize=12)
axes[1].set_title(f'Top {n_bar_dn} Down-regulated', fontsize=12, fontweight='bold', color='#4575B4')
axes[1].axvline(0, color='black', lw=0.8)
axes[1].spines[['top', 'right']].set_visible(False)
axes[1].tick_params(labelsize=9)
for i, (_, row) in enumerate(top_bar_d.iterrows()):
    fdr_str = f"FDR={row['FDR']:.3f}" if pd.notna(row['FDR']) and row['FDR'] < 1 else f"p={row['P.value']:.3f}"
    axes[1].text(row['log2FC'] - 0.1, i, fdr_str, va='center', ha='right', fontsize=8, color='#555555')

plt.suptitle('Top Differentially Expressed Proteins – PTX3 vs Control (60 h)\n'
             '(Quantitative DEPs only, excluding On/Off artifacts)',
             fontsize=13, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(f'{OUT}/03_bar_chart_top_DEPs.pdf', dpi=300, bbox_inches='tight')
plt.savefig(f'{OUT}/03_bar_chart_top_DEPs.png', dpi=300, bbox_inches='tight')
plt.close()
print("✓ Bar chart saved")


# ══════════════════════════════════════════════════════════════════════════════
# 4. SAMPLE CORRELATION HEATMAP
# ══════════════════════════════════════════════════════════════════════════════
all_sample_cols = ctrl_cols + ptx3_cols
sample_mat = df[all_sample_cols].dropna().apply(np.log2).replace([np.inf, -np.inf], np.nan).dropna()
corr_mat = sample_mat.corr(method='pearson')
labels = ['Ctrl-1', 'Ctrl-2', 'Ctrl-3', 'PTX3-1', 'PTX3-2', 'PTX3-3']

fig, ax = plt.subplots(figsize=(7, 6))
mask = np.triu(np.ones_like(corr_mat, dtype=bool), k=1)
sns.heatmap(corr_mat, mask=mask, annot=True, fmt='.4f',
            cmap='YlOrRd', vmin=0.95, vmax=1.0,
            xticklabels=labels, yticklabels=labels,
            linewidths=0.5, ax=ax,
            cbar_kws={'label': 'Pearson r', 'shrink': 0.7})
ax.set_title('Sample-to-Sample Pearson Correlation\n(log₂ protein intensities)',
             fontsize=12, fontweight='bold')
ax.tick_params(labelsize=10)
plt.tight_layout()
plt.savefig(f'{OUT}/04_sample_correlation.pdf', dpi=300, bbox_inches='tight')
plt.savefig(f'{OUT}/04_sample_correlation.png', dpi=300, bbox_inches='tight')
plt.close()
print("✓ Correlation heatmap saved")


# ══════════════════════════════════════════════════════════════════════════════
# 5. P-VALUE DISTRIBUTION + FC DISTRIBUTION
# ══════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

axes[0].hist(df['P.value'].dropna(), bins=50, color='#5B8DB8', edgecolor='white', lw=0.3)
axes[0].set_xlabel('P value', fontsize=12)
axes[0].set_ylabel('Number of proteins', fontsize=12)
axes[0].set_title('P-value Distribution', fontsize=12, fontweight='bold')
axes[0].axvline(0.05, ls='--', color='red', lw=1.2, label='p = 0.05')
axes[0].axvline(0.01, ls='--', color='orange', lw=1.0, label='p = 0.01')
axes[0].legend(fontsize=10)
axes[0].spines[['top', 'right']].set_visible(False)

axes[1].hist(df['log2FC'].clip(-6, 6), bins=60, color='#7DBD8E', edgecolor='white', lw=0.3)
axes[1].axvline(0, color='black', lw=0.8)
axes[1].axvline(FC_CUTOFF, ls=':', color='red', lw=1.0, label=f'|log2FC| = {FC_CUTOFF}')
axes[1].axvline(-FC_CUTOFF, ls=':', color='red', lw=1.0)
axes[1].set_xlabel('log2(Fold Change)', fontsize=12)
axes[1].set_ylabel('Number of proteins', fontsize=12)
axes[1].set_title('log2(FC) Distribution', fontsize=12, fontweight='bold')
axes[1].legend(fontsize=10)
axes[1].spines[['top', 'right']].set_visible(False)

plt.suptitle('Global Proteome Statistics – PTX3 vs Control', fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig(f'{OUT}/05_distribution_plots.pdf', dpi=300, bbox_inches='tight')
plt.savefig(f'{OUT}/05_distribution_plots.png', dpi=300, bbox_inches='tight')
plt.close()
print("✓ Distribution plots saved")


# ══════════════════════════════════════════════════════════════════════════════
# 6. SCATTER PLOT – mean CTRL vs mean PTX3
# ══════════════════════════════════════════════════════════════════════════════
df_scatter = df[df['mean c'] > 0.01].copy()
log_c = np.log2(df_scatter['mean c'])
log_p = np.log2(df_scatter['mean ptx3'].clip(lower=0.001))

fig, ax = plt.subplots(figsize=(7, 7))
color_map = {'NS': '#CCCCCC',
             'Stringent_Up': '#D73027', 'Exploratory_Up': '#E8927C',
             'Stringent_Down': '#4575B4', 'Exploratory_Down': '#7FAED4',
             'On/Off': '#888888'}
colors = df_scatter['Class'].map(color_map).fillna('#CCCCCC')
ax.scatter(log_c, log_p, c=colors, s=8, alpha=0.6, linewidths=0)

lim = [min(log_c.min(), log_p.min()) - 0.5, max(log_c.max(), log_p.max()) + 0.5]
ax.plot(lim, lim, 'k--', lw=1, label='y = x')
ax.set_xlim(lim)
ax.set_ylim(lim)

r, p_r = pearsonr(log_c, log_p)
ax.set_xlabel('log2(mean Intensity) - Control', fontsize=12)
ax.set_ylabel('log2(mean Intensity) - PTX3', fontsize=12)
ax.set_title(f'Protein Intensity Correlation\n(r = {r:.4f})', fontsize=12, fontweight='bold')
patches = [mpatches.Patch(color='#D73027', label=f'Stringent Up ({len(stringent_up)})'),
           mpatches.Patch(color='#4575B4', label=f'Stringent Down ({len(stringent_down)})'),
           mpatches.Patch(color='#E8927C', label=f'Exploratory Up ({len(exploratory_up)})'),
           mpatches.Patch(color='#7FAED4', label=f'Exploratory Down ({len(exploratory_down)})'),
           mpatches.Patch(color='#CCCCCC', label='NS')]
ax.legend(handles=patches, fontsize=8)
ax.spines[['top', 'right']].set_visible(False)
plt.tight_layout()
plt.savefig(f'{OUT}/06_scatter_ctrl_vs_ptx3.pdf', dpi=300, bbox_inches='tight')
plt.savefig(f'{OUT}/06_scatter_ctrl_vs_ptx3.png', dpi=300, bbox_inches='tight')
plt.close()
print("✓ Scatter plot saved")


# ══════════════════════════════════════════════════════════════════════════════
# 7. ENRICHR PATHWAY ENRICHMENT (online API)
# ══════════════════════════════════════════════════════════════════════════════
def enrichr_query(gene_list, lib='KEGG_2021_Human', name='gene_set', min_genes=3):
    """Query Enrichr API and filter results by minimum gene overlap."""
    if len(gene_list) < 2:
        print(f"  Skipping {lib} — only {len(gene_list)} genes")
        return pd.DataFrame()
    try:
        url_add = 'https://maayanlab.cloud/Enrichr/addList'
        gene_str = '\n'.join(gene_list)
        boundary = '----PythonFormBoundary'
        body = (
            f'--{boundary}\r\n'
            f'Content-Disposition: form-data; name="list"\r\n\r\n'
            f'{gene_str}\r\n'
            f'--{boundary}\r\n'
            f'Content-Disposition: form-data; name="description"\r\n\r\n'
            f'{name}\r\n'
            f'--{boundary}--\r\n'
        ).encode('utf-8')
        req = urllib.request.Request(url_add, data=body,
            headers={'Content-Type': f'multipart/form-data; boundary={boundary}'})
        resp = json.loads(urllib.request.urlopen(req, timeout=15).read())
        user_list_id = resp['userListId']

        url_enr = (f'https://maayanlab.cloud/Enrichr/enrich?'
                   f'userListId={user_list_id}&backgroundType={lib}')
        result = json.loads(urllib.request.urlopen(url_enr, timeout=15).read())
        results = result[lib]
        out = []
        for r in results[:30]:
            genes = r[5]
            n_genes = len(genes)
            if n_genes < min_genes:
                continue
            out.append({
                'Rank': r[0], 'Term': r[1], 'P.value': r[2],
                'Adj.P': r[6], 'Combined_Score': r[4],
                'n_genes': n_genes,
                'Genes': ';'.join(genes)
            })
        return pd.DataFrame(out)
    except Exception as e:
        print(f"  Enrichr API error ({lib}): {e}")
        return pd.DataFrame()

# Gene lists for enrichment (quantitative DEPs only, excluding on/off)
up_genes_list   = enrich_up['Gene'].dropna().tolist()
down_genes_list = enrich_down['Gene'].dropna().tolist()
print(f"\nEnrichment gene lists: {len(up_genes_list)} up, {len(down_genes_list)} down")

# Query all libraries
print("Querying Enrichr (KEGG) for up-regulated genes...")
kegg_up = enrichr_query(up_genes_list, 'KEGG_2021_Human', 'PTX3_up')
print("Querying Enrichr (KEGG) for down-regulated genes...")
kegg_dn = enrichr_query(down_genes_list, 'KEGG_2021_Human', 'PTX3_down')

print("Querying Enrichr (GO BP) for up-regulated genes...")
gobp_up = enrichr_query(up_genes_list, 'GO_Biological_Process_2021', 'PTX3_up_GO')
print("Querying Enrichr (GO BP) for down-regulated genes...")
gobp_dn = enrichr_query(down_genes_list, 'GO_Biological_Process_2021', 'PTX3_down_GO')

print("Querying Enrichr (GO MF) ...")
gomf_up = enrichr_query(up_genes_list, 'GO_Molecular_Function_2021', 'PTX3_up_MF')

print("Querying Enrichr (GO CC) ...")
gocc_up = enrichr_query(up_genes_list, 'GO_Cellular_Component_2021', 'PTX3_up_CC')

print("Querying Enrichr (WikiPathways) ...")
wiki_up = enrichr_query(up_genes_list, 'WikiPathway_2021_Human', 'PTX3_up_Wiki')


# ── Plot enrichment bar charts ───────────────────────────────────────────────
def plot_enrichment_bar(enrich_df, title, color, filename, top_n=15):
    """Plot enrichment using -log10(Adj.P) on x-axis with significance line."""
    if enrich_df.empty:
        print(f"  No enrichment data for: {title}")
        return
    df_plot = enrich_df.head(top_n).sort_values('Adj.P', ascending=True)
    df_plot = df_plot.iloc[::-1]  # reverse for barh (bottom = most significant)
    df_plot = df_plot.copy()
    df_plot['Term_short'] = df_plot['Term'].str.split(' Homo sapiens').str[0].str[:55]
    df_plot['-log10AdjP'] = -np.log10(df_plot['Adj.P'].clip(lower=1e-10))

    fig, ax = plt.subplots(figsize=(10, max(4, len(df_plot) * 0.5)))
    bars = ax.barh(df_plot['Term_short'], df_plot['-log10AdjP'],
                   color=color, edgecolor='none', height=0.65)

    # Fade bars that don't pass Adj.P < 0.05
    for bar, adjp in zip(bars, df_plot['Adj.P']):
        if adjp >= 0.05:
            bar.set_alpha(0.4)

    # Significance threshold line
    ax.axvline(-np.log10(0.05), ls='--', lw=1.0, c='#333333', label='Adj.P = 0.05')

    ax.set_xlabel('-log10(Adjusted P value)', fontsize=12)
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(labelsize=9)
    ax.legend(fontsize=9)

    for i, (_, row) in enumerate(df_plot.iterrows()):
        label_text = f"n={row['n_genes']}  adj.p={row['Adj.P']:.2e}"
        ax.text(0.3, i, label_text, va='center', fontsize=7, color='white', fontweight='bold')

    plt.tight_layout()
    plt.savefig(f'{OUT}/{filename}.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{OUT}/{filename}.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✓ {filename} saved")

plot_enrichment_bar(kegg_up, 'KEGG Pathway Enrichment – Up-regulated Proteins\n(PTX3 vs Control, 60 h)',
                    '#C0392B', '07_KEGG_up')
plot_enrichment_bar(kegg_dn, 'KEGG Pathway Enrichment – Down-regulated Proteins\n(PTX3 vs Control, 60 h)',
                    '#2980B9', '08_KEGG_down')
plot_enrichment_bar(gobp_up, 'GO Biological Process – Up-regulated Proteins\n(PTX3 vs Control, 60 h)',
                    '#E67E22', '09_GOBP_up', top_n=20)
plot_enrichment_bar(gobp_dn, 'GO Biological Process – Down-regulated Proteins\n(PTX3 vs Control, 60 h)',
                    '#8E44AD', '10_GOBP_down', top_n=20)
plot_enrichment_bar(gomf_up, 'GO Molecular Function – Up-regulated Proteins\n(PTX3 vs Control, 60 h)',
                    '#16A085', '11_GOMF_up', top_n=15)
plot_enrichment_bar(gocc_up, 'GO Cellular Component – Up-regulated Proteins\n(PTX3 vs Control, 60 h)',
                    '#2E86C1', '12_GOCC_up', top_n=15)


# ══════════════════════════════════════════════════════════════════════════════
# 8. SAVE FULL DEP TABLES
# ══════════════════════════════════════════════════════════════════════════════
export_cols = ['Protein', 'Protein ID', 'Gene', 'log2FC', 'FC', 'P.value',
               'FDR', 'mean c', 'mean ptx3', 'Class', 'is_onoff']

# All quantitative DEPs
sig_export = all_sig[export_cols].sort_values('log2FC', ascending=False)
sig_export.to_csv(f'{OUT}/DEPs_table.csv', index=False, encoding='utf-8-sig')
print(f"✓ DEPs table saved ({len(sig_export)} proteins)")

# On/off proteins table
onoff_export = onoff_df[export_cols].sort_values('log2FC', ascending=False)
onoff_export.to_csv(f'{OUT}/OnOff_proteins_table.csv', index=False, encoding='utf-8-sig')
print(f"✓ On/Off table saved ({len(onoff_export)} proteins)")

# Full table with all classifications
df[export_cols].to_csv(f'{OUT}/full_protein_table.csv', index=False, encoding='utf-8-sig')
print("✓ Full protein table saved")

# Save enrichment results
for name, dff in [('KEGG_up', kegg_up), ('KEGG_down', kegg_dn),
                  ('GOBP_up', gobp_up), ('GOBP_down', gobp_dn),
                  ('GOMF_up', gomf_up), ('GOCC_up', gocc_up),
                  ('WikiPath_up', wiki_up)]:
    if not dff.empty:
        dff.to_csv(f'{OUT}/enrichment_{name}.csv', index=False, encoding='utf-8-sig')
print("✓ Enrichment tables saved")


# ══════════════════════════════════════════════════════════════════════════════
# 9. COMPREHENSIVE TEXT REPORT
# ══════════════════════════════════════════════════════════════════════════════
report_lines = []
R = report_lines.append

R("=" * 72)
R("PTX3 PROTEOMICS ANALYSIS – COMPREHENSIVE REPORT (v2)")
R("Experiment: PTX3 Treatment vs Control (60 h), 3+3 Replicates")
R("=" * 72)
R("")

R("─" * 60)
R("1. DATA OVERVIEW")
R("─" * 60)
R(f"  Total proteins quantified  : {len(df):,}")
R(f"  log2FC range               : {df['log2FC'].min():.2f} to {df['log2FC'].max():.2f}")
R(f"  Median log2FC              : {df['log2FC'].median():.4f}")
R("")

R("─" * 60)
R("2. DEP CLASSIFICATION (dual-tier)")
R("─" * 60)
R(f"  FC cutoff: |log2FC| >= {FC_CUTOFF}")
R(f"  On/Off threshold: max replicate intensity < {ONOFF_THRESH}")
R(f"  FDR method: Benjamini-Hochberg")
R("")
R(f"  STRINGENT (p < {P_STRINGENT}, |log2FC| >= {FC_CUTOFF}, excl on/off):")
R(f"    Up-regulated   : {len(stringent_up)}")
R(f"    Down-regulated : {len(stringent_down)}")
R(f"  EXPLORATORY (p < {P_EXPLORATORY}, |log2FC| >= {FC_CUTOFF}, excl on/off):")
R(f"    Up-regulated   : {len(exploratory_up)}")
R(f"    Down-regulated : {len(exploratory_down)}")
R(f"  TOTAL quantitative DEPs: {len(all_sig)} ({len(all_up)} up, {len(all_down)} down)")
R(f"  ON/OFF proteins (p < 0.05): {len(onoff_df)}")
R("")
R("  NOTE: With n=3 per group, BH FDR correction is very conservative.")
R(f"  Only {(df['FDR'] < 0.05).sum()} proteins pass FDR < 0.05 (all on/off artifacts).")
R("  The dual-tier approach uses raw p-value thresholds with FC cutoff,")
R("  which is appropriate for discovery-phase proteomics with small n.")
R("")

# Sample correlation
log_mat = df[ctrl_cols + ptx3_cols].clip(lower=0.001).apply(np.log2)
corr = log_mat.corr()
ctrl_intra = [corr.loc[a, b] for i, a in enumerate(ctrl_cols)
              for j, b in enumerate(ctrl_cols) if i < j]
ptx3_intra = [corr.loc[a, b] for i, a in enumerate(ptx3_cols)
              for j, b in enumerate(ptx3_cols) if i < j]

R("─" * 60)
R("3. SAMPLE QUALITY CONTROL")
R("─" * 60)
R(f"  Intra-Control Pearson r    : {np.mean(ctrl_intra):.4f} ± {np.std(ctrl_intra):.4f}")
R(f"  Intra-PTX3 Pearson r       : {np.mean(ptx3_intra):.4f} ± {np.std(ptx3_intra):.4f}")
R("")

R("─" * 60)
R("4. TOP STRINGENT DEPs (p < 0.01, |log2FC| >= 1.0, excl on/off)")
R("─" * 60)
R(f"  {'Gene':<14} {'log2FC':>8}  {'P value':>10}  {'FDR':>10}  {'FC':>8}")
R(f"  {'-'*60}")

top_str_up = stringent_up.nlargest(15, 'log2FC')
for _, row in top_str_up.iterrows():
    fdr_s = f"{row['FDR']:.4f}" if pd.notna(row['FDR']) else "NA"
    R(f"  {row['Gene']:<14} {row['log2FC']:>8.3f}  {row['P.value']:>10.4f}  {fdr_s:>10}  {row['FC']:>8.2f}")
R("")
top_str_dn = stringent_down.nsmallest(15, 'log2FC')
for _, row in top_str_dn.iterrows():
    fdr_s = f"{row['FDR']:.4f}" if pd.notna(row['FDR']) else "NA"
    R(f"  {row['Gene']:<14} {row['log2FC']:>8.3f}  {row['P.value']:>10.4f}  {fdr_s:>10}  {row['FC']:>8.2f}")
R("")

R("─" * 60)
R("5. ON/OFF PROTEINS (p < 0.05, one group near-zero)")
R("─" * 60)
R("  These proteins are detected in one condition but not the other.")
R("  Their extreme log2FC values are artifacts of near-zero denominators.")
R("  They are excluded from enrichment analysis but may be biologically")
R("  relevant if validated by orthogonal methods.")
R("")
R(f"  {'Gene':<14} {'log2FC':>8}  {'P value':>10}  {'Direction':<10}")
R(f"  {'-'*50}")
for _, row in onoff_df.sort_values('log2FC', ascending=False).head(15).iterrows():
    direction = 'Gained' if row['log2FC'] > 0 else 'Lost'
    R(f"  {row['Gene']:<14} {row['log2FC']:>8.2f}  {row['P.value']:>10.4f}  {direction:<10}")
R("")

R("─" * 60)
R("6. PATHWAY ENRICHMENT SUMMARY")
R("─" * 60)
R(f"  Gene lists: {len(up_genes_list)} up-regulated, {len(down_genes_list)} down-regulated")
R(f"  Minimum gene overlap filter: n >= 3")
R("")
for label, edf in [('KEGG – Up', kegg_up), ('KEGG – Down', kegg_dn),
                   ('GO BP – Up', gobp_up), ('GO BP – Down', gobp_dn),
                   ('GO MF – Up', gomf_up), ('GO CC – Up', gocc_up),
                   ('WikiPathways – Up', wiki_up)]:
    R(f"\n  [{label}]")
    if edf.empty:
        R("  (No terms with n >= 3 genes)")
    else:
        n_sig = (edf['Adj.P'] < 0.05).sum()
        R(f"  ({len(edf)} terms total, {n_sig} with Adj.P < 0.05)")
        for _, row in edf.head(8).iterrows():
            sig_mark = "**" if row['Adj.P'] < 0.05 else "  "
            R(f"  {sig_mark}{row['Rank']:>2}. {row['Term'][:50]:<50} adj.p={row['Adj.P']:.2e} n={row['n_genes']}")
R("")

R("─" * 60)
R("7. STATISTICAL NOTES")
R("─" * 60)
R(f"""
  - Quantification: label-free protein intensity (normalized)
  - DEP classification uses dual-tier thresholds:
    Stringent : raw p < {P_STRINGENT} AND |log2FC| >= {FC_CUTOFF} (excl on/off)
    Exploratory: raw p < {P_EXPLORATORY} AND |log2FC| >= {FC_CUTOFF} (excl on/off)
  - BH FDR correction computed but not used as primary threshold
    (too conservative with n=3; only {(df['FDR'] < 0.05).sum()} proteins pass FDR < 0.05)
  - On/Off proteins (max replicate < {ONOFF_THRESH} in either group)
    are flagged and excluded from enrichment analysis
  - Enrichment terms filtered to require >= 3 overlapping genes
  - Enrichment significance shown as Adj.P (Enrichr BH-corrected)
""")

R("=" * 72)
R("Files generated in: " + OUT)
R("=" * 72)

report_text = "\n".join(report_lines)
with open(f'{OUT}/ANALYSIS_REPORT.txt', 'w', encoding='utf-8') as f:
    f.write(report_text)
print("\n" + report_text)
print(f"\n✓ All outputs saved to: {OUT}")
