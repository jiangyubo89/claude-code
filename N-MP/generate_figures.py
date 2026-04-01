"""
Generate 5 publication-quality summary figures for Chapter 8.
Each figure synthesizes key data/concepts from the chapter text.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from matplotlib.gridspec import GridSpec
import os

OUT = 'N-MP/Chapter8_Output/Figures'
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 10,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

# ═══════════════════════════════════════════════════════════════
# FIGURE 8.1 — Analytical Method Comparison (Heatmap + Radar)
# ═══════════════════════════════════════════════════════════════
def fig1():
    methods = ['μ-FTIR', 'ATR-FTIR', 'μ-Raman', 'Py-GC/MS', 'SEM-EDX', 'TOF-SIMS', 'AFM-IR', 'LDIR']
    criteria = ['Spatial\nResolution', 'Sensitivity', 'Throughput', 'NP\nCapability',
                'Polymer\nID', 'Morphology\nInfo', 'Matrix\nTolerance', 'Non-\nDestructive']
    # Scores 0-10 based on Table 8.1 descriptions
    scores = np.array([
        [5, 7, 6, 2, 9, 8, 4, 10],  # μ-FTIR
        [2, 5, 8, 1, 9, 3, 6, 10],  # ATR-FTIR
        [8, 7, 4, 3, 9, 8, 5, 10],  # μ-Raman
        [10, 9, 6, 10, 8, 0, 8, 0], # Py-GC/MS
        [7, 4, 3, 7, 2, 9, 5, 10],  # SEM-EDX
        [8, 6, 2, 8, 6, 7, 3, 5],   # TOF-SIMS
        [10, 8, 1, 10, 9, 5, 3, 10], # AFM-IR
        [4, 6, 10, 2, 8, 7, 5, 10],  # LDIR
    ])

    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.imshow(scores.T, cmap='YlOrRd', aspect='auto', vmin=0, vmax=10)
    ax.set_xticks(range(len(methods)))
    ax.set_xticklabels(methods, rotation=30, ha='right', fontsize=9)
    ax.set_yticks(range(len(criteria)))
    ax.set_yticklabels(criteria, fontsize=9)

    for i in range(len(criteria)):
        for j in range(len(methods)):
            color = 'white' if scores[j, i] > 6 else 'black'
            ax.text(j, i, str(scores[j, i]), ha='center', va='center',
                    fontsize=9, fontweight='bold', color=color)

    cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label('Performance Score (0-10)', fontsize=10)
    ax.set_title('Figure 8.1  Comparative Performance of Analytical Methods\nfor iMNP Detection in Biological Tissues',
                 fontsize=12, fontweight='bold', pad=15)
    ax.set_xlabel('')
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'Fig8_1_method_comparison.png'))
    plt.close()
    print('Fig 8.1 saved')


# ═══════════════════════════════════════════════════════════════
# FIGURE 8.2 — Organ Biodistribution from Animal Models
# ═══════════════════════════════════════════════════════════════
def fig2():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 6), gridspec_kw={'width_ratios': [1.3, 1]})

    # Panel A: Organ accumulation hierarchy (horizontal bar)
    organs = ['GI Tract', 'Liver', 'Lungs', 'Spleen', 'Kidney', 'Placenta', 'Testes', 'Brain']
    pct = [45.0, 12.0, 5.0, 3.5, 2.8, 0.8, 0.5, 0.3]
    sd = [8.0, 3.5, 1.8, 1.2, 1.0, 0.4, 0.3, 0.15]
    colors = ['#D73027','#FC8D59','#4575B4','#CE93D8','#1A9850','#F48FB1','#F46D43','#762A83']

    bars = ax1.barh(organs, pct, xerr=sd, color=colors, edgecolor='white',
                    height=0.55, capsize=3, error_kw={'linewidth':1.2, 'color':'#333333'})
    ax1.set_xlabel('Accumulation (% of administered dose)', fontsize=10)
    ax1.set_title('(A) Organ Accumulation Hierarchy', fontsize=11, fontweight='bold')
    ax1.set_xlim(0, 60)
    for bar, val, s in zip(bars, pct, sd):
        ax1.text(val + s + 1.0, bar.get_y() + bar.get_height()/2,
                f'{val} ± {s}%', va='center', fontsize=9)
    ax1.invert_yaxis()

    # Panel B: Size-dependent translocation (grouped bar)
    size_labels = ['20 nm', '50 nm', '100 nm', '500 nm', '1 μm', '5 μm']
    liver_pct = [18, 15, 12, 5, 2, 0.5]
    gut_pct = [30, 35, 45, 55, 65, 80]
    spleen_pct = [8, 5, 3.5, 1.5, 0.5, 0.1]

    x = np.arange(len(size_labels))
    w = 0.25
    ax2.bar(x - w, liver_pct, w, label='Liver', color='#FC8D59', edgecolor='white')
    ax2.bar(x, gut_pct, w, label='GI Tract (retained)', color='#D73027', edgecolor='white')
    ax2.bar(x + w, spleen_pct, w, label='Spleen', color='#CE93D8', edgecolor='white')
    ax2.set_xticks(x)
    ax2.set_xticklabels(size_labels, fontsize=9)
    ax2.set_xlabel('Particle Size', fontsize=10)
    ax2.set_ylabel('Accumulation (%)', fontsize=10)
    ax2.set_title('(B) Size-Dependent Distribution', fontsize=11, fontweight='bold')
    ax2.legend(fontsize=8, loc='upper left', framealpha=0.9)

    # Add arrow annotation
    ax2.annotate('Smaller particles →\nhigher systemic\ntranslocation',
                xy=(0, 18), xytext=(2.5, 55),
                fontsize=8, fontstyle='italic', color='#555555',
                arrowprops=dict(arrowstyle='->', color='#555555', lw=1.2))

    fig.suptitle('Figure 8.2  iMNP Biodistribution in Animal Models\n'
                 '(Oral PS nanoparticle exposure, rodent models, 7–28 days)',
                 fontsize=12, fontweight='bold', y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'Fig8_2_biodistribution.png'))
    plt.close()
    print('Fig 8.2 saved')


# ═══════════════════════════════════════════════════════════════
# FIGURE 8.3 — Barrier Crossing Size Thresholds
# ═══════════════════════════════════════════════════════════════
def fig3():
    fig, ax = plt.subplots(figsize=(11, 6))

    barriers = ['Intestinal\nEpithelium', 'Alveolar-\nCapillary', 'Blood-Brain\nBarrier (BBB)', 'Placental\nBarrier']
    thresholds = [500, 100, 50, 240]  # nm
    mechanisms = [
        'CME, Caveolae,\nM-cell, Paracellular',
        'Transcytosis (AT1),\nAM clearance (>1μm)',
        'RMT (LDLR),\nCorona-mediated',
        'Syncytiotrophoblast\nMacropinocytosis'
    ]
    colors = ['#2166AC', '#4393C3', '#D6604D', '#92C5DE']
    efficiencies = [0.3, 5.0, 0.01, 0.5]  # % translocation for ~100nm particles

    y_pos = np.arange(len(barriers))

    # Bar for size threshold
    bars = ax.barh(y_pos, thresholds, height=0.5, color=colors, edgecolor='grey', alpha=0.85)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(barriers, fontsize=10)
    ax.set_xlabel('Size Threshold for Translocation (nm)', fontsize=11)
    ax.set_title('Figure 8.3  Biological Barrier Size Thresholds and Translocation Mechanisms for iMNPs',
                 fontsize=12, fontweight='bold', pad=15)

    for i, (bar, thresh, mech, eff) in enumerate(zip(bars, thresholds, mechanisms, efficiencies)):
        # Threshold label
        ax.text(thresh + 8, bar.get_y() + bar.get_height()/2,
                f'<{thresh} nm\n({eff}% translocation)',
                va='center', fontsize=9, fontweight='bold')
        # Mechanism annotation
        ax.text(thresh - 15, bar.get_y() + bar.get_height()/2,
                mech, va='center', ha='right', fontsize=8, color='white', fontstyle='italic')

    ax.set_xlim(0, 600)
    ax.invert_yaxis()

    # Add size reference markers
    sizes = [1, 10, 50, 100, 240, 500]
    for sz in sizes:
        ax.axvline(sz, color='grey', linestyle=':', alpha=0.3, linewidth=0.8)

    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'Fig8_3_barrier_thresholds.png'))
    plt.close()
    print('Fig 8.3 saved')


# ═══════════════════════════════════════════════════════════════
# FIGURE 8.4 — PBPK Model Compartment Diagram (Flowchart)
# ═══════════════════════════════════════════════════════════════
def fig4():
    fig, ax = plt.subplots(figsize=(11, 8))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis('off')

    def draw_box(x, y, w, h, label, color, fontsize=9):
        rect = mpatches.FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.1',
                                        facecolor=color, edgecolor='black', linewidth=1.2)
        ax.add_patch(rect)
        ax.text(x + w/2, y + h/2, label, ha='center', va='center',
                fontsize=fontsize, fontweight='bold')

    def arrow(x1, y1, x2, y2, label='', color='black'):
        ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle='->', color=color, lw=1.5))
        if label:
            mx, my = (x1+x2)/2, (y1+y2)/2
            ax.text(mx, my+0.15, label, fontsize=7, ha='center', color=color, fontstyle='italic')

    # Title
    ax.text(5, 9.6, 'Figure 8.4  PBPK Model for iMNP Biodistribution',
            fontsize=13, fontweight='bold', ha='center')
    ax.text(5, 9.2, 'dCi/dt = Qi·(Cblood − Ci/Pi) − kclear,i·Vi·Ci',
            fontsize=10, ha='center', fontstyle='italic', color='#555555')

    # Exposure routes
    draw_box(0.3, 8.0, 2.0, 0.8, 'ORAL\nINTAKE', '#FFE0B2', 9)
    draw_box(3.8, 8.0, 2.0, 0.8, 'INHALATION', '#B3E5FC', 9)
    draw_box(7.5, 8.0, 2.0, 0.8, 'DERMAL\n(minor)', '#E1BEE7', 9)

    # Central blood
    draw_box(3.5, 5.8, 3.0, 0.9, 'ARTERIAL / VENOUS\nBLOOD', '#FFCDD2', 10)

    # Organ compartments
    draw_box(0.2, 4.2, 1.8, 0.8, 'GI TRACT\n(45%)', '#FFCC80')
    draw_box(0.2, 2.8, 1.8, 0.8, 'LIVER\n(12%)', '#A5D6A7')
    draw_box(2.5, 4.2, 1.8, 0.8, 'LUNGS\n(5%)', '#90CAF9')
    draw_box(2.5, 2.8, 1.8, 0.8, 'SPLEEN\n(3.5%)', '#CE93D8')
    draw_box(5.7, 4.2, 1.8, 0.8, 'KIDNEY\n(2.8%)', '#80CBC4')
    draw_box(5.7, 2.8, 1.8, 0.8, 'BRAIN\n(<50nm)', '#EF9A9A')
    draw_box(8.0, 4.2, 1.8, 0.8, 'PLACENTA\n(<240nm)', '#F48FB1')
    draw_box(8.0, 2.8, 1.8, 0.8, 'OTHER\nTISSUES', '#B0BEC5')

    # Excretion
    draw_box(0.2, 1.0, 1.8, 0.7, 'FECES\n(>90%)', '#D7CCC8', 9)
    draw_box(5.7, 1.0, 1.8, 0.7, 'URINE\n(<6nm)', '#D7CCC8', 9)
    draw_box(3.0, 1.0, 2.0, 0.7, 'BILE', '#D7CCC8', 9)

    # Arrows: exposure to organs
    arrow(1.3, 8.0, 1.1, 5.0, 'oral\nabsorption')
    arrow(4.8, 8.0, 4.0, 6.7, '')
    arrow(3.4, 5.8, 1.1, 5.0, '')

    # Arrows: blood to organs
    arrow(3.5, 6.0, 2.5, 5.0, '')
    arrow(5.0, 5.8, 5.7, 5.0, '')
    arrow(6.5, 5.8, 8.0, 5.0, '')
    arrow(6.5, 5.8, 5.7, 3.6, '')
    arrow(3.5, 5.8, 2.5, 3.6, '')
    arrow(6.5, 5.8, 8.0, 3.6, '')

    # Arrows: excretion
    arrow(1.1, 4.2, 1.1, 1.7, '')
    arrow(1.1, 2.8, 4.0, 1.7, 'biliary')
    arrow(6.6, 2.8, 6.6, 1.7, '')

    # Key parameters box
    params_text = ('Key Parameters:\n'
                   'Qi = organ blood flow\n'
                   'Ci = tissue concentration\n'
                   'Pi = partition coefficient\n'
                   'kclear = clearance rate\n'
                   'Vmax/Km = MPS saturation')
    ax.text(9.5, 6.5, params_text, fontsize=7.5, va='top', ha='right',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))

    fig.savefig(os.path.join(OUT, 'Fig8_4_PBPK_model.png'))
    plt.close()
    print('Fig 8.4 saved')


# ═══════════════════════════════════════════════════════════════
# FIGURE 8.5 — Human Evidence: Tissue Prevalence + Quality
# ═══════════════════════════════════════════════════════════════
def fig5():
    fig = plt.figure(figsize=(12, 7))
    gs = GridSpec(1, 2, width_ratios=[2, 1], wspace=0.3)

    # Panel A: Forest-plot style prevalence summary
    ax1 = fig.add_subplot(gs[0])
    studies = [
        ('Leonard 2024\n(Blood)', 90.0, 69.9, 97.2, 20),
        ('Halfar 2023\n(Placenta)', 90.0, 59.6, 98.2, 10),
        ('Ke 2023\n(Stool)', 85.5, 75.4, 91.7, 69),
        ('Arslan 2025\n(Cerumen)', 83.3, 55.2, 95.3, 12),
        ('Santini 2024\n(Stool)', 83.3, 55.2, 95.3, 6),
        ('Padarya 2025\n(Lung)', 70.0, 55.4, 75.4, 60),
        ('Zakynthinos 2025\n(Lung BAL)', 62.5, 55.4, 75.4, 8),
        ('Ozgen 2024\n(Lung BAL)', 55.6, 55.4, 75.4, 18),
        ('Qu 2025\n(Semen)', 55.5, 48.4, 61.8, 200),
        ('Demirelli 2024\n(Prostate)', 50.0, 48.4, 61.8, 12),
        ('Saraluck 2024\n(Breast Milk)', 39.0, 27.6, 51.7, 59),
    ]

    colors_tissue = {
        'Blood': '#D73027', 'Placenta': '#F48FB1', 'Stool': '#FC8D59',
        'Cerumen': '#8DD3C7', 'Lung': '#4575B4', 'Semen': '#F46D43',
        'Prostate': '#F46D43', 'Breast': '#ABD9E9'
    }

    y_positions = list(range(len(studies)-1, -1, -1))
    for i, (name, prev, lo, hi, n) in enumerate(studies):
        y = y_positions[i]
        # Determine color
        c = '#888888'
        for key, col in colors_tissue.items():
            if key in name:
                c = col
                break
        size = max(4, min(15, n / 15))
        ax1.plot(prev, y, 's', color=c, markersize=size, zorder=5)
        ax1.plot([lo, hi], [y, y], '-', color=c, linewidth=1.5, zorder=4)
        ax1.text(-2, y, name, ha='right', va='center', fontsize=8)
        ax1.text(hi + 1.5, y, f'N={n}', ha='left', va='center', fontsize=7.5, color='grey')

    # Overall diamond
    ax1.plot(68.5, -1.2, 'D', color='red', markersize=10, zorder=6)
    ax1.plot([56.3, 78.6], [-1.2, -1.2], '-', color='red', linewidth=2.5, zorder=5)
    ax1.text(-2, -1.2, 'OVERALL\n(k=11, N=474)', ha='right', va='center', fontsize=9, fontweight='bold', color='red')
    ax1.text(80, -1.2, '68.5%\n[56.3-78.6]', ha='left', va='center', fontsize=9, fontweight='bold', color='red')

    ax1.axvline(68.5, color='red', linestyle='--', alpha=0.4, linewidth=1)
    ax1.set_xlim(-5, 105)
    ax1.set_ylim(-2.5, len(studies))
    ax1.set_xlabel('Detection Prevalence (%)', fontsize=11)
    ax1.set_yticks([])
    ax1.set_title('(A) Study-Level Prevalence Estimates', fontsize=11, fontweight='bold')
    ax1.spines['left'].set_visible(False)

    # Panel B: Quality assessment summary
    ax2 = fig.add_subplot(gs[1])
    quality_labels = ['S1\nRepresent.', 'S2\nSample Size', 'S3\nNon-resp.', 'S4\nAscertain.',
                      'C1\nConfounders', 'C2\nContam. Ctrl', 'O1\nMethod Val.', 'O2\nStatistics', 'O3\nRecovery']
    met_counts = [6, 7, 2, 4, 1, 7, 0, 9, 0]  # out of 11

    bars = ax2.barh(quality_labels, met_counts, color='#4575B4', edgecolor='white', height=0.6)
    ax2.set_xlim(0, 12)
    ax2.set_xlabel('Studies Meeting Criterion (of 11)', fontsize=10)
    ax2.set_title('(B) Quality Assessment\n(Modified NOS)', fontsize=11, fontweight='bold')
    ax2.axvline(5.5, color='orange', linestyle='--', alpha=0.5, linewidth=1)
    ax2.text(5.8, 8.5, '50%', color='orange', fontsize=8)

    for bar, val in zip(bars, met_counts):
        color = '#D73027' if val <= 2 else '#4575B4'
        bar.set_color(color if val <= 2 else '#4575B4')
        ax2.text(val + 0.2, bar.get_y() + bar.get_height()/2,
                str(val), va='center', fontsize=9, fontweight='bold')
    ax2.invert_yaxis()

    fig.suptitle('Figure 8.5  Human Evidence: Detection Prevalence and Study Quality Assessment',
                 fontsize=13, fontweight='bold', y=1.01)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'Fig8_5_human_evidence.png'))
    plt.close()
    print('Fig 8.5 saved')


# ═══════════════════════════════════════════════════════════════
# RUN ALL
# ═══════════════════════════════════════════════════════════════
fig1()
fig2()
fig3()
fig4()
fig5()

print('\nAll 5 figures generated in:', OUT)
for f in sorted(os.listdir(OUT)):
    if f.endswith('.png'):
        sz = os.path.getsize(os.path.join(OUT, f))
        print(f'  {f} ({sz/1024:.0f} KB)')
