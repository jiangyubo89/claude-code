"""
Fig.8-1: iMNPs检测技术比较雷达图
Radar chart comparing µ-FTIR, µ-Raman, Py-GC/MS, ATR-FTIR, TEM-EDX
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

plt.rcParams['font.family'] = 'Microsoft YaHei'
plt.rcParams['axes.unicode_minus'] = False

# 评估维度（中英双语）
categories = [
    '灵敏度\n(Sensitivity)',
    '空间分辨率\n(Spatial Resolution)',
    '检测通量\n(Throughput)',
    '纳米级适用性\n(Nano Applicability)',
    '形态信息\n(Morphology Info)',
    '定量能力\n(Quantification)',
    '生物基质耐受\n(Matrix Tolerance)'
]

N = len(categories)
angles = np.linspace(0, 2 * np.pi, N, endpoint=False).tolist()
angles += angles[:1]

# 各技术评分（0–10）
tech_data = {
    'µ-FTIR':       [7, 5, 8, 2, 9, 7, 4],
    'µ-Raman':      [8, 9, 5, 7, 9, 6, 3],
    'Py-GC/MS':     [9, 1, 3, 10, 1, 10, 9],
    'ATR-FTIR':     [6, 3, 9, 1, 5, 5, 6],
    'TEM-EDX':      [10, 10, 1, 10, 10, 3, 2],
}

colors = ['#2196F3', '#4CAF50', '#FF5722', '#9C27B0', '#FF9800']
linestyles = ['-', '-', '-', '--', '--']

fig, ax = plt.subplots(figsize=(10, 9), subplot_kw=dict(polar=True))
fig.patch.set_facecolor('#FAFAFA')
ax.set_facecolor('#F5F5F5')

for (label, values), color, ls in zip(tech_data.items(), colors, linestyles):
    vals = values + values[:1]
    ax.plot(angles, vals, color=color, linewidth=2.2, linestyle=ls, label=label)
    ax.fill(angles, vals, color=color, alpha=0.08)

ax.set_xticks(angles[:-1])
ax.set_xticklabels(categories, size=9.5, fontweight='bold')
ax.set_yticks([2, 4, 6, 8, 10])
ax.set_yticklabels(['2', '4', '6', '8', '10'], size=8, color='gray')
ax.set_ylim(0, 10)
ax.set_title('Fig.8-1  iMNPs检测技术综合性能比较雷达图\n'
             'Comprehensive Performance Comparison of iMNPs Detection Techniques',
             fontsize=12, fontweight='bold', pad=20, color='#1A237E')

ax.legend(loc='upper right', bbox_to_anchor=(1.38, 1.18),
          fontsize=10, framealpha=0.9, edgecolor='#BDBDBD')

ax.grid(color='gray', linestyle='--', linewidth=0.5, alpha=0.5)

fig.text(0.5, 0.01,
         '评分说明：各维度满分10分；灵敏度指检测下限表现，空间分辨率指最小可分辨粒径，\n'
         '纳米级适用性指<100 nm颗粒的检测能力，矩阵耐受指对复杂生物基质的适应性。\n'
         'Score description: 10-point scale per dimension; Nano applicability refers to <100 nm detection; '
         'Matrix tolerance refers to performance in complex biological matrices.',
         ha='center', fontsize=7.5, color='#546E7A', style='italic',
         wrap=True)

plt.tight_layout(rect=[0, 0.06, 1, 1])
plt.savefig('D:/桌面/claude/N-MP/figures/Fig8-1_radar.png', dpi=300, bbox_inches='tight')
plt.savefig('D:/桌面/claude/N-MP/figures/Fig8-1_radar.pdf', bbox_inches='tight')
print("Fig.8-1 saved.")
plt.show()
