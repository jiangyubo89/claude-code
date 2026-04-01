"""
Fig.8-5: 人体样本iMNPs检测里程碑时间线（2016-2024）
Timeline of key milestones in human iMNPs detection studies
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np

plt.rcParams['font.family'] = 'Microsoft YaHei'
plt.rcParams['axes.unicode_minus'] = False

fig, ax = plt.subplots(figsize=(18, 9))
fig.patch.set_facecolor('#FAFAFA')
ax.set_facecolor('#F9FBE7')

# 里程碑事件
milestones = [
    {'year': 2015, 'month': 6,  'tissue': '粪便\n(Feces)',        'study': 'Feldman et al.',
     'finding': '首次人类粪便MP检测\n(First human fecal MP)',       'color': '#8BC34A', 'y': 0.7},
    {'year': 2018, 'month': 10, 'tissue': '粪便\n(Feces)',        'study': 'Schwabl et al.',
     'finding': '8国志愿者粪便全部检出\n(All 10 subjects positive)','color': '#8BC34A', 'y': 0.3},
    {'year': 2019, 'month': 3,  'tissue': '肺\n(Lung)',           'study': 'Amato-Lourenço et al.',
     'finding': '人肺尸检组织MP检出\n(MP in autopsy lung)',          'color': '#03A9F4', 'y': 0.7},
    {'year': 2020, 'month': 5,  'tissue': '胎盘\n(Placenta)',     'study': 'Ragusa et al.',
     'finding': '人胎盘首次MP检出\n(First MP in human placenta)',    'color': '#E91E63', 'y': 0.3},
    {'year': 2021, 'month': 1,  'tissue': '胎盘\n(Placenta)',     'study': 'Ragusa et al. 2021',
     'finding': '正式发表胎盘MP证据\n(Confirmed placental MPs)',     'color': '#E91E63', 'y': 0.8},
    {'year': 2021, 'month': 8,  'tissue': '肺\n(Lung)',           'study': 'Jenner et al.',
     'finding': '手术患者肺组织检出\n(Surgical patients lung MPs)',  'color': '#03A9F4', 'y': 0.4},
    {'year': 2022, 'month': 3,  'tissue': '血液\n(Blood)',        'study': 'Leslie et al.',
     'finding': '健康人血液MP\n(MPs in blood, 77%)',                 'color': '#F44336', 'y': 0.7},
    {'year': 2022, 'month': 6,  'tissue': '肺\n(Lung)',           'study': 'Jenner et al.',
     'finding': 'BAL及肺组织双重确认\n(BAL + tissue confirmed)',     'color': '#03A9F4', 'y': 0.25},
    {'year': 2022, 'month': 9,  'tissue': '母乳\n(Breast Milk)', 'study': 'Ragusa et al.',
     'finding': '人母乳中MP检出\n(MPs in human breast milk)',        'color': '#FF9800', 'y': 0.7},
    {'year': 2022, 'month': 11, 'tissue': '肝脏\n(Liver)',        'study': 'Horvatits et al.',
     'finding': '肝活检MP检出(58%)\n(Liver biopsy, 58%)',           'color': '#9C27B0', 'y': 0.4},
    {'year': 2023, 'month': 2,  'tissue': '胎盘\n(Placenta)',     'study': 'Braun et al.',
     'finding': '62例胎盘全部检出\n(All 62 placentas positive)',     'color': '#E91E63', 'y': 0.15},
    {'year': 2024, 'month': 3,  'tissue': '睾丸\n(Testis)',       'study': 'Garcia et al.',
     'finding': '睾丸MP检出(高浓度)\n(Testis: ~327 µg/g DW)',       'color': '#607D8B', 'y': 0.65},
    {'year': 2024, 'month': 7,  'tissue': '动脉壁\n(Artery)',     'study': 'Marfella et al.',
     'finding': '颈动脉粥样斑块中MP\n(MPs in carotid plaques)',      'color': '#795548', 'y': 0.35},
]

years = [m['year'] + (m['month']-1)/12 for m in milestones]
min_year, max_year = 2014.5, 2025

# 主时间轴
ax.axhline(y=0.5, xmin=0, xmax=1, color='#37474F', linewidth=3.5, zorder=2)
ax.set_xlim(min_year, max_year)
ax.set_ylim(0, 1)

# 年份刻度
for yr in range(2015, 2025):
    xpos = yr
    ax.axvline(x=xpos, ymin=0.47, ymax=0.53, color='#455A64', linewidth=1.5, zorder=3)
    ax.text(xpos, 0.44, str(yr), ha='center', va='top', fontsize=9,
            fontweight='bold', color='#37474F')

# 绘制各里程碑
for i, (m, yr) in enumerate(zip(milestones, years)):
    y = m['y']
    color = m['color']

    # 连接线
    ax.plot([yr, yr], [0.5, y], color=color, linewidth=1.5,
            linestyle='--', alpha=0.7, zorder=3)

    # 圆圈节点
    ax.scatter([yr], [0.5], s=80, color=color, zorder=5, edgecolors='white', linewidths=1.5)

    # 信息框
    box_width = 0.85
    box_height = 0.14
    bx = yr - box_width/2
    by = y - box_height/2

    box = FancyBboxPatch((bx, by), box_width, box_height,
                         boxstyle='round,pad=0.015',
                         facecolor=color, edgecolor='white',
                         alpha=0.85, linewidth=1.2, zorder=4)
    ax.add_patch(box)

    # 文字
    ax.text(yr, y + 0.03, m['tissue'], ha='center', va='center',
            fontsize=8.5, fontweight='bold', color='white', zorder=6)
    ax.text(yr, y - 0.015, m['finding'], ha='center', va='center',
            fontsize=7, color='white', zorder=6, linespacing=1.3)
    ax.text(yr, y - 0.055, m['study'], ha='center', va='center',
            fontsize=6.5, color='#F5F5F5', style='italic', zorder=6)

# 图例
tissue_colors = {
    '粪便 Feces': '#8BC34A', '肺 Lung': '#03A9F4',
    '血液 Blood': '#F44336', '胎盘 Placenta': '#E91E63',
    '母乳 Breast Milk': '#FF9800', '肝脏 Liver': '#9C27B0',
    '睾丸 Testis': '#607D8B', '动脉 Artery': '#795548'
}
legend_handles = [mpatches.Patch(color=c, label=l) for l, c in tissue_colors.items()]
ax.legend(handles=legend_handles, loc='upper left', bbox_to_anchor=(0.01, 0.99),
          fontsize=8.5, ncol=2, framealpha=0.9, edgecolor='#BDBDBD', title='组织类型 Tissue Type')

ax.set_xlim(min_year, max_year)
ax.set_ylim(0, 1.0)
ax.axis('off')

ax.set_title('Fig.8-5  人体组织iMNPs检测里程碑时间线（2015–2024）\n'
             'Timeline of Key Milestones in Human Tissue iMNPs Detection Studies (2015–2024)',
             fontsize=13, fontweight='bold', color='#0D47A1', pad=15)

fig.text(0.5, 0.01,
         '注：时间轴位置代表研究发表时间；文献信息经文献核实整理；BAL=支气管肺泡灌洗液；DW=干重\n'
         'Note: Positions represent publication dates; BAL=bronchoalveolar lavage; DW=dry weight.',
         ha='center', fontsize=8, color='#546E7A', style='italic')

plt.tight_layout(rect=[0, 0.04, 1, 1])
plt.savefig('D:/桌面/claude/N-MP/figures/Fig8-5_timeline.png', dpi=300, bbox_inches='tight')
plt.savefig('D:/桌面/claude/N-MP/figures/Fig8-5_timeline.pdf', bbox_inches='tight')
print("Fig.8-5 saved.")
plt.show()
