"""
Fig.8-2: iMNPs体内分布解剖示意图（器官积累热力图）
Anatomical biodistribution of iMNPs with organ accumulation heatmap
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch, Circle, Ellipse
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.cm as cm

plt.rcParams['font.family'] = 'Microsoft YaHei'
plt.rcParams['axes.unicode_minus'] = False

fig = plt.figure(figsize=(16, 11))
fig.patch.set_facecolor('#FAFAFA')

# ── 左侧：解剖示意图 ────────────────────────────────────────────────
ax_body = fig.add_axes([0.03, 0.05, 0.48, 0.88])
ax_body.set_xlim(0, 10)
ax_body.set_ylim(0, 14)
ax_body.set_aspect('equal')
ax_body.axis('off')
ax_body.set_facecolor('#F0F4FF')

# 积累量数据（单位：% of given dose，综合文献数据）
organ_data = {
    'GI道\n(GI Tract)':    {'center': (5.0, 6.5), 'accumulation': 45.0, 'size': (2.4, 2.8)},
    '肝脏\n(Liver)':       {'center': (4.0, 8.2), 'accumulation': 12.0, 'size': (1.8, 1.1)},
    '肺\n(Lung)':          {'center': (4.2, 10.5),'accumulation': 5.0,  'size': (1.3, 1.0)},
    '脾脏\n(Spleen)':      {'center': (6.2, 8.2), 'accumulation': 3.5,  'size': (0.7, 0.9)},
    '肾脏\n(Kidney)':      {'center': (3.5, 7.2), 'accumulation': 2.8,  'size': (0.65, 0.9)},
    '脑\n(Brain)':         {'center': (5.0, 12.8),'accumulation': 0.3,  'size': (1.2, 0.9)},
    '胎盘\n(Placenta)':    {'center': (5.0, 4.5), 'accumulation': 0.8,  'size': (1.0, 0.7)},
    '骨髓\n(Bone Marrow)':{  'center': (7.5, 7.0), 'accumulation': 1.2,  'size': (0.8, 0.8)},
}

# 颜色映射（积累量 → 颜色强度）
cmap_heat = LinearSegmentedColormap.from_list(
    'heat', ['#E3F2FD', '#1565C0', '#B71C1C'], N=256)
norm = plt.Normalize(vmin=0, vmax=50)

# 绘制人体轮廓（简化线条）
body_outline_x = [4.5, 4.0, 3.5, 3.2, 3.0, 3.0, 3.5, 4.0, 4.2, 4.5, 5.5, 5.8, 6.0, 7.0, 7.0, 6.8, 6.5, 6.0, 5.5, 5.5, 4.5]
body_outline_y = [13.5,12.5,11.5,10.0, 8.5, 5.5, 3.5, 2.0, 1.5, 1.0, 1.0, 1.5, 2.0, 3.5, 5.5, 8.5,10.0,11.5,12.5,13.5,13.5]

ax_body.fill(body_outline_x, body_outline_y, color='#ECEFF1', alpha=0.7, zorder=1)
ax_body.plot(body_outline_x, body_outline_y, color='#78909C', linewidth=1.5, zorder=2)

# 绘制各器官（椭圆热力图）
for organ_name, data in organ_data.items():
    cx, cy = data['center']
    rx, ry = data['size']
    acc = data['accumulation']
    color = cmap_heat(norm(acc))

    ellipse = Ellipse((cx, cy), rx, ry, angle=0,
                      facecolor=color, edgecolor='white',
                      linewidth=1.5, alpha=0.9, zorder=3)
    ax_body.add_patch(ellipse)

    # 标注器官名称
    fontsize = 7.5 if len(organ_name) < 8 else 6.5
    ax_body.text(cx, cy, f'{organ_name}\n{acc}%',
                 ha='center', va='center', fontsize=fontsize,
                 fontweight='bold', color='white', zorder=4,
                 bbox=dict(boxstyle='round,pad=0.1', facecolor='none', edgecolor='none'))

# 颈部/头部连接
ax_body.plot([4.8, 5.2], [13.5, 12.8], color='#78909C', linewidth=1.5, zorder=2)

ax_body.set_title('(A) 器官积累分布示意图\nOrgan Accumulation Diagram',
                  fontsize=11, fontweight='bold', color='#1A237E', pad=8)

# 颜色标尺
ax_cb = fig.add_axes([0.03, 0.03, 0.44, 0.025])
sm = cm.ScalarMappable(cmap=cmap_heat, norm=norm)
sm.set_array([])
cb = plt.colorbar(sm, cax=ax_cb, orientation='horizontal')
cb.set_label('器官积累量（% of administered dose）\nOrgan Accumulation (% of Administered Dose)',
             fontsize=9)
cb.ax.tick_params(labelsize=8)

# ── 右侧：条形图 ───────────────────────────────────────────────────
ax_bar = fig.add_axes([0.56, 0.12, 0.42, 0.75])
ax_bar.set_facecolor('#F9FBE7')

organs = ['GI道\n(GI Tract)', '肝脏\n(Liver)', '肺\n(Lung)', '脾脏\n(Spleen)',
          '肾脏\n(Kidney)', '骨髓\n(Bone Marrow)', '胎盘\n(Placenta)', '脑\n(Brain)']
values = [45.0, 12.0, 5.0, 3.5, 2.8, 1.2, 0.8, 0.3]
errors = [8.0, 3.5, 1.8, 1.2, 1.0, 0.5, 0.4, 0.15]
colors_bar = [cmap_heat(norm(v)) for v in values]

bars = ax_bar.barh(range(len(organs)), values, xerr=errors, color=colors_bar,
                   capsize=4, ecolor='#546E7A', edgecolor='white', linewidth=1.2,
                   height=0.65, alpha=0.9)

ax_bar.set_yticks(range(len(organs)))
ax_bar.set_yticklabels(organs, fontsize=9.5)
ax_bar.set_xlabel('积累量（% of administered dose）\nOrgan Accumulation (%)', fontsize=10)
ax_bar.set_xlim(0, 60)
ax_bar.set_title('(B) 各器官iMNPs积累量定量比较\n'
                 'Quantitative Comparison of iMNPs Accumulation Across Organs',
                 fontsize=11, fontweight='bold', color='#1A237E', pad=10)
ax_bar.axvline(x=0, color='gray', linewidth=0.8)
ax_bar.grid(axis='x', alpha=0.4, linestyle='--')

# 数值标注
for i, (v, e) in enumerate(zip(values, errors)):
    ax_bar.text(v + e + 0.5, i, f'{v}±{e}%', va='center', fontsize=8.5, color='#37474F')

# 注释（粒径说明）
note_text = ('数据代表<200 nm PS颗粒经口暴露（100 mg/kg）7天后啮齿类器官分布（综合多项研究中位估计值）\n'
             'Data represent organ distribution of <200 nm PS particles after oral exposure (100 mg/kg) '
             'for 7 days in rodents\n(median estimates from multiple studies)')
fig.text(0.56, 0.03, note_text, fontsize=7.5, color='#546E7A', style='italic',
         wrap=True, va='bottom')

fig.suptitle('Fig.8-2  内源性纳米/微塑料（iMNPs）体内分布解剖示意图\n'
             'Anatomical Biodistribution of Endogenous Micro/Nanoplastics (iMNPs)',
             fontsize=13, fontweight='bold', color='#0D47A1', y=0.98)

plt.savefig('D:/桌面/claude/N-MP/figures/Fig8-2_biodistribution.png', dpi=300, bbox_inches='tight')
plt.savefig('D:/桌面/claude/N-MP/figures/Fig8-2_biodistribution.pdf', bbox_inches='tight')
print("Fig.8-2 saved.")
plt.show()
