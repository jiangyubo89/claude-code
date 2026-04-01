"""
Fig.8-4: iMNPs PBPK模型隔室结构图
PBPK compartment model for iMNPs with key equations
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.patheffects as pe

plt.rcParams['font.family'] = 'Microsoft YaHei'
plt.rcParams['axes.unicode_minus'] = False

fig = plt.figure(figsize=(16, 11))
fig.patch.set_facecolor('#FAFAFA')

# ── 左侧：PBPK模型结构图 ─────────────────────────────────────────────
ax = fig.add_axes([0.02, 0.05, 0.58, 0.90])
ax.set_xlim(0, 12)
ax.set_ylim(0, 14)
ax.axis('off')
ax.set_facecolor('#F0F4FF')

# 隔室定义
compartments = {
    # name: (x_center, y_center, width, height, color, text)
    'oral':      (6.0, 13.0, 3.5, 0.9,  '#FFF8E1', '经口摄入 (Oral Intake)\n↓胃肠道吸收 (GI Absorption)'),
    'inhal':     (1.5, 13.0, 3.0, 0.9,  '#E3F2FD', '经肺吸入 (Inhalation)\n↓肺泡沉积 (Alveolar Deposition)'),
    'lung':      (1.5, 11.0, 2.8, 1.1,  '#BBDEFB', '肺\n(Lung)'),
    'gi':        (6.0, 11.0, 3.0, 1.1,  '#FFF9C4', '胃肠道\n(GI Tract)'),
    'art_blood': (4.0,  8.8, 3.2, 1.0,  '#FFCDD2', '动脉血\n(Arterial Blood)'),
    'ven_blood': (4.0,  6.8, 3.2, 1.0,  '#EF9A9A', '静脉血\n(Venous Blood)'),
    'liver':     (8.5,  8.3, 2.5, 1.2,  '#C8E6C9', '肝脏 (Liver)\n胆汁排泄→粪便'),
    'kidney':    (8.5,  5.8, 2.5, 1.2,  '#DCEDC8', '肾脏 (Kidney)\n尿液排泄→尿'),
    'spleen':    (1.5,  8.3, 2.5, 1.2,  '#D1C4E9', '脾脏\n(Spleen)'),
    'brain':     (1.5,  5.8, 2.5, 1.2,  '#B2EBF2', '脑 (Brain)\nBBB: <50 nm'),
    'placenta':  (5.5,  3.8, 2.8, 1.2,  '#FCE4EC', '胎盘 (Placenta)\n胎儿: <240 nm'),
    'other':     (1.5,  3.8, 2.5, 1.2,  '#F5F5F5', '其他组织\n(Other Tissues)'),
    'feces':     (10.5, 11.5, 1.6, 0.8, '#EFEBE9', '粪便\nFeces'),
    'urine':     (10.5, 5.8, 1.6, 0.8,  '#E8EAF6', '尿液\nUrine'),
}

def draw_box(ax, cx, cy, w, h, color, text, fontsize=8.2):
    box = FancyBboxPatch((cx - w/2, cy - h/2), w, h,
                         boxstyle='round,pad=0.12',
                         facecolor=color, edgecolor='#455A64', linewidth=1.3)
    ax.add_patch(box)
    ax.text(cx, cy, text, ha='center', va='center', fontsize=fontsize,
            fontweight='bold', color='#212121', linespacing=1.35)

for key, val in compartments.items():
    cx, cy, w, h, col, txt = val
    draw_box(ax, cx, cy, w, h, col, txt)

def arrow(ax, x1, y1, x2, y2, label='', color='#455A64', style='->'):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle=style, color=color,
                                lw=1.8, connectionstyle='arc3,rad=0'))
    if label:
        mx, my = (x1+x2)/2, (y1+y2)/2
        ax.text(mx+0.1, my+0.1, label, fontsize=7, color=color,
                ha='center', style='italic')

# 暴露→吸收
arrow(ax, 1.5, 12.55, 1.5, 11.55, 'Q_lung')
arrow(ax, 6.0, 12.55, 6.0, 11.55, 'Q_gi')
# GI → 门静脉 → 肝 (portal circulation)
arrow(ax, 6.0, 10.45, 8.5, 9.0, 'portal vein')
# 肺 → 动脉血
arrow(ax, 1.5, 10.45, 3.2, 9.3, 'Q_pul')
# 动脉血 → 各器官
arrow(ax, 5.2, 8.8, 1.8, 9.0, '')
arrow(ax, 5.2, 8.8, 1.8, 6.4, 'Q_br')
arrow(ax, 5.0, 8.3, 5.0, 7.8, 'Q_circ', color='#B71C1C')
# 各器官 → 静脉血
arrow(ax, 3.2, 8.3, 4.0, 7.3, '')
arrow(ax, 3.2, 6.4, 4.0, 7.2, '')
arrow(ax, 8.5, 7.7, 5.6, 7.3, '')
# 静脉血 → 肺 (循环)
arrow(ax, 4.0, 6.3, 2.0, 10.8, 'venous return', color='#1565C0')
# 肾 → 尿液
arrow(ax, 9.75, 6.4, 9.75, 6.2, '')
arrow(ax, 9.75, 6.1, 10.35, 6.0, '')
# 肝 → 胆汁→粪便
arrow(ax, 9.75, 9.0, 9.75, 11.5, '')
arrow(ax, 9.75, 11.5, 10.4, 11.6, '')
# 动脉血 → 胎盘
arrow(ax, 5.0, 8.3, 5.5, 4.4, 'Q_pl', color='#AD1457')
# 动脉血 → 其他
arrow(ax, 3.5, 8.5, 1.8, 4.6, '')
# 脾
arrow(ax, 4.3, 8.8, 2.8, 8.8, '')
arrow(ax, 2.8, 7.7, 3.5, 7.3, '')

# 肾脏 ← 动脉血
arrow(ax, 5.6, 8.3, 8.5, 6.9, 'Q_kid')

# 图例
ax.text(6.0, 0.8, 'Q = 器官血流量（mL/min）\nC = 隔室颗粒浓度（µg/mL）\nP = 组织-血液分配系数\nk_clear = 清除速率常数（1/min）',
        ha='center', va='bottom', fontsize=8, color='#37474F',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='#FFFDE7',
                  edgecolor='#F9A825', linewidth=1.2))

ax.set_title('(A) iMNPs PBPK模型隔室结构图\n'
             'Compartmental Structure of iMNPs PBPK Model',
             fontsize=11, fontweight='bold', color='#1A237E', pad=8)

# ── 右侧：关键方程 ────────────────────────────────────────────────────
ax_eq = fig.add_axes([0.62, 0.05, 0.36, 0.90])
ax_eq.axis('off')
ax_eq.set_facecolor('#F9FBE7')

equations_text = [
    ('通用质量平衡方程', 'General Mass Balance Equation',
     r'$\frac{dC_i}{dt} = Q_i \cdot C_{blood} - Q_i \cdot \frac{C_i}{P_i} - k_{clear,i} \cdot V_i \cdot C_i$',
     '#1A237E'),
    ('MPS饱和动力学', 'MPS Saturation Kinetics (Michaelis-Menten)',
     r'$V_{uptake} = \frac{V_{max} \cdot C_{blood}}{K_m + C_{blood}}$',
     '#1B5E20'),
    ('经口吸收', 'GI Absorption Rate',
     r'$\dot{m}_{abs} = k_{abs} \cdot f_{abs}(d_p) \cdot m_{GI}$',
     '#4A148C'),
    ('肺泡沉积', 'Pulmonary Deposition',
     r'$\dot{m}_{dep} = V_E \cdot C_{air} \cdot DF(d_p)$',
     '#0D47A1'),
    ('蛋白冠动力学', 'Protein Corona Kinetics',
     r'$\frac{d\Gamma}{dt} = k_{on} \cdot C_{prot} \cdot (1-\theta) - k_{off} \cdot \theta$',
     '#BF360C'),
    ('粒径分布加权', 'PSD-weighted Total Distribution',
     r'$C_{organ} = \int_{d_{min}}^{d_{max}} C_{organ}(d_p) \cdot f(d_p) \, dd_p$',
     '#004D40'),
]

y_pos = 0.95
ax_eq.text(0.5, y_pos, '(B) PBPK关键方程汇总\nKey PBPK Equations', ha='center',
           va='top', transform=ax_eq.transAxes, fontsize=11,
           fontweight='bold', color='#1A237E')

y_pos -= 0.08
for cn, en, eq, col in equations_text:
    ax_eq.add_patch(FancyBboxPatch((0.02, y_pos-0.105), 0.96, 0.12,
                                   boxstyle='round,pad=0.01',
                                   facecolor='white', edgecolor=col,
                                   linewidth=1.5, transform=ax_eq.transAxes, zorder=1))
    ax_eq.text(0.5, y_pos - 0.01, f'{cn}（{en}）', ha='center', va='top',
               transform=ax_eq.transAxes, fontsize=8, fontweight='bold', color=col)
    ax_eq.text(0.5, y_pos - 0.06, eq, ha='center', va='top',
               transform=ax_eq.transAxes, fontsize=10.5, color='#212121')
    y_pos -= 0.145

# 参数表
param_y = y_pos - 0.03
ax_eq.text(0.5, param_y, '主要参数说明 Key Parameters', ha='center',
           va='top', transform=ax_eq.transAxes, fontsize=9, fontweight='bold', color='#37474F')
params = [
    ('$C_i$',     '隔室颗粒浓度 (µg/mL)'),
    ('$Q_i$',     '器官血流量 (mL/min/kg)'),
    ('$P_i$',     '组织-血液分配系数 (dimensionless)'),
    ('$V_{max}$', 'MPS最大摄取速率 (µg/min)'),
    ('$K_m$',     'MPS半饱和浓度 (µg/mL)'),
    ('$k_{abs}$', '胃肠道吸收速率常数 (1/min)'),
    ('$f_{abs}(d_p)$', '粒径依赖性吸收分数'),
    ('$DF(d_p)$', '粒径依赖性肺泡沉积率'),
    ('$\\Gamma$',  '蛋白冠表面覆盖率'),
    ('$f(d_p)$',  '粒径频率分布函数'),
]
for i, (sym, desc) in enumerate(params):
    y_cur = param_y - 0.055 - i * 0.042
    ax_eq.text(0.05, y_cur, sym, va='center', transform=ax_eq.transAxes,
               fontsize=8.5, color='#1565C0')
    ax_eq.text(0.25, y_cur, f': {desc}', va='center', transform=ax_eq.transAxes,
               fontsize=7.8, color='#424242')

fig.suptitle('Fig.8-4  内源性微/纳米塑料（iMNPs）生理药代动力学（PBPK）模型\n'
             'Physiologically Based Pharmacokinetic (PBPK) Model for iMNPs',
             fontsize=13, fontweight='bold', color='#0D47A1', y=0.99)

plt.savefig('D:/桌面/claude/N-MP/figures/Fig8-4_PBPK.png', dpi=300, bbox_inches='tight')
plt.savefig('D:/桌面/claude/N-MP/figures/Fig8-4_PBPK.pdf', bbox_inches='tight')
print("Fig.8-4 saved.")
plt.show()
