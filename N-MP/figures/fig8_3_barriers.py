"""
Fig.8-3: 四大生物屏障iMNPs跨越机制图
Four biological barriers crossed by iMNPs
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import matplotlib.patheffects as pe

plt.rcParams['font.family'] = 'Microsoft YaHei'
plt.rcParams['axes.unicode_minus'] = False

fig, axes = plt.subplots(2, 2, figsize=(16, 12))
fig.patch.set_facecolor('#FAFAFA')

barrier_configs = [
    {
        'title': '(A) 肠上皮屏障\nIntestinal Epithelial Barrier',
        'threshold': '粒径阈值 <500 nm',
        'bg_top': '#FFF8E1', 'bg_bot': '#E8F5E9',
        'top_label': '肠腔 (Intestinal Lumen)',
        'bot_label': '固有层 → 血/淋巴 (Lamina Propria → Blood/Lymph)',
        'cell_color': '#FFB74D',
        'mechanisms': ['M细胞吞噬\n(M cell phagocytosis)', 'Clathrin介导内吞\n(CME, <200 nm)',
                       'Caveolae介导内吞\n(CavME, ~50 nm)', '炎症旁细胞\n(Paracellular)'],
        'size_note': '< 500 nm：淋巴吸收\n< 200 nm：血液循环\n< 50 nm：高效跨越',
        'ax_idx': (0, 0)
    },
    {
        'title': '(B) 肺泡-毛细血管屏障\nAlveolar-Capillary Barrier',
        'threshold': '粒径阈值 <100 nm',
        'bg_top': '#E3F2FD', 'bg_bot': '#FCE4EC',
        'top_label': '肺泡腔 (Alveolar Space)',
        'bot_label': '毛细血管 (Pulmonary Capillary)',
        'cell_color': '#64B5F6',
        'mechanisms': ['AT1细胞跨细胞转运\n(Transcellular, AT1 cell)', '表面活性物质层吸附\n(Surfactant adsorption)',
                       '巨噬细胞清除\n(AM clearance, >1 µm)', '损伤后旁细胞\n(Post-injury paracellular)'],
        'size_note': '< 100 nm：系统性分布\n< 50 nm：高效跨越\n> 1 µm：主要被AM清除',
        'ax_idx': (0, 1)
    },
    {
        'title': '(C) 血脑屏障\nBlood-Brain Barrier (BBB)',
        'threshold': '粒径阈值 <50 nm',
        'bg_top': '#FCE4EC', 'bg_bot': '#F3E5F5',
        'top_label': '脑毛细血管腔 (Brain Capillary Lumen)',
        'bot_label': '脑实质 (Brain Parenchyma)',
        'cell_color': '#EF9A9A',
        'mechanisms': ['受体介导内吞（RMT）\n(Receptor-mediated transcytosis)', 'ApoE/ApoJ蛋白冠\n(Lipoprotein corona)',
                       '嗅觉神经轴浆转运\n(Olfactory nerve route)', '神经炎症旁细胞\n(Neuroinflammation)'],
        'size_note': '< 50 nm：跨越完整BBB\n< 100 nm：炎症状态下跨越\n关键：BBB TEER最高',
        'ax_idx': (1, 0)
    },
    {
        'title': '(D) 胎盘屏障\nPlacental Barrier',
        'threshold': '粒径阈值 <240 nm',
        'bg_top': '#E8EAF6', 'bg_bot': '#E0F2F1',
        'top_label': '母体血液 (Maternal Blood)',
        'bot_label': '胎儿循环 (Fetal Circulation)',
        'cell_color': '#9FA8DA',
        'mechanisms': ['合体滋养层胞饮\n(Syncytiotrophoblast macropinocytosis)', 'FcRn受体途径\n(FcRn-mediated)',
                       'Clathrin/Caveolae内吞\n(CME/CavME)', '病理妊娠屏障损伤\n(Pathological disruption)'],
        'size_note': '< 240 nm：跨越证据\n50 nm：最高效转运\n妊娠晚期：屏障变薄(<4 µm)',
        'ax_idx': (1, 1)
    }
]

for cfg in barrier_configs:
    ax = axes[cfg['ax_idx']]
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis('off')

    # 背景分区
    ax.add_patch(FancyBboxPatch((0, 5.2), 10, 4.8, boxstyle='round,pad=0',
                                facecolor=cfg['bg_top'], edgecolor='none', zorder=0))
    ax.add_patch(FancyBboxPatch((0, 0), 10, 5.0, boxstyle='round,pad=0',
                                facecolor=cfg['bg_bot'], edgecolor='none', zorder=0))

    # 屏障细胞层（梯形代表细胞单层）
    cell_y_center = 5.1
    for xi in [1.0, 2.8, 4.6, 6.4, 8.2]:
        cell = plt.Polygon([[xi, cell_y_center+0.65], [xi+1.5, cell_y_center+0.65],
                             [xi+1.7, cell_y_center-0.55], [xi-0.2, cell_y_center-0.55]],
                           facecolor=cfg['cell_color'], edgecolor='white', linewidth=1.5,
                           alpha=0.85, zorder=2)
        ax.add_patch(cell)

    # 紧密连接指示线
    ax.plot([0, 10], [cell_y_center-0.6, cell_y_center-0.6], '--',
            color='#757575', linewidth=0.8, alpha=0.5, zorder=3)

    # 区域标签
    ax.text(5, 9.4, cfg['top_label'], ha='center', va='center',
            fontsize=8.5, color='#37474F', fontweight='bold')
    ax.text(5, 0.6, cfg['bot_label'], ha='center', va='center',
            fontsize=8.5, color='#37474F', fontweight='bold')

    # 转运箭头和颗粒
    particle_colors = ['#E53935', '#7B1FA2', '#1565C0', '#2E7D32']
    particle_sizes_nm = ['20 nm', '50 nm', '100 nm', '500 nm']
    particle_radii = [0.12, 0.16, 0.22, 0.32]

    arrow_x_positions = [1.8, 3.6, 5.4, 7.2]
    for i, (x, pc, ps, pr) in enumerate(zip(arrow_x_positions, particle_colors,
                                            particle_sizes_nm, particle_radii)):
        if i < 3:  # 小颗粒：画向下箭头（跨越）
            ax.annotate('', xy=(x, 4.0), xytext=(x, 6.2),
                        arrowprops=dict(arrowstyle='->', color=pc, lw=2.0))
            circle = plt.Circle((x, 6.4), pr, color=pc, alpha=0.8, zorder=5)
            ax.add_patch(circle)
        else:  # 大颗粒：画阻挡符号（被截留）
            ax.plot(x, 5.1, 'x', color=pc, markersize=14, markeredgewidth=2.5, zorder=5)
            circle = plt.Circle((x, 6.4), pr, color=pc, alpha=0.8, zorder=5)
            ax.add_patch(circle)

        ax.text(x, 7.0, ps, ha='center', fontsize=7.5, color=pc, fontweight='bold')

    # 机制列表
    for j, mech in enumerate(cfg['mechanisms']):
        ax.text(0.2, 4.5 - j*0.95, f'• {mech}', va='center', fontsize=7.2,
                color='#424242', wrap=True)

    # 粒径阈值框
    ax.add_patch(FancyBboxPatch((6.3, 0.1), 3.55, 2.2,
                                boxstyle='round,pad=0.15',
                                facecolor='#FFFDE7', edgecolor='#FFA000',
                                linewidth=1.5, zorder=4))
    ax.text(8.08, 1.6, '粒径临界值', ha='center', fontsize=7.5,
            color='#E65100', fontweight='bold', zorder=5)
    ax.text(8.08, 1.1, cfg['size_note'], ha='center', fontsize=6.8,
            color='#424242', zorder=5, linespacing=1.4)

    # 标题和阈值标注
    ax.set_title(cfg['title'], fontsize=10.5, fontweight='bold',
                 color='#1A237E', pad=6)
    ax.text(5, 9.0, cfg['threshold'], ha='center', fontsize=9,
            color='#B71C1C', fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#FFEBEE',
                      edgecolor='#C62828', linewidth=1.2))

fig.suptitle('Fig.8-3  iMNPs跨越四大生物屏障的机制与粒径阈值\n'
             'Mechanisms and Size Thresholds of iMNPs Crossing Four Biological Barriers',
             fontsize=13, fontweight='bold', color='#0D47A1', y=0.99)

plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.savefig('D:/桌面/claude/N-MP/figures/Fig8-3_barriers.png', dpi=300, bbox_inches='tight')
plt.savefig('D:/桌面/claude/N-MP/figures/Fig8-3_barriers.pdf', bbox_inches='tight')
print("Fig.8-3 saved.")
plt.show()
