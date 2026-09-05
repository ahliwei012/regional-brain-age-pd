# -*- coding: utf-8 -*-
"""修正补充图的排版与数据问题:
S1  标题与首个方框重叠 -> 扩大上方留白
S7  (a) 图例压住 Sensorimotor 数据; (b) 仍是重建前的旧纵向结果(n≈195)
    -> 用 revision_longitudinal.py 的新结果重绘(n=349 / n=440)
S5  图例标 ComBat 为蓝色菱形, 实际绘制为灰色 -> 统一配色
"""
import os
os.environ["MPLCONFIGDIR"] = r"G:/Temp/mpl"
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrow
from matplotlib.lines import Line2D
from PIL import Image
import warnings
warnings.filterwarnings("ignore")
Image.MAX_IMAGE_PIXELS = None

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False,
                     'axes.spines.right': False, 'axes.linewidth': 0.9})

BASE = r"F:/PD_brainage_data"
FIG = BASE + "/SUBMISSION_BrainComms_R1/Figures"
INK, GREY, ORANGE, BLUE, RED = '#1A1A1A', '#9A9A9A', '#D55E00', '#0F4D92', '#B22222'
MAXW, MAXH = 170.0, 210.0


def report(name):
    im = Image.open(FIG + "/%s.png" % name)
    dpi = im.info.get('dpi', (400, 400))[0]
    w, h = im.size[0] / dpi * 25.4, im.size[1] / dpi * 25.4
    k = max(w / MAXW, h / MAXH, 1.0)
    if k > 1.0:
        im.save(FIG + "/%s.png" % name, dpi=(dpi * k, dpi * k))
    print("  %-12s %5.0f x %5.0f mm" % (name, w / k, h / k))


# ==================== S1 参与者流程图 ====================
fig = plt.figure(figsize=(6.6, 3.6))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis('off')


def box(x, y, w, h, text, ec, fc='white', fs=8.4, tc=None):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle="round,pad=0,rounding_size=1.5",
                                facecolor=fc, edgecolor=ec, linewidth=1.4))
    ax.text(x + w / 2, y + h / 2, text, ha='center', va='center',
            fontsize=fs, color=tc or INK, linespacing=1.35)


ax.text(50, 96.5, "Participant flow", fontsize=11.5, fontweight='bold',
        ha='center', va='top', color=INK)

box(14, 71.0, 52, 15.0,
    "PPMI participants with a baseline T1-weighted scan\n"
    "submitted for CIVET cortical processing\n(560 patients; 282 controls; n = 842)", INK)
box(14, 42.0, 52, 15.0,
    "Successfully processed with valid cortical thickness\n"
    "and grey/white matter features\n(550 patients; 265 controls; n = 815)", INK)
box(14, 13.0, 52, 15.0,
    "Final analytic sample\n542 patients with Parkinson's disease\n"
    "261 healthy controls (n = 803)", BLUE)

box(69, 60.0, 30, 13.0,
    "Excluded: CIVET processing\nfailed to yield cortical thickness\n"
    "(10 patients; 17 controls; n = 27)", RED, '#FCF2F2', 7.6, RED)
box(69, 31.0, 30, 13.0,
    "Excluded: missing MDS-UPDRS III,\nMoCA or other required clinical\n"
    "variables (8 patients; 4 controls; n = 12)", RED, '#FCF2F2', 7.6, RED)

for y0, y1 in [(71.0, 57.0), (42.0, 28.0)]:
    ax.add_patch(FancyArrow(40, y0, 0, y1 - y0 + 0.6, width=0.35, head_width=1.6,
                            head_length=1.8, length_includes_head=True,
                            color=INK, linewidth=0))
for ys, ye in [(66.5, 66.5), (37.5, 37.5)]:
    ax.add_patch(FancyArrow(40, ys, 27.5, 0, width=0.3, head_width=1.4,
                            head_length=1.5, length_includes_head=True,
                            color=RED, linewidth=0))

ax.text(50, 6.5, "Numbers reconciled across CIVET processing logs, the feature-extraction output\n"
                 "directory (815 files) and the final analytic dataset (803 rows).",
        ha='center', va='center', fontsize=7.2, style='italic', color='#555')
for f in ['png', 'pdf']:
    fig.savefig(FIG + "/Figure_S1." + f, dpi=400 if f == 'png' else None)
plt.close(fig)
print("S1 重绘 (标题与方框分离):")
report("Figure_S1")

# ==================== S7 纵向(用重建后的结果) ====================
U = pd.read_csv(BASE + "/rev_long_updrs3_score.csv")
M = pd.read_csv(BASE + "/rev_long_moca.csv")
nU, nM = int(U["n_subj"].max()), int(M["n_subj"].max())

fig, axes = plt.subplots(1, 2, figsize=(6.6, 3.5),
                         gridspec_kw={'wspace': 0.10, 'left': 0.20, 'right': 0.98,
                                      'top': 0.74, 'bottom': 0.20})
for ax, (D, ttl, xl) in zip(axes, [
        (U, "A  Motor progression (MDS-UPDRS III, n = %d)" % nU,
         "\u03b2  time \u00d7 brain-PAD  (points per SD per year)"),
        (M, "B  Cognitive progression (MoCA, n = %d)" % nM,
         "\u03b2  time \u00d7 brain-PAD  (points per SD per year)")]):
    D = D.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(D))
    for i in range(len(D)):
        nom = D['p'][i] < 0.05
        col = ORANGE if nom else GREY
        ax.plot([D['lo'][i], D['hi'][i]], [i, i], color=col, lw=1.6, solid_capstyle='round')
        ax.scatter(D['beta_time_x_PAD'][i], i, s=30 if nom else 22, color=col,
                   zorder=3, edgecolor='white', lw=0.5)
    ax.axvline(0, color='#444', ls=(0, (4, 3)), lw=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels(D['Network'], fontsize=7.6)
    ax.set_xlabel(xl, fontsize=7.4, color='black')
    ax.set_title(ttl, fontsize=8.4, loc='left', fontweight='bold', pad=6, color='black')
    ax.tick_params(colors='black', labelsize=7.4)
axes[1].set_yticklabels([])
axes[0].tick_params(axis='y', length=0)
axes[1].tick_params(axis='y', length=0)

fig.suptitle("Baseline regional brain age does not predict the rate of clinical change",
             fontsize=9.6, fontweight='bold', x=0.02, ha='left', y=0.975, color='black')
# 图例移到图外, 不再压住数据
fig.legend([Line2D([0], [0], color=ORANGE, marker='o', lw=1.6),
            Line2D([0], [0], color=GREY, marker='o', lw=1.6)],
           ['nominal P < 0.05 (not significant after correction)',
            'P \u2265 0.05'],
           loc='upper center', bbox_to_anchor=(0.55, 0.925), ncol=2, fontsize=6.8,
           frameon=False, handletextpad=0.4, columnspacing=1.4, labelcolor='black')
for f in ['png', 'pdf']:
    fig.savefig(FIG + "/Figure_S7." + f, dpi=400 if f == 'png' else None)
plt.close(fig)
print("\nS7 重绘 (新纵向结果 n=%d / n=%d, 图例移出绘图区):" % (nU, nM))
print("   所有网络 FDR: UPDRS-III min=%.3f, MoCA min=%.3f"
      % (U['fdr'].min(), M['fdr'].min()))
report("Figure_S7")
