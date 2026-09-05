# -*- coding: utf-8 -*-
"""Brain Communications 图形摘要
按 120 x 120 mm 实际展示尺寸设计(而非缩放大画布), 字号 Arial 12-16 pt,
渲染为 3000 x 3000 px (1:1)。按期刊要求不使用缩写。
"""
import os
os.environ["MPLCONFIGDIR"] = r"G:/Temp/mpl"
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrow
from PIL import Image

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['pdf.fonttype'] = 42

OUT = r"F:/PD_brainage_data/SUBMISSION_BrainComms_R1/Figures"
MM = 120.0                      # 展示尺寸(mm), 期刊规定
IN = MM / 25.4                  # 4.724 in
DPI = int(round(3000 / IN))     # -> 3000 px
INK, GREY = '#1A1A1A', '#5A5A5A'
RED, BLUE, GREEN, TEAL = '#B22222', '#0F4D92', '#1E6B34', '#0F7B7B'

fig = plt.figure(figsize=(IN, IN), facecolor='white')
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis('off')


def box(x, y, w, h, fc, ec, r=1.6, lw=1.4, alpha=1.0):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=%.1f" % r,
                                facecolor=fc, edgecolor=ec, linewidth=lw, alpha=alpha))


# ---------------- 标题 ----------------
ax.text(50, 97.5, "Regional cortical brain age in", ha='center', va='top',
        fontsize=16, fontweight='bold', color=INK)
ax.text(50, 92.0, "Parkinson's disease", ha='center', va='top',
        fontsize=16, fontweight='bold', color=INK)
ax.text(50, 86.4, "803 participants  ·  542 patients, 261 controls",
        ha='center', va='top', fontsize=12, color=GREY)

# ---------------- 流程 ----------------
yb, hb = 72.8, 8.0
for x, w, t, c in [(4, 26, "Structural\nbrain scan", GREY),
                   (37, 26, "Graph neural\nnetwork", TEAL),
                   (70, 26, "Brain age for\n10 networks", BLUE)]:
    box(x, yb, w, hb, 'white', c)
    ax.text(x + w / 2, yb + hb / 2, t, ha='center', va='center',
            fontsize=12, color=INK, linespacing=1.15)
for x in (31.2, 64.2):
    ax.add_patch(FancyArrow(x, yb + hb / 2, 4.2, 0, width=0.5, head_width=2.0,
                            head_length=1.6, length_includes_head=True,
                            color=GREY, linewidth=0))

# ---------------- 发现 ----------------
ax.text(50, 68.6, "Older brain age, worse clinical status",
        ha='center', va='top', fontsize=13, fontstyle='italic', color=INK)

box(3.5, 41.5, 45, 21.5, '#FBF0F0', RED, lw=1.6)
ax.text(26, 61.2, "Motor severity", ha='center', va='top',
        fontsize=13.5, fontweight='bold', color=RED)
ax.text(26, 49.5, "visual, ventral attention,\nsalience and language\nnetworks",
        ha='center', va='center', fontsize=12, color=INK, linespacing=1.3)

box(51.5, 41.5, 45, 21.5, '#EEF3FA', BLUE, lw=1.6)
ax.text(74, 61.2, "Cognition", ha='center', va='top',
        fontsize=13.5, fontweight='bold', color=BLUE)
ax.text(74, 49.5, "frontoparietal\nnetwork", ha='center', va='center',
        fontsize=12, color=INK, linespacing=1.3)

ax.text(50, 38.6, "Whole-brain brain age was associated with neither",
        ha='center', va='top', fontsize=12, fontstyle='italic', color=GREY)

# ---------------- 独立性 ----------------
box(3.5, 16.0, 93, 17.8, '#F2F8F3', GREEN, lw=1.6)
ax.text(50, 32.4, "Independent of", ha='center', va='top',
        fontsize=13.5, fontweight='bold', color=GREEN)
ax.text(50, 24.3, "cortical thickness  ·  dopamine transporter binding",
        ha='center', va='center', fontsize=12, color=INK)
ax.text(50, 19.6, "alpha-synuclein seeding  ·  amyloid and tau",
        ha='center', va='center', fontsize=12, color=INK)

# ---------------- 结论 ----------------
ax.text(50, 11.6, "A marker of current cortical status,", ha='center', va='top',
        fontsize=13, fontweight='bold', color=INK)
ax.text(50, 6.4, "not of subsequent progression", ha='center', va='top',
        fontsize=13, fontweight='bold', color=INK)

fig.savefig(OUT + "/Graphical_Abstract.png", dpi=DPI, facecolor='white')
fig.savefig(OUT + "/Graphical_Abstract.pdf", facecolor='white')
plt.close(fig)

im = Image.open(OUT + "/Graphical_Abstract.png").convert("RGB")
im.save(OUT + "/Graphical_Abstract.tif", format="TIFF", compression="tiff_lzw", dpi=(DPI, DPI))
print("Graphical abstract:", im.size, "aspect", im.size[0] / im.size[1],
      "| 展示尺寸 %.0f x %.0f mm | 有效字号 12-16 pt" % (MM, MM))
