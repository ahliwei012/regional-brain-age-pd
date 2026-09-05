# -*- coding: utf-8 -*-
"""Figure 1: 模型验证图。panel B 同时给出原始与年龄偏倚校正后的 MAE 并标明纵轴含义,
内嵌小图显示未校正的原始预测, 使预测斜率的压缩一目了然。
"""
import os
os.environ["MPLCONFIGDIR"] = r"G:/Temp/mpl"
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from PIL import Image
import warnings
warnings.filterwarnings("ignore")

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False,
                     'axes.spines.right': False, 'axes.linewidth': 0.9})

BASE = r"F:/PD_brainage_data"
FIG = BASE + "/SUBMISSION_BrainComms_R1/Figures"
P = {'grey': '#9A9A9A', 'teal': '#0F7B7B', 'blue': '#0F4D92', 'red': '#B22222',
     'dark': '#1A1A1A', 'blue_soft': '#7BA7D4'}

df = pd.read_csv(BASE + "/analysis_revision_base.csv")
hc = df[df.group == "HC"]
x = hc["age_at_visit"].values.astype(float)
y = hc["global"].values.astype(float)
a, b = np.polyfit(x, y, 1)
y_c = y - (a * x + b) + x
r = np.corrcoef(x, y)[0, 1]
mae_raw = np.abs(y - x).mean()
mae_c = np.abs(y_c - x).mean()
print("global: raw MAE %.2f | bias-corrected MAE %.2f | r %.3f | slope %.3f"
      % (mae_raw, mae_c, r, a))

fig = plt.figure(figsize=(6.6, 5.9))
gs = fig.add_gridspec(2, 1, height_ratios=[0.60, 1.40], hspace=0.20,
                      left=0.12, right=0.97, top=0.96, bottom=0.08)

# ---------- A 流程 ----------
axa = fig.add_subplot(gs[0])
axa.axis('off')
axa.set_xlim(0, 10)
axa.set_ylim(0, 3)
steps = ["T1-weighted\nMRI", "CIVET cortical\nfeatures", "Graph conv.\nnetwork",
         "10 network\nbrain ages", "Regional\nbrain-PAD"]
cols = [P['grey'], P['teal'], P['blue'], P['blue'], P['red']]
xs = np.linspace(1.15, 8.85, 5)
for xi, s, c in zip(xs, steps, cols):
    axa.add_patch(FancyBboxPatch((xi - 0.76, 1.00), 1.52, 1.10,
                                 boxstyle="round,pad=0.06", fc='white', ec=c, lw=1.7))
    axa.text(xi, 1.55, s, ha='center', va='center', fontsize=7.2, color=P['dark'])
for x0, x1 in zip(xs[:-1], xs[1:]):
    axa.annotate('', xy=(x1 - 0.80, 1.55), xytext=(x0 + 0.80, 1.55),
                 arrowprops=dict(arrowstyle='-|>', color=P['dark'], lw=1.3))
axa.text(5, 2.72, 'n = 542 patients + 261 controls   (PPMI, multi-site; 1.5 T and 3 T)',
         ha='center', fontsize=8.2, style='italic', color=P['dark'])
axa.text(-0.02, 1.02, 'A', transform=axa.transAxes, fontsize=15,
         fontweight='bold', va='top', color='black')

# ---------- B 校正后散点 ----------
axb = fig.add_subplot(gs[1])
axb.scatter(x, y_c, s=13, alpha=0.5, color=P['blue_soft'], edgecolor='none')
lim = [x.min() - 3, x.max() + 3]
axb.plot(lim, lim, color=P['dark'], ls='--', lw=1)
axb.set_xlim(lim)
axb.set_ylim(lim)
axb.set_xlabel("Chronological age (years)", fontsize=8.6, color='black')
axb.set_ylabel("Age-bias-corrected predicted\nbrain age (years)", fontsize=8.6, color='black')
axb.text(0.035, 0.965,
         "Healthy controls, n = %d\nr = %.2f\nMean absolute error\n"
         "   raw  %.2f years\n   after bias correction  %.2f years"
         % (len(hc), r, mae_raw, mae_c),
         transform=axb.transAxes, va='top', fontsize=8.2, color='black', linespacing=1.4)
axb.text(-0.13, 1.02, 'B', transform=axb.transAxes, fontsize=15,
         fontweight='bold', va='top', color='black')
axb.set_title("Brain-age model validation in healthy controls", fontsize=9.5,
              loc='left', fontweight='bold', color='black')
axb.tick_params(colors='black', labelsize=8)

# 内嵌小图: 未校正的原始预测, 显示斜率压缩
axi = axb.inset_axes([0.615, 0.115, 0.365, 0.375])
axi.scatter(x, y, s=5, alpha=0.45, color=P['grey'], edgecolor='none')
axi.plot(lim, lim, color=P['dark'], ls='--', lw=0.8)
xx = np.array(lim)
axi.plot(xx, a * xx + b, color=P['red'], lw=1.2)
axi.set_xlim(lim)
axi.set_ylim(lim)
axi.set_title("Before bias correction", fontsize=7, color='black', pad=2)
axi.set_xlabel("Chronological age", fontsize=6.2, color='black', labelpad=1)
axi.set_ylabel("Predicted", fontsize=6.2, color='black', labelpad=1)
axi.tick_params(labelsize=5.8, colors='black', length=2)
axi.text(0.04, 0.96, "slope %.2f" % a, transform=axi.transAxes, va='top',
         fontsize=6.2, color=P['red'])

for f in ['png', 'pdf']:
    fig.savefig(FIG + "/Figure_1." + f, dpi=400 if f == 'png' else None)
plt.close(fig)
im = Image.open(FIG + "/Figure_1.png")
d = im.info.get('dpi', (400, 400))[0]
print("Figure_1: %.0f x %.0f mm" % (im.size[0] / d * 25.4, im.size[1] / d * 25.4))
