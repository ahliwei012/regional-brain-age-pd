# -*- coding: utf-8 -*-
"""补充图图例修正:
S5  图例把 ComBat 标为蓝色菱形, 但因无一结果达 FDR<0.05, 实际绘制全为灰色 -> 图例与实绘一致
S7  重绘后无网络达名义显著, 橙色图例项已无对应数据 -> 去掉, 改为文字说明
"""
import os
os.environ["MPLCONFIGDIR"] = r"G:/Temp/mpl"
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
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
RED, BLUE, GREY, LGREY, ORANGE = '#B22222', '#0F4D92', '#9A9A9A', '#C8C8C8', '#D55E00'
REG = ["sns", "fpn", "drs", "vnt", "dft", "slt", "lng", "adt", "vsl", "lmb"]
EN = {"sns": "Sensorimotor", "fpn": "Frontoparietal", "drs": "Dorsal attention",
      "vnt": "Ventral attention", "dft": "Default mode", "slt": "Salience",
      "lng": "Language", "adt": "Auditory", "vsl": "Visual", "lmb": "Limbic",
      "global": "Global"}


def report(name):
    im = Image.open(FIG + "/%s.png" % name)
    dpi = im.info.get('dpi', (400, 400))[0]
    print("  %-12s %5.0f x %5.0f mm" % (name, im.size[0] / dpi * 25.4, im.size[1] / dpi * 25.4))


# ==================== S5 谐和敏感性 ====================
df = pd.read_csv(BASE + "/analysis_revision_base.csv")
hc = df[df.group == "HC"]
pdf = df[df.group == "PD"]
U = pd.read_csv(BASE + "/results_final_updrs.csv").set_index("r")
rows = []
for c in REG + ["global"]:
    a, b = np.polyfit(hc["age_at_visit"], hc[c + "_h"], 1)
    pad = pdf[c + "_h"] - (a * pdf["age_at_visit"] + b)
    d = pdf.assign(z=(pad - pad.mean()) / pad.std()).dropna(
        subset=["z", "updrs3_score", "age_at_visit", "SEX", "EDUCYRS"])
    m = smf.ols("updrs3_score ~ z + age_at_visit + C(SEX) + EDUCYRS", d).fit()
    rows.append([c, EN[c], U.loc[c, "beta"], U.loc[c, "fdr"], m.params["z"], m.pvalues["z"]])
S = pd.DataFrame(rows, columns=["r", "name", "b_lme", "fdr_lme", "b_cb", "p_cb"])
S["fdr_cb"] = multipletests(S["p_cb"], method="fdr_bh")[1]
S = S.sort_values("b_lme").reset_index(drop=True)
n_cb_sig = int((S.fdr_cb < 0.05).sum())

fig, ax = plt.subplots(figsize=(6.4, 4.0))
for i in range(len(S)):
    ax.plot([S['b_lme'][i], S['b_cb'][i]], [i, i], color=GREY, lw=1.2, zorder=1)
    ax.scatter(S['b_lme'][i], i, s=38, color=RED if S['fdr_lme'][i] < 0.05 else GREY,
               zorder=3, edgecolor='white', lw=0.6)
    ax.scatter(S['b_cb'][i], i, s=38, marker='D',
               color=BLUE if S['fdr_cb'][i] < 0.05 else LGREY,
               zorder=3, edgecolor='white', lw=0.6)
ax.axvline(0, color=GREY, ls=(0, (4, 3)), lw=0.9)
ax.set_yticks(range(len(S)))
ax.set_yticklabels(S['name'], fontsize=8.5)
ax.tick_params(axis='y', length=0, colors='black')
ax.tick_params(axis='x', colors='black', labelsize=8)
ax.set_xlabel("Standardized \u03b2  (regional brain-PAD \u2192 MDS-UPDRS III, points per SD)",
              fontsize=8.2, color='black')
# 图例与实际绘制一致
handles = [Line2D([0], [0], marker='o', color=RED, ls='', label='Mixed model, scanner random effect (FDR < 0.05)'),
           Line2D([0], [0], marker='o', color=GREY, ls='', label='Mixed model (FDR \u2265 0.05)'),
           Line2D([0], [0], marker='D', color=BLUE if n_cb_sig else LGREY, ls='',
                  label='ComBat harmonized' + ('' if n_cb_sig else ' (none reached FDR < 0.05)'))]
ax.legend(handles=handles, loc='lower right', fontsize=7, frameon=False, labelcolor='black')
ax.set_title("Sensitivity of the motor associations to scanner harmonization",
             fontsize=9.5, loc='left', fontweight='bold', color='black')
for f in ['png', 'pdf']:
    fig.savefig(FIG + "/Figure_S5." + f, dpi=400 if f == 'png' else None, bbox_inches='tight')
plt.close(fig)
print("S5 重绘 (图例与实绘一致; ComBat 达 FDR<0.05 的网络数 = %d):" % n_cb_sig)
report("Figure_S5")

# ==================== S7 去掉无对应数据的图例项 ====================
Ud = pd.read_csv(BASE + "/rev_long_updrs3_score.csv")
Md = pd.read_csv(BASE + "/rev_long_moca.csv")
nU, nM = int(Ud["n_subj"].max()), int(Md["n_subj"].max())
n_nom = int((Ud["p"] < 0.05).sum() + (Md["p"] < 0.05).sum())

fig, axes = plt.subplots(1, 2, figsize=(6.6, 3.4),
                         gridspec_kw={'wspace': 0.10, 'left': 0.20, 'right': 0.98,
                                      'top': 0.78, 'bottom': 0.21})
for ax, (D, ttl, xl) in zip(axes, [
        (Ud, "A  Motor progression (MDS-UPDRS III, n = %d)" % nU,
         "\u03b2  time \u00d7 brain-PAD  (points per SD per year)"),
        (Md, "B  Cognitive progression (MoCA, n = %d)" % nM,
         "\u03b2  time \u00d7 brain-PAD  (points per SD per year)")]):
    D = D.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(D))
    for i in range(len(D)):
        col = ORANGE if D['p'][i] < 0.05 else GREY
        ax.plot([D['lo'][i], D['hi'][i]], [i, i], color=col, lw=1.6, solid_capstyle='round')
        ax.scatter(D['beta_time_x_PAD'][i], i, s=24, color=col, zorder=3,
                   edgecolor='white', lw=0.5)
    ax.axvline(0, color='#444', ls=(0, (4, 3)), lw=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels(D['Network'], fontsize=7.6)
    ax.set_xlabel(xl, fontsize=7.4, color='black')
    ax.set_title(ttl, fontsize=8.4, loc='left', fontweight='bold', pad=6, color='black')
    ax.tick_params(colors='black', labelsize=7.4)
axes[1].set_yticklabels([])
for a in axes:
    a.tick_params(axis='y', length=0)
fig.suptitle("Baseline regional brain age does not predict the rate of clinical change",
             fontsize=9.6, fontweight='bold', x=0.02, ha='left', y=0.97, color='black')
note = ("Points are time-by-brain-PAD interaction coefficients with 95% confidence intervals. "
        "No network reached nominal significance for either outcome."
        if n_nom == 0 else
        "Orange marks nominal P < 0.05, none significant after correction.")
fig.text(0.02, 0.895, note, fontsize=7, color='#444', style='italic', va='top')
for f in ['png', 'pdf']:
    fig.savefig(FIG + "/Figure_S7." + f, dpi=400 if f == 'png' else None)
plt.close(fig)
print("\nS7 重绘 (名义显著网络数 = %d, 已去掉无对应数据的图例项):" % n_nom)
report("Figure_S7")
