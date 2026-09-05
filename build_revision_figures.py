# -*- coding: utf-8 -*-
"""主图生成: 重建超出 170x210 mm 排版区的图, 并新增补充图。
Figure 4  重建(内容更新: 加入 GWR 三阶段分析 + 采集方差分解), 2x2 版式
Figure 5  重排为 2 行(原 1x3 宽 260 mm)
Figure S5 重排(原 211 mm)
Figure 2 / S6 保持内容, 通过 dpi 元数据缩放至排版区内
其余合规图件直接复制
"""
import os, shutil
os.environ["MPLCONFIGDIR"] = r"G:/Temp/mpl"
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import statsmodels.formula.api as smf
from PIL import Image
import warnings
warnings.filterwarnings("ignore")
Image.MAX_IMAGE_PIXELS = None

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False,
                     'axes.spines.right': False, 'axes.linewidth': 0.9})

BASE = r"F:/PD_brainage_data"
SRC = BASE + "/figures_previous"
OUT = BASE + "/SUBMISSION_BrainComms_R1/Figures"
os.makedirs(OUT, exist_ok=True)
BLUE, ORANGE, GREY, GREEN = '#0F4D92', '#D55E00', '#9A9A9A', '#009E73'
MAXW, MAXH = 170.0, 210.0   # mm, Brain Communications 排版区


def save(fig, name, tight=True):
    kw = {'bbox_inches': 'tight'} if tight else {}
    for f in ['png', 'pdf']:
        fig.savefig("%s/%s.%s" % (OUT, name, f), dpi=400 if f == 'png' else None, **kw)
    plt.close(fig)
    im = Image.open("%s/%s.png" % (OUT, name))
    dpi = im.info.get('dpi', (400, 400))[0]
    print("  %-14s %5.0f x %5.0f mm" % (name, im.size[0] / dpi * 25.4, im.size[1] / dpi * 25.4))


def star(p):
    return '***' if p < .001 else '**' if p < .01 else '*' if p < .05 else 'n.s.'


# ---------------- 数据 ----------------
df = pd.read_csv(BASE + "/analysis_revision_base.csv").merge(
    pd.read_csv(BASE + "/network_thickness_gwr.csv"), on="PATNO", how="left")
REG = ["sns", "fpn", "drs", "vnt", "dft", "slt", "lng", "adt", "vsl", "lmb"]


def combat(vals, batch, mod):
    y = vals.astype(float).copy()
    ok = np.isfinite(y)
    X = np.column_stack([np.ones(ok.sum())] + [m[ok] for m in mod])
    beta = np.linalg.lstsq(X, y[ok], rcond=None)[0]
    fit = np.full(len(y), np.nan)
    fit[ok] = X @ beta
    res = y - fit
    sd = np.nanstd(res[ok])
    out = res / sd
    for b in pd.unique(batch):
        m = (batch == b) & ok
        if m.sum() >= 2:
            out[m] = (out[m] - np.nanmean(out[m])) / (np.nanstd(out[m]) or 1.0)
    return out * sd + fit


mod = [df["age_at_visit"].values, (df["SEX"] == df["SEX"].iloc[0]).astype(float).values,
       (df["group"] == "PD").astype(float).values]
for r in REG:
    df[r + "_gwr_h"] = combat(df[r + "_gwr"].values, df["batch"].values, mod)
pdf = df[df.group == "PD"].copy()
for r in REG:
    for s in ["_thick", "_gwr_h"]:
        v = pdf[r + s]
        pdf[r + s + "z"] = (v - v.mean()) / v.std()


def lme(outcome, terms):
    d = pdf.dropna(subset=[outcome, "age_at_visit", "SEX", "EDUCYRS", "batch"] + terms).copy()
    m = smf.mixedlm("%s ~ %s + age_at_visit + C(SEX) + EDUCYRS" % (outcome, " + ".join(terms)),
                    d, groups=d["batch"]).fit()
    k = terms[0]
    ci = m.conf_int()
    return m.params[k], ci.loc[k, 0], ci.loc[k, 1], m.pvalues[k], len(d)


# ==================== Figure 4 (重建) ====================
print("Figure 4 (重建, 2x2):")
ITEMS = {"updrs3_score": [("vsl", "Visual"), ("vnt", "Ventral attention"),
                          ("slt", "Salience"), ("lng", "Language")],
         "moca": [("fpn", "Frontoparietal"), ("vnt", "Ventral attention")]}

fig = plt.figure(figsize=(6.6, 7.4))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], hspace=0.45, wspace=0.42,
                      top=0.885, bottom=0.135, left=0.15, right=0.975)

for j, (out, ttl, xl) in enumerate([
        ("updrs3_score", "A  Motor severity", "MDS-UPDRS III  \u03b2  (points per SD)"),
        ("moca", "B  Global cognition", "MoCA  \u03b2  (points per SD)")]):
    ax = fig.add_subplot(gs[0, j])
    items = ITEMS[out]
    yp = np.arange(len(items))[::-1]
    for yi, (r, name) in zip(yp, items):
        for off, col, mk, terms in [
                (0.24, BLUE, 'o', [r + "_z"]),
                (0.00, ORANGE, 'D', [r + "_z", r + "_thickz"]),
                (-0.24, GREEN, 's', [r + "_z", r + "_thickz", r + "_gwr_hz"])]:
            b, lo, hi, p, n = lme(out, terms)
            ax.plot([lo, hi], [yi + off] * 2, color=col, lw=1.6, solid_capstyle='round')
            ax.scatter(b, yi + off, s=30, color=col, marker=mk, zorder=3,
                       edgecolor='white', lw=0.5)
    ax.axvline(0, color=GREY, ls=(0, (4, 3)), lw=0.9)
    ax.set_yticks(yp)
    ax.set_yticklabels([n for _, n in items], fontsize=8)
    ax.set_ylim(-0.65, len(items) - 0.35)
    ax.set_xlabel(xl, fontsize=8)
    ax.set_title(ttl, fontsize=9.5, loc='left', fontweight='bold', color='black')
    ax.tick_params(colors='black', labelsize=8)
    ax.xaxis.label.set_color('black')

fig.legend(handles=[
    Line2D([0], [0], color=BLUE, marker='o', lw=1.6, label='brain-PAD alone'),
    Line2D([0], [0], color=ORANGE, marker='D', lw=1.6, label='+ cortical thickness'),
    Line2D([0], [0], color=GREEN, marker='s', lw=1.6, label='+ thickness and harmonized GWR')],
    loc='upper center', bbox_to_anchor=(0.5, 0.935), ncol=3, fontsize=7.2,
    frameon=False, handletextpad=0.35, columnspacing=1.1, labelcolor='black')

# C 采集方差分解
axC = fig.add_subplot(gs[1, 0])
V = pd.read_csv(BASE + "/rev_B2c_gwr_scanner_variance.csv")
x = np.arange(len(V))
w = 0.27
axC.bar(x - w, V["GWR_scanner_var_%"], w, color=ORANGE, label='GWR (input)')
axC.bar(x, V["Thickness_scanner_var_%"], w, color=GREY, label='Thickness (input)')
axC.bar(x + w, V["brainPAD_scanner_var_%"], w, color=BLUE, label='brain-PAD (output)')
axC.set_xticks(x)
axC.set_xticklabels(V["Network"], rotation=55, ha='right', fontsize=6.6)
axC.set_ylabel('Scanner-attributable\nvariance (%)', fontsize=8)
axC.set_title('C  Acquisition robustness', fontsize=9.5, loc='left',
              fontweight='bold', color='black')
axC.legend(fontsize=6.6, frameon=False, loc='upper left', labelcolor='black')
axC.tick_params(colors='black', labelsize=7.5)
axC.yaxis.label.set_color('black')

# D GWR-运动关联 谐和前后
axD = fig.add_subplot(gs[1, 1])
G = pd.read_csv(BASE + "/rev_B2c_updrs3_score_gwr_site.csv")
yp = np.arange(len(G))[::-1]
for yi, (_, r) in zip(yp, G.iterrows()):
    axD.scatter(r["b_GWR_LME"], yi + 0.18, s=26, color=ORANGE, zorder=3,
                edgecolor='white', lw=0.5)
    axD.scatter(r["b_GWRh_LME"], yi - 0.18, s=26, color=BLUE, marker='D', zorder=3,
                edgecolor='white', lw=0.5)
    axD.plot([r["b_GWR_LME"], r["b_GWRh_LME"]], [yi + 0.18, yi - 0.18],
             color='#CCC', lw=0.7, zorder=1)
axD.axvline(0, color=GREY, ls=(0, (4, 3)), lw=0.9)
axD.set_yticks(yp)
axD.set_yticklabels(G["Network"], fontsize=6.6)
axD.set_xlabel('MDS-UPDRS III  \u03b2  (points per SD of GWR)', fontsize=7.5)
axD.set_title('D  GWR effect is acquisition-driven', fontsize=9.2, loc='left',
              fontweight='bold', color='black')
axD.legend(handles=[
    Line2D([0], [0], color=ORANGE, marker='o', ls='', label='raw GWR (9/10 FDR < 0.05)'),
    Line2D([0], [0], color=BLUE, marker='D', ls='', label='harmonized GWR (0/10)')],
    fontsize=6.4, frameon=False, loc='upper center', bbox_to_anchor=(0.45, -0.20),
    labelcolor='black')
axD.tick_params(colors='black', labelsize=7.5)
axD.xaxis.label.set_color('black')

fig.suptitle("Regional brain age carries clinical information beyond both model inputs",
             fontsize=10.5, fontweight='bold', x=0.02, ha='left', y=0.975, color='black')
save(fig, "Figure_4", tight=False)

# ==================== Figure 5 (重排为 2 行) ====================
print("Figure 5 (重排 2 行):")
A = pd.read_csv(BASE + "/DAT_coupling.csv")
B = pd.read_csv(BASE + "/DAT_incremental_validity.csv")
V5 = pd.read_csv(BASE + "/DAT_variance_partition.csv")

fig = plt.figure(figsize=(6.6, 7.2))
gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 0.9], hspace=0.5, wspace=0.55)

ax = fig.add_subplot(gs[0, 0])
Ar = A.iloc[::-1].reset_index(drop=True)
for i in range(len(Ar)):
    ax.plot([Ar['lo'][i], Ar['hi'][i]], [i, i], color=GREY, lw=1.7, solid_capstyle='round')
    ax.scatter(Ar['beta'][i], i, s=24, color=GREY, zorder=3, edgecolor='white', lw=0.5)
ax.axvline(0, color='black', ls=(0, (4, 3)), lw=0.9)
ax.set_yticks(np.arange(len(Ar)))
ax.set_yticklabels(Ar['name'], fontsize=7.2)
ax.tick_params(axis='y', length=0, colors='black')
ax.tick_params(axis='x', labelsize=7.5, colors='black')
ax.set_xlabel('\u03b2  brain-PAD (years per unit putamen SBR)', fontsize=7.5, color='black')
ax.set_title('A  Decoupled from dopaminergic\n    denervation', fontsize=9.5, loc='left',
             fontweight='bold', color='black')
ax.text(0.98, 0.03, 'all FDR > 0.47', transform=ax.transAxes, ha='right',
        fontsize=7, color='black', style='italic')

ax = fig.add_subplot(gs[0, 1])
OCN = {'updrs3_score': 'UPDRS-III', 'moca': 'MoCA', 'UPDRS-III': 'UPDRS-III', 'MoCA': 'MoCA'}
B['ocn'] = B['outcome'].map(OCN)
yp = np.arange(len(B))[::-1]
for yi, (_, r) in zip(yp, B.iterrows()):
    ax.scatter(abs(r['b_PADonly']), yi + 0.16, s=32, color=BLUE, zorder=3,
               edgecolor='white', lw=0.5)
    ax.scatter(abs(r['b_PADadjDAT']), yi - 0.16, s=32, color=ORANGE, marker='D',
               zorder=3, edgecolor='white', lw=0.5)
    ax.plot([abs(r['b_PADonly']), abs(r['b_PADadjDAT'])], [yi + 0.16, yi - 0.16],
            color='#CCC', lw=0.7, zorder=1)
    if r['p_PADadjDAT'] < 0.05:
        ax.text(1.02, yi, '*', transform=ax.get_yaxis_transform(), va='center',
                fontsize=10, color='black')
ax.set_yticks(yp)
ax.set_yticklabels(["%s\n(%s)" % (r['network'], r['ocn']) for _, r in B.iterrows()], fontsize=6.6)
ax.tick_params(axis='y', length=0, colors='black')
ax.tick_params(axis='x', labelsize=7.5, colors='black')
ax.set_xlabel('|\u03b2|  clinical points per SD of brain-PAD', fontsize=7.5, color='black')
ax.set_title('B  Association survives adjustment\n    for DAT binding', fontsize=9.5,
             loc='left', fontweight='bold', color='black')
ax.legend(handles=[Line2D([0], [0], marker='o', color=BLUE, ls='', label='brain-PAD alone'),
                   Line2D([0], [0], marker='D', color=ORANGE, ls='', label='brain-PAD | DAT')],
          loc='lower right', fontsize=6.6, frameon=False, labelcolor='black')

ax = fig.add_subplot(gs[1, :])
grp = V5.groupby('outcome')[['unique_PAD', 'unique_SBR', 'shared']].mean() * 100
grp = grp.reindex(['UPDRS-III', 'MoCA'])
xs = np.arange(len(grp))
ax.bar(xs, grp['unique_PAD'], 0.45, color=BLUE, label='Brain age (unique)')
ax.bar(xs, grp['unique_SBR'], 0.45, bottom=grp['unique_PAD'], color=ORANGE,
       label='DAT-SPECT (unique)')
ax.bar(xs, grp['shared'].clip(lower=0), 0.45,
       bottom=grp['unique_PAD'] + grp['unique_SBR'], color=GREY, label='Shared')
ax.set_xticks(xs)
ax.set_xticklabels(grp.index, fontsize=8.5)
ax.set_ylabel('Unique variance explained\n(\u0394R\u00b2, percentage points)', fontsize=8)
ax.set_title('C  Variance partition', fontsize=9.5, loc='left', fontweight='bold', color='black')
ax.legend(fontsize=7, frameon=False, loc='upper right', labelcolor='black')
for xi, (pa, sb) in enumerate(zip(grp['unique_PAD'], grp['unique_SBR'])):
    ax.text(xi, pa / 2, "%.2f" % pa, ha='center', va='center', fontsize=7.5,
            color='white', fontweight='bold')
    if sb > 0.15:
        ax.text(xi, pa + sb / 2, "%.2f" % sb, ha='center', va='center', fontsize=7.5,
                color='white', fontweight='bold')
    else:
        ax.text(xi, pa + sb + 0.08, "%.2f" % sb, ha='center', va='bottom', fontsize=7.5,
                color=ORANGE)
ax.tick_params(colors='black', labelsize=8)
ax.yaxis.label.set_color('black')

fig.suptitle("Regional cortical brain age is independent of dopaminergic denervation",
             fontsize=10.5, fontweight='bold', x=0.02, ha='left', y=0.99, color='black')
save(fig, "Figure_5")

# ==================== 复制合规图件 / 缩放超限图件 ====================
print("\n其余图件:")


def rescale(name):
    """保持像素不变, 调整 dpi 元数据使物理尺寸落入排版区。"""
    im = Image.open("%s/%s.png" % (SRC, name))
    dpi = im.info.get('dpi', (400, 400))[0]
    w, h = im.size[0] / dpi * 25.4, im.size[1] / dpi * 25.4
    k = max(w / MAXW, h / MAXH, 1.0)
    nd = dpi * k
    im.save("%s/%s.png" % (OUT, name), dpi=(nd, nd))
    if os.path.exists("%s/%s.pdf" % (SRC, name)):
        shutil.copy("%s/%s.pdf" % (SRC, name), "%s/%s.pdf" % (OUT, name))
    print("  %-14s %5.0f x %5.0f mm  (dpi %.0f -> %.0f)"
          % (name, w / k, h / k, dpi, nd))


for n in ["Figure_1", "Figure_3", "Figure_S1", "Figure_S2", "Figure_S3", "Figure_S4"]:
    for ext in ["png", "pdf"]:
        p = "%s/%s.%s" % (SRC, n, ext)
        if os.path.exists(p):
            shutil.copy(p, "%s/%s.%s" % (OUT, n, ext))
    im = Image.open("%s/%s.png" % (OUT, n))
    d = im.info.get('dpi', (400, 400))[0]
    print("  %-14s %5.0f x %5.0f mm  (原样复制)" % (n, im.size[0] / d * 25.4, im.size[1] / d * 25.4))

for n in ["Figure_2", "Figure_S5", "Figure_S6"]:
    rescale(n)
print("\n完成 ->", OUT)
