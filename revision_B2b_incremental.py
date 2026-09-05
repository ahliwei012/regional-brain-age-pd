# -*- coding: utf-8 -*-
"""B2b 增量效度与方差分解: brain-PAD vs 皮层厚度 vs GWR。
起因: 控制同网络平均 GWR 后, 运动关联由 FDR 0.016 降至 0.085。
需判定这是"过度校正"(GWR 本就是脑龄模型的输入)还是 brain-PAD 缺乏独立信息。
方法与稿件中 brain-PAD vs DAT 的方差分解一致: 嵌套模型 F 检验 + 独有 R^2。
"""
import numpy as np, pandas as pd, warnings
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings("ignore")

BASE = r"F:/PD_brainage_data"
REG = ["sns", "fpn", "drs", "vnt", "dft", "slt", "lng", "adt", "vsl", "lmb"]
EN = {"sns": "Sensorimotor", "fpn": "Frontoparietal", "drs": "Dorsal attention",
      "vnt": "Ventral attention", "dft": "Default mode", "slt": "Salience",
      "lng": "Language", "adt": "Auditory", "vsl": "Visual", "lmb": "Limbic"}

df = pd.read_csv(BASE + "/analysis_revision_base.csv").merge(
    pd.read_csv(BASE + "/network_thickness_gwr.csv"), on="PATNO", how="left")
pdf = df[df.group == "PD"].copy()

# 厚度与 GWR 一律 z 标准化, 使系数与 brain-PAD 可直接比较
for r in REG:
    for s in ["_thick", "_gwr"]:
        v = pdf[r + s]
        pdf[r + s + "z"] = (v - v.mean()) / v.std()

COV = "age_at_visit + C(SEX) + EDUCYRS"


def r2_adj(m):
    return 1 - (1 - m.rsquared) * (m.nobs - 1) / (m.nobs - m.df_model - 1)


def analyse(outcome, tag):
    rows = []
    for r in REG:
        d = pdf.dropna(subset=[r + "_z", r + "_thickz", r + "_gwrz", outcome,
                               "age_at_visit", "SEX", "EDUCYRS"]).copy()
        f0 = smf.ols("%s ~ %s" % (outcome, COV), d).fit()
        fp = smf.ols("%s ~ %s_z + %s" % (outcome, r, COV), d).fit()
        fg = smf.ols("%s ~ %s_gwrz + %s" % (outcome, r, COV), d).fit()
        ft = smf.ols("%s ~ %s_thickz + %s" % (outcome, r, COV), d).fit()
        fall = smf.ols("%s ~ %s_z + %s_thickz + %s_gwrz + %s" % (outcome, r, r, r, COV), d).fit()
        # 独有贡献 = 全模型 R^2 - 去掉该项后的 R^2
        f_nopad = smf.ols("%s ~ %s_thickz + %s_gwrz + %s" % (outcome, r, r, COV), d).fit()
        f_nogwr = smf.ols("%s ~ %s_z + %s_thickz + %s" % (outcome, r, r, COV), d).fit()
        f_noth = smf.ols("%s ~ %s_z + %s_gwrz + %s" % (outcome, r, r, COV), d).fit()
        # 嵌套 F 检验: 在 厚度+GWR 之上加入 brain-PAD 是否显著改善拟合
        Fp = fall.compare_f_test(f_nopad)
        Fg = fall.compare_f_test(f_nogwr)
        rows.append([EN[r],
                     fp.params[r + "_z"], fp.pvalues[r + "_z"],
                     fg.params[r + "_gwrz"], fg.pvalues[r + "_gwrz"],
                     ft.params[r + "_thickz"], ft.pvalues[r + "_thickz"],
                     (fp.rsquared - f0.rsquared) * 100,
                     (fg.rsquared - f0.rsquared) * 100,
                     (ft.rsquared - f0.rsquared) * 100,
                     (fall.rsquared - f_nopad.rsquared) * 100, Fp[1],
                     (fall.rsquared - f_nogwr.rsquared) * 100, Fg[1],
                     int(fall.nobs)])
    out = pd.DataFrame(rows, columns=[
        "Network", "b_PAD_alone", "P_PAD_alone", "b_GWR_alone", "P_GWR_alone",
        "b_thick_alone", "P_thick_alone", "R2_PAD_%", "R2_GWR_%", "R2_thick_%",
        "R2_PAD_unique_%", "P_F_addPAD", "R2_GWR_unique_%", "P_F_addGWR", "n"])
    for c in ["P_PAD_alone", "P_GWR_alone", "P_thick_alone", "P_F_addPAD", "P_F_addGWR"]:
        out["fdr_" + c[2:]] = multipletests(out[c], method="fdr_bh")[1]
    out.to_csv(BASE + "/rev_B2b_%s_incremental.csv" % outcome, index=False, encoding="utf-8-sig")

    print("=" * 92)
    print("%s  --  单独效应 (标准化 β, 每 1 SD)" % tag)
    print("=" * 92)
    print(out[["Network", "b_PAD_alone", "P_PAD_alone", "b_GWR_alone", "P_GWR_alone",
               "b_thick_alone", "P_thick_alone"]].round(4).to_string(index=False))
    print("\n%s  --  解释方差 (%% , 在协变量之上)" % tag)
    print(out[["Network", "R2_PAD_%", "R2_GWR_%", "R2_thick_%",
               "R2_PAD_unique_%", "P_F_addPAD", "R2_GWR_unique_%", "P_F_addGWR"]]
          .round(3).to_string(index=False))
    print("\n  在 厚度+GWR 之上加入 brain-PAD 仍显著改善拟合 (FDR<0.05): %s"
          % out[out["fdr_F_addPAD"] < .05]["Network"].tolist())
    print("  在 厚度+brain-PAD 之上加入 GWR 仍显著改善拟合 (FDR<0.05): %s\n"
          % out[out["fdr_F_addGWR"] < .05]["Network"].tolist())
    return out


for outcome, tag in [("updrs3_score", "UPDRS-III"), ("moca", "MoCA")]:
    analyse(outcome, tag)
