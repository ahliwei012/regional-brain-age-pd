# -*- coding: utf-8 -*-
"""尺度稳健性检验: 各网络预测年龄的动态范围差异极大(尤以 ventral attention 被压缩),
需确认主关联不是由 brain-PAD 的尺度/分布特性驱动。
做三件事: (1) 秩变换 brain-PAD 重跑; (2) Spearman 偏相关; (3) 稳健(Huber)回归。
"""
import pandas as pd, numpy as np, warnings
import statsmodels.formula.api as smf
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from scipy import stats
warnings.filterwarnings("ignore")

BASE = r"F:/PD_brainage_data"
REG = ["sns", "fpn", "drs", "vnt", "dft", "slt", "lng", "adt", "vsl", "lmb"]
EN = {"sns": "Sensorimotor", "fpn": "Frontoparietal", "drs": "Dorsal attention",
      "vnt": "Ventral attention", "dft": "Default mode", "slt": "Salience",
      "lng": "Language", "adt": "Auditory", "vsl": "Visual", "lmb": "Limbic",
      "global": "Global"}

df = pd.read_csv(BASE + "/analysis_revision_base.csv")
pdf = df[df.group == "PD"].copy()

# 秩变换(逆正态), 完全去除尺度与分布形状的影响
for c in REG + ["global"]:
    v = pdf[c + "_PAD"]
    pdf[c + "_rank"] = stats.norm.ppf((v.rank() - 0.5) / v.notna().sum())


def run(outcome, suf):
    rows = []
    for c in REG + ["global"]:
        d = pdf.dropna(subset=[c + suf, outcome, "age_at_visit", "SEX",
                               "EDUCYRS", "batch"]).rename(columns={c + suf: "x"})
        m = smf.mixedlm("%s ~ x + age_at_visit + C(SEX) + EDUCYRS" % outcome,
                        d, groups=d["batch"]).fit()
        # 稳健回归(Huber), 对离群值不敏感
        rb = smf.rlm("%s ~ x + age_at_visit + C(SEX) + EDUCYRS" % outcome,
                     d, M=sm.robust.norms.HuberT()).fit()
        # 偏 Spearman: 先各自对协变量取残差再求秩相关
        ry = smf.ols("%s ~ age_at_visit + C(SEX) + EDUCYRS" % outcome, d).fit().resid
        rx = smf.ols("x ~ age_at_visit + C(SEX) + EDUCYRS", d).fit().resid
        rho, prho = stats.spearmanr(rx, ry)
        rows.append([EN[c], m.params["x"], m.pvalues["x"], rb.params["x"],
                     rb.pvalues["x"], rho, prho, len(d)])
    r = pd.DataFrame(rows, columns=["Network", "beta_LME", "P_LME", "beta_robust",
                                    "P_robust", "partial_rho", "P_spearman", "n"])
    for col in ["P_LME", "P_robust", "P_spearman"]:
        r["fdr_" + col.split("_")[1]] = multipletests(r[col], method="fdr_bh")[1]
    return r


for outcome, tag in [("updrs3_score", "UPDRS-III"), ("moca", "MoCA")]:
    print("=" * 78)
    print("%s  --  秩变换 brain-PAD (逆正态), 尺度无关" % tag)
    print("=" * 78)
    r = run(outcome, "_rank")
    r.to_csv(BASE + "/rev_A2_%s_rank.csv" % outcome, index=False, encoding="utf-8-sig")
    print(r.round(4).to_string(index=False))
    print("  FDR<0.05 (LME): %s" % r[r.fdr_LME < 0.05]["Network"].tolist())
    print("  FDR<0.05 (Spearman偏相关): %s\n" % r[r.fdr_spearman < 0.05]["Network"].tolist())

# 参考: 原始 z 标准化结果
print("=" * 78)
print("对照 -- 原始 z 标准化 brain-PAD")
print("=" * 78)
for outcome, tag in [("updrs3_score", "UPDRS-III"), ("moca", "MoCA")]:
    r = run(outcome, "_z")
    print("[%s] FDR<0.05: %s" % (tag, r[r.fdr_LME < 0.05]["Network"].tolist()))
