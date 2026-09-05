# -*- coding: utf-8 -*-
"""Longitudinal analyses.
基于 PPMI curated 纵向表构建随访数据集; 满足 >=3 次随访的 PD 为 349 例(UPDRS-III) / 440 例(MoCA)。
问题: 基线区域 brain-PAD 是否预测随后的临床变化速率? (time x brain-PAD 交互)
另加: H&Y 进展 与 MCI 转化 的生存/logistic 分析。
"""
import numpy as np, pandas as pd, warnings
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings("ignore")

BASE = r"F:/PD_brainage_data"
REG = ["sns", "fpn", "drs", "vnt", "dft", "slt", "lng", "adt", "vsl", "lmb"]
EN = {"sns": "Sensorimotor", "fpn": "Frontoparietal", "drs": "Dorsal attention",
      "vnt": "Ventral attention", "dft": "Default mode", "slt": "Salience",
      "lng": "Language", "adt": "Auditory", "vsl": "Visual", "lmb": "Limbic",
      "global": "Global"}

base = pd.read_csv(BASE + "/analysis_revision_base.csv")
base = base[base.group == "PD"]
lg = pd.read_csv(BASE + "/master_clinical_longitudinal.csv")

# 访视 -> 随访年数
YEAR = {"BL": 0.0, "V01": 0.25, "V02": 0.5, "V03": 0.75, "V04": 1.0, "V05": 1.5,
        "V06": 2.0, "V07": 2.5, "V08": 3.0, "V09": 3.5, "V10": 4.0, "V11": 5.0,
        "V12": 6.0, "V13": 7.0, "V14": 8.0, "V15": 9.0, "V16": 10.0, "V17": 11.0,
        "V18": 12.0, "V19": 13.0, "V20": 14.0}
lg["YEAR"] = lg["EVENT_ID"].map(YEAR)
lg = lg.dropna(subset=["YEAR"])

cols = ["PATNO", "SEX", "EDUCYRS", "batch"] + [r + "_z" for r in REG] + ["global_z"]
bl_age = lg[lg.EVENT_ID == "BL"][["PATNO", "age_at_visit"]].rename(
    columns={"age_at_visit": "age_bl"})
d = lg.merge(base[cols], on="PATNO", how="inner").merge(bl_age, on="PATNO", how="left")


def star(p):
    return "***" if p < .001 else "**" if p < .01 else "*" if p < .05 else ""


def progression(outcome, tag):
    dd = d.dropna(subset=[outcome, "YEAR", "age_bl", "SEX", "EDUCYRS"]).copy()
    vc = dd.groupby("PATNO").size()
    dd = dd[dd.PATNO.isin(vc[vc >= 3].index)]
    fu = dd.groupby("PATNO").YEAR.max()
    print("\n%s: %d 例, %d 次访视, 中位随访 %.1f 年 (IQR %.1f-%.1f)"
          % (tag, dd.PATNO.nunique(), len(dd), fu.median(), fu.quantile(.25), fu.quantile(.75)))
    rows = []
    for r in REG + ["global"]:
        x = dd.rename(columns={r + "_z": "PAD"}).dropna(subset=["PAD"])
        m = smf.mixedlm("%s ~ YEAR*PAD + age_bl + C(SEX) + EDUCYRS" % outcome, x,
                        groups=x["PATNO"], re_formula="~YEAR").fit(method="lbfgs")
        t = "YEAR:PAD"
        ci = m.conf_int().loc[t]
        rows.append([EN[r], m.params[t], ci[0], ci[1], m.pvalues[t],
                     m.params["YEAR"], x.PATNO.nunique()])
    o = pd.DataFrame(rows, columns=["Network", "beta_time_x_PAD", "lo", "hi", "p",
                                    "beta_time", "n_subj"])
    o["fdr"] = multipletests(o["p"], method="fdr_bh")[1]
    o["sig"] = o["fdr"].map(star)
    o.to_csv(BASE + "/rev_long_%s.csv" % outcome, index=False, encoding="utf-8-sig")
    print(o.round(4).to_string(index=False))
    print("  年均变化率 (time 主效应): %.3f 分/年" % o["beta_time"].mean())
    print("  FDR<0.05: %s" % (o[o.fdr < .05]["Network"].tolist() or "无"))
    return o


print("=" * 84)
print("基线区域脑龄 -> 随后临床变化速率 (time x brain-PAD 交互)")
print("=" * 84)
progression("updrs3_score", "MDS-UPDRS-III")
progression("moca", "MoCA")
progression("updrs2_score", "MDS-UPDRS-II (日常生活运动体验)")

# ---- H&Y 进展: 首次达到 H&Y >=3 ----
print("\n" + "=" * 84)
print("基线区域脑龄 -> 疾病分期进展 (首次达到 Hoehn & Yahr >= 3)")
print("=" * 84)
hy = d.dropna(subset=["NHY"]).copy()
ev = hy[hy.NHY >= 3].groupby("PATNO").YEAR.min().rename("t_event")
last = hy.groupby("PATNO").YEAR.max().rename("t_last")
surv = pd.concat([ev, last], axis=1)
surv["event"] = surv["t_event"].notna().astype(int)
surv["time"] = surv["t_event"].fillna(surv["t_last"])
surv = surv[surv["time"] > 0].reset_index().merge(
    base[cols + ["age_at_visit"]], on="PATNO", how="inner")
print("  n=%d, 事件数=%d (%.1f%%), 中位随访 %.1f 年"
      % (len(surv), surv.event.sum(), 100 * surv.event.mean(), surv.time.median()))
try:
    from lifelines import CoxPHFitter
    rows = []
    for r in REG + ["global"]:
        x = surv.dropna(subset=[r + "_z", "EDUCYRS"])[
            ["time", "event", r + "_z", "age_at_visit", "EDUCYRS"]].copy()
        x["male"] = (surv.loc[x.index, "SEX"] == surv["SEX"].iloc[0]).astype(int)
        cph = CoxPHFitter().fit(x, "time", "event")
        rows.append([EN[r], np.exp(cph.params_[r + "_z"]),
                     cph.summary.loc[r + "_z", "p"], len(x)])
    o = pd.DataFrame(rows, columns=["Network", "HR_per_SD", "p", "n"])
    o["fdr"] = multipletests(o["p"], method="fdr_bh")[1]
    o["sig"] = o["fdr"].map(star)
    o.to_csv(BASE + "/rev_long_HY3_cox.csv", index=False, encoding="utf-8-sig")
    print(o.round(4).to_string(index=False))
    print("  FDR<0.05: %s" % (o[o.fdr < .05]["Network"].tolist() or "无"))
except ImportError:
    print("  [跳过 Cox] lifelines 未安装, 改用 logistic")
    rows = []
    for r in REG + ["global"]:
        x = surv.dropna(subset=[r + "_z", "EDUCYRS"]).rename(columns={r + "_z": "PAD"})
        m = smf.logit("event ~ PAD + age_at_visit + C(SEX) + EDUCYRS + time", x).fit(disp=0)
        rows.append([EN[r], np.exp(m.params["PAD"]), m.pvalues["PAD"], len(x)])
    o = pd.DataFrame(rows, columns=["Network", "OR_per_SD", "p", "n"])
    o["fdr"] = multipletests(o["p"], method="fdr_bh")[1]
    o["sig"] = o["fdr"].map(star)
    o.to_csv(BASE + "/rev_long_HY3_logit.csv", index=False, encoding="utf-8-sig")
    print(o.round(4).to_string(index=False))
    print("  FDR<0.05: %s" % (o[o.fdr < .05]["Network"].tolist() or "无"))

# ---- MCI 转化: 基线认知正常者随访中首次 cogstate>1 ----
print("\n" + "=" * 84)
print("基线区域脑龄 -> MCI/痴呆转化 (基线认知正常者)")
print("=" * 84)
cg = d.dropna(subset=["cogstate"]).copy()
bl_norm = cg[(cg.EVENT_ID == "BL") & (cg.cogstate == 1)]["PATNO"].unique()
cg = cg[cg.PATNO.isin(bl_norm)]
conv = cg[cg.cogstate > 1].groupby("PATNO").YEAR.min().rename("t_event")
lastc = cg.groupby("PATNO").YEAR.max().rename("t_last")
s2 = pd.concat([conv, lastc], axis=1)
s2["event"] = s2["t_event"].notna().astype(int)
s2["time"] = s2["t_event"].fillna(s2["t_last"])
s2 = s2[s2["time"] > 0].reset_index().merge(base[cols + ["age_at_visit"]], on="PATNO", how="inner")
print("  基线认知正常且有随访: n=%d, 转化 %d 例 (%.1f%%), 中位随访 %.1f 年"
      % (len(s2), s2.event.sum(), 100 * s2.event.mean(), s2.time.median()))
rows = []
for r in REG + ["global"]:
    x = s2.dropna(subset=[r + "_z", "EDUCYRS"]).rename(columns={r + "_z": "PAD"})
    m = smf.logit("event ~ PAD + age_at_visit + C(SEX) + EDUCYRS + time", x).fit(disp=0)
    rows.append([EN[r], np.exp(m.params["PAD"]), m.pvalues["PAD"], len(x)])
o = pd.DataFrame(rows, columns=["Network", "OR_per_SD", "p", "n"])
o["fdr"] = multipletests(o["p"], method="fdr_bh")[1]
o["sig"] = o["fdr"].map(star)
o.to_csv(BASE + "/rev_long_MCI_conversion.csv", index=False, encoding="utf-8-sig")
print(o.round(4).to_string(index=False))
print("  FDR<0.05: %s" % (o[o.fdr < .05]["Network"].tolist() or "无"))
