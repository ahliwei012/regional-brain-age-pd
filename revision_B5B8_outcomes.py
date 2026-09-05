# -*- coding: utf-8 -*-
"""B8 (R2): MDS-UPDRS-III 四个运动亚域(震颤/强直/运动迟缓/中轴-PIGD)与区域脑龄的关联
   B5 (R1-4, R2-2): 超出 UPDRS-III 与 MoCA 的临床结局
       —— 日常生活活动(UPDRS-I/II, 改良 Schwab-England)、疾病分期(H&Y)、
          MCI 转化(cogstate)、以及域特异认知(记忆/视空间/工作记忆/加工速度/语义流畅)
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

df = pd.read_csv(BASE + "/analysis_revision_base.csv")

# =========================== B8: 运动亚域 ===========================
SUB = {
    "Tremor": ["NP3PTRMR", "NP3PTRML", "NP3KTRMR", "NP3KTRML", "NP3RTARU",
               "NP3RTALU", "NP3RTARL", "NP3RTALL", "NP3RTALJ", "NP3RTCON"],
    "Rigidity": ["NP3RIGN", "NP3RIGRU", "NP3RIGLU", "NP3RIGRL", "NP3RIGLL"],
    "Bradykinesia": ["NP3FTAPR", "NP3FTAPL", "NP3HMOVR", "NP3HMOVL", "NP3PRSPR",
                     "NP3PRSPL", "NP3TTAPR", "NP3TTAPL", "NP3LGAGR", "NP3LGAGL", "NP3BRADY"],
    "Axial_PIGD": ["NP3SPCH", "NP3FACXP", "NP3RISNG", "NP3GAIT", "NP3FRZGT",
                   "NP3PSTBL", "NP3POSTR"],
}
p3 = pd.read_csv(BASE + "/MDS-UPDRS_Part/MDS-UPDRS_Part_III_04Sep2026.csv", low_memory=False)
print("Part III 原始: %d 行, PDSTATE 分布 %s" % (len(p3), p3["PDSTATE"].value_counts(dropna=False).to_dict()))

# 基线、OFF 态(或未治疗)的检查, 与主分析使用的 updrs3_score(OFF) 保持一致
bl = p3[(p3["EVENT_ID"] == "BL") & (p3["PDSTATE"].isin(["OFF"]) | p3["PDSTATE"].isna())].copy()
bl = bl.sort_values("PDSTATE").drop_duplicates("PATNO", keep="first")
for k, items in SUB.items():
    have = [c for c in items if c in bl.columns]
    bl[k] = bl[have].sum(axis=1, min_count=len(have))
sub = bl[["PATNO"] + list(SUB)].copy()
print("基线 OFF/未治疗 记录: %d 例" % len(sub))

d8 = df.merge(sub, on="PATNO", how="left")
pdf8 = d8[d8.group == "PD"].copy()
print("与影像数据匹配的 PD: %d 例\n" % pdf8["Tremor"].notna().sum())
print("各亚域得分 (PD 基线, 均值±SD):")
for k in SUB:
    print("  %-14s %.2f +/- %.2f  (n=%d)" % (k, pdf8[k].mean(), pdf8[k].std(), pdf8[k].notna().sum()))


def star(p):
    return "***" if p < .001 else "**" if p < .01 else "*" if p < .05 else ""


def fam(data, outcome, label, cov="age_at_visit + C(SEX) + EDUCYRS"):
    rows = []
    for r in REG + ["global"]:
        need = [r + "_z", outcome, "age_at_visit", "SEX", "EDUCYRS", "batch"]
        d = data.dropna(subset=need).rename(columns={r + "_z": "x"})
        if len(d) < 40 or d[outcome].nunique() < 3:
            rows.append([EN[r], np.nan, np.nan, len(d)])
            continue
        m = smf.mixedlm("%s ~ x + %s" % (outcome, cov), d, groups=d["batch"]).fit()
        rows.append([EN[r], m.params["x"], m.pvalues["x"], len(d)])
    out = pd.DataFrame(rows, columns=["Network", "beta", "p", "n"])
    ok = out["p"].notna()
    out.loc[ok, "fdr"] = multipletests(out.loc[ok, "p"], method="fdr_bh")[1]
    out["sig"] = out["fdr"].map(lambda v: star(v) if pd.notna(v) else "")
    out.to_csv(BASE + "/rev_%s.csv" % label, index=False, encoding="utf-8-sig")
    return out


print("\n" + "=" * 84)
print("B8  区域脑龄与 MDS-UPDRS-III 四个运动亚域 (标准化 β, 每 1 SD brain-PAD)")
print("=" * 84)
tab = pd.DataFrame({"Network": [EN[r] for r in REG + ["global"]]})
for k in SUB:
    o = fam(pdf8, k, "B8_%s" % k.lower())
    tab[k] = o["beta"].round(3).astype(str) + o["sig"]
    tab[k + "_fdr"] = o["fdr"].round(4)
print(tab[["Network"] + list(SUB)].to_string(index=False))
print()
for k in SUB:
    o = pd.read_csv(BASE + "/rev_B8_%s.csv" % k.lower())
    print("  %-14s FDR<0.05: %s" % (k, o[o.fdr < .05]["Network"].tolist()))

# =========================== B5: 扩展临床结局 ===========================
pdf = df[df.group == "PD"].copy()
OUT = [
    ("updrs1_score", "MDS-UPDRS I (非运动日常体验)"),
    ("updrs2_score", "MDS-UPDRS II (运动性日常生活体验)"),
    ("MSEADLG", "改良 Schwab-England ADL (%)"),
    ("NHY", "Hoehn & Yahr 分期"),
    ("hvlt_immediaterecall", "HVLT 即刻回忆 (言语记忆)"),
    ("hvlt_retention", "HVLT 保持率 (延迟记忆)"),
    ("bjlot", "Benton 线方向判断 (视空间)"),
    ("DVS_LNS", "字母-数字排序 (工作记忆)"),
    ("DVT_SDM", "符号数字模式 (加工速度)"),
    ("DVS_SFTANIM", "语义流畅 (动物)"),
    ("COG_COMPOSITE_INT", "认知综合分"),
]
print("\n" + "=" * 84)
print("B5  区域脑龄与扩展临床结局 (标准化 β, 每 1 SD brain-PAD)")
print("=" * 84)
summary = []
for col, name in OUT:
    if col not in pdf.columns or pdf[col].notna().sum() < 100:
        print("  [跳过] %s (n=%s)" % (name, pdf[col].notna().sum() if col in pdf else 0))
        continue
    o = fam(pdf, col, "B5_%s" % col)
    hit = o[o.fdr < .05]
    print("\n  %s   n=%d" % (name, int(o["n"].max())))
    print("    FDR<0.05: %s" % (hit["Network"].tolist() if len(hit) else "无"))
    if len(hit):
        for _, x in hit.iterrows():
            print("       %-18s beta=%+.3f  FDR=%.4f" % (x.Network, x.beta, x.fdr))
    summary.append([name, col, int(o["n"].max()), "; ".join(hit["Network"])])

pd.DataFrame(summary, columns=["Outcome", "variable", "n", "networks_FDR<0.05"]).to_csv(
    BASE + "/rev_B5_summary.csv", index=False, encoding="utf-8-sig")

# MCI: cogstate 1=正常, >1 为 MCI/痴呆
if "cogstate" in pdf.columns:
    pdf["MCI"] = (pdf["cogstate"] > 1).astype(float)
    pdf.loc[pdf["cogstate"].isna(), "MCI"] = np.nan
    print("\n  基线认知状态: %s (1=正常认知)" % pdf["cogstate"].value_counts().to_dict())
    rows = []
    for r in REG + ["global"]:
        d = pdf.dropna(subset=[r + "_z", "MCI", "age_at_visit", "SEX", "EDUCYRS"]).rename(
            columns={r + "_z": "x"})
        m = smf.logit("MCI ~ x + age_at_visit + C(SEX) + EDUCYRS", d).fit(disp=0)
        rows.append([EN[r], np.exp(m.params["x"]), m.pvalues["x"], len(d)])
    o = pd.DataFrame(rows, columns=["Network", "OR_per_SD", "p", "n"])
    o["fdr"] = multipletests(o["p"], method="fdr_bh")[1]
    o["sig"] = o["fdr"].map(star)
    o.to_csv(BASE + "/rev_B5_MCI_baseline.csv", index=False, encoding="utf-8-sig")
    print("\n  基线 MCI/痴呆 (logistic, OR per 1 SD brain-PAD):")
    print(o.round(4).to_string(index=False))
