# -*- coding: utf-8 -*-
"""B7 (R1 最后一点): 区域脑龄与共病理/神经退变生物标志的关系。
R1 指出 PD 的认知衰退部分由黑质纹状体以外的共病理驱动(淀粉样蛋白、tau),
要求讨论并纳入 PPMI 中可用的 amyloid 与 tau 生物标志。
可用变量(取每人最早一次非缺失测量):
  CSFSAA          α-突触核蛋白种子扩增试验 (0=阴性, 1=阳性)   PD n≈522
  IU_ABeta42_CSF / IU_ABeta40_CSF / IU_pTau181_CSF          PD n≈279
  衍生: Abeta42/40 比值, pTau181/Abeta42 比值(AD 病理负荷代理)
传统 Elecsys abeta/tau/ptau 在 PD 组仅 25 例, 不纳入。
"""
import numpy as np, pandas as pd, warnings
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings("ignore")

BASE = r"F:/PD_brainage_data"
CUR = r"F:/Deep_learning_derived_ MRI regional brain age/PPMI_Curated_Data_Cut_Public_20260511.xlsx"
REG = ["sns", "fpn", "drs", "vnt", "dft", "slt", "lng", "adt", "vsl", "lmb"]
EN = {"sns": "Sensorimotor", "fpn": "Frontoparietal", "drs": "Dorsal attention",
      "vnt": "Ventral attention", "dft": "Default mode", "slt": "Salience",
      "lng": "Language", "adt": "Auditory", "vsl": "Visual", "lmb": "Limbic",
      "global": "Global"}
BIO = ["CSFSAA", "IU_ABeta42_CSF", "IU_ABeta40_CSF", "IU_pTau181_CSF"]

df = pd.read_csv(BASE + "/analysis_revision_base.csv")
# 基础表中已带基线一次的生物标志(缺失多), 此处改取"最早一次非缺失", 先移除避免列名冲突
df = df.drop(columns=[c for c in BIO if c in df.columns])
cur = pd.read_excel(CUR, sheet_name="20260511")
cur = cur[cur["PATNO"].isin(df["PATNO"])]

# 每人取最早一次非缺失测量
ORDER = ["BL", "V01", "V02", "V03", "V04", "V05", "V06", "V07", "V08", "V09",
         "V10", "V11", "V12", "V13", "V14", "V15", "SC"]
cur["_o"] = cur["EVENT_ID"].map({e: i for i, e in enumerate(ORDER)}).fillna(99)
cur = cur.sort_values("_o")
bio = cur.groupby("PATNO")[BIO].first().reset_index()
for c in BIO:
    bio[c] = cur[cur[c].notna()].groupby("PATNO")[c].first().reindex(bio["PATNO"]).values

bio["Abeta_ratio"] = bio["IU_ABeta42_CSF"] / bio["IU_ABeta40_CSF"]
bio["pTau_Abeta42"] = bio["IU_pTau181_CSF"] / bio["IU_ABeta42_CSF"]
d = df.merge(bio, on="PATNO", how="left")
pdf = d[d.group == "PD"].copy()

print("=" * 84)
print("生物标志在 PD 影像队列中的可用例数")
print("=" * 84)
for c in BIO + ["Abeta_ratio", "pTau_Abeta42"]:
    print("  %-18s n = %d" % (c, pdf[c].notna().sum()))
pdf["SAA_pos"] = (pdf["CSFSAA"] == 1).astype(float)
pdf.loc[~pdf["CSFSAA"].isin([0, 1]), "SAA_pos"] = np.nan
print("\n  α-syn SAA: 阳性 %d / 阴性 %d"
      % ((pdf.SAA_pos == 1).sum(), (pdf.SAA_pos == 0).sum()))


def star(p):
    return "***" if p < .001 else "**" if p < .01 else "*" if p < .05 else ""


def assoc(data, marker, label, logit=False):
    rows = []
    for r in REG + ["global"]:
        dd = data.dropna(subset=[r + "_z", marker, "age_at_visit", "SEX", "EDUCYRS", "batch"]).rename(
            columns={r + "_z": "x"})
        if len(dd) < 60:
            rows.append([EN[r], np.nan, np.nan, len(dd)])
            continue
        if logit:
            m = smf.logit("%s ~ x + age_at_visit + C(SEX) + EDUCYRS" % marker, dd).fit(disp=0)
            rows.append([EN[r], np.exp(m.params["x"]), m.pvalues["x"], len(dd)])
        else:
            m = smf.mixedlm("%s ~ x + age_at_visit + C(SEX) + EDUCYRS" % marker,
                            dd, groups=dd["batch"]).fit()
            rows.append([EN[r], m.params["x"], m.pvalues["x"], len(dd)])
    o = pd.DataFrame(rows, columns=["Network", "OR_per_SD" if logit else "beta", "p", "n"])
    ok = o["p"].notna()
    o.loc[ok, "fdr"] = multipletests(o.loc[ok, "p"], method="fdr_bh")[1]
    o["sig"] = o["fdr"].map(lambda v: star(v) if pd.notna(v) else "")
    o.to_csv(BASE + "/rev_B7_%s.csv" % label, index=False, encoding="utf-8-sig")
    return o


print("\n" + "=" * 84)
print("B7-1  区域脑龄 与 共病理生物标志 的关联")
print("=" * 84)
for marker, name, lg in [("SAA_pos", "α-syn SAA 阳性 (logistic, OR per 1 SD)", True),
                         ("Abeta_ratio", "CSF Aβ42/Aβ40 比值", False),
                         ("IU_pTau181_CSF", "CSF p-tau181", False),
                         ("pTau_Abeta42", "CSF p-tau181/Aβ42 (AD 病理负荷)", False)]:
    o = assoc(pdf, marker, marker, logit=lg)
    hit = o[o.fdr < .05]
    print("\n  %s   n=%d" % (name, int(o["n"].max())))
    print(o.round(4).to_string(index=False))
    print("    FDR<0.05: %s" % (hit["Network"].tolist() if len(hit) else "无"))

# ---- B7-2: 认知关联在校正共病理后是否存活 (R1 的核心关切) ----
print("\n" + "=" * 84)
print("B7-2  MoCA 关联在校正 AD 共病理标志后是否存活 (核心关切)")
print("=" * 84)
for adj, name in [("pTau_Abeta42", "p-tau181/Aβ42"), ("Abeta_ratio", "Aβ42/Aβ40")]:
    rows = []
    for r in REG + ["global"]:
        dd = pdf.dropna(subset=[r + "_z", "moca", adj, "age_at_visit", "SEX",
                                "EDUCYRS", "batch"]).rename(columns={r + "_z": "x"})
        m0 = smf.mixedlm("moca ~ x + age_at_visit + C(SEX) + EDUCYRS", dd, groups=dd["batch"]).fit()
        m1 = smf.mixedlm("moca ~ x + %s + age_at_visit + C(SEX) + EDUCYRS" % adj,
                         dd, groups=dd["batch"]).fit()
        rows.append([EN[r], m0.params["x"], m0.pvalues["x"], m1.params["x"],
                     m1.pvalues["x"], m1.params[adj], m1.pvalues[adj], len(dd)])
    o = pd.DataFrame(rows, columns=["Network", "b_unadj", "P_unadj", "b_adj",
                                    "P_adj", "b_marker", "P_marker", "n"])
    o["fdr_unadj"] = multipletests(o["P_unadj"], method="fdr_bh")[1]
    o["fdr_adj"] = multipletests(o["P_adj"], method="fdr_bh")[1]
    o.to_csv(BASE + "/rev_B7_moca_adj_%s.csv" % adj, index=False, encoding="utf-8-sig")
    print("\n  校正 %s 后 (同一子样本 n=%d):" % (name, int(o["n"].max())))
    print(o[["Network", "b_unadj", "fdr_unadj", "b_adj", "fdr_adj",
             "b_marker", "P_marker"]].round(4).to_string(index=False))
    print("    校正前 FDR<0.05: %s" % o[o.fdr_unadj < .05]["Network"].tolist())
    print("    校正后 FDR<0.05: %s" % o[o.fdr_adj < .05]["Network"].tolist())
