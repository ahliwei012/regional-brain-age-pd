# -*- coding: utf-8 -*-
"""B2 (R2 两条意见):
1) 把灰白质信号强度比(GWR)——模型的另一个输入特征——作为协变量重跑主关联;
2) 报告各网络 brain-PAD 与该网络平均皮层厚度、平均 GWR 的相关;
3) 报告共线性统计量(VIF), 排除上述模型的共线性问题。
特征文件 features_v2/<PATNO>_features_20k.txt: 20484 顶点 x 2 列 = [厚度, GWR]
"""
import os, numpy as np, pandas as pd, warnings
import statsmodels.formula.api as smf
from statsmodels.stats.outliers_influence import variance_inflation_factor
from statsmodels.stats.multitest import multipletests
from scipy import stats
warnings.filterwarnings("ignore")

BASE = r"F:/PD_brainage_data"
REPO = r"F:/Deep_learning_derived_ MRI regional brain age/regional_Brain_age-main"
EN = {"sns": "Sensorimotor", "fpn": "Frontoparietal", "drs": "Dorsal attention",
      "vnt": "Ventral attention", "dft": "Default mode", "slt": "Salience",
      "lng": "Language", "adt": "Auditory", "vsl": "Visual", "lmb": "Limbic",
      "global": "Global"}
# 网络 -> AAL 区域(取自模型源码), 与主分析完全一致
NET = {"sns": [1, 2, 19, 20, 57, 58, 69, 70], "fpn": [5, 6, 7, 8, 9, 10, 65, 66],
       "drs": [3, 4, 59, 60], "vnt": [11, 12, 13, 14, 63, 64],
       "dft": [21, 22, 25, 26, 27, 28, 35, 36, 65, 66, 67, 68, 85, 86],
       "slt": [29, 30, 31, 32], "lng": [11, 12, 13, 17, 63],
       "adt": [79, 80, 81, 82],
       "vsl": [43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 89, 90],
       "lmb": [15, 16, 23, 24, 33, 34, 39, 40, 83, 84, 87, 88]}
REG = list(NET.keys())

aal = pd.read_table(REPO + "/surface_information/aal_atlas_20k.txt",
                    header=None).values.ravel().astype(int)
assert aal.shape[0] == 20484
masks = {r: np.isin(aal, NET[r]) for r in REG}

df = pd.read_csv(BASE + "/analysis_revision_base.csv")

# ---- 逐被试计算各网络平均厚度与平均 GWR ----
cache = BASE + "/network_thickness_gwr.csv"
if os.path.exists(cache):
    feat = pd.read_csv(cache)
    print("[复用] %s" % cache)
else:
    rows = []
    for p in df["PATNO"].astype(int):
        f = BASE + "/features_v2/%d_features_20k.txt" % p
        if not os.path.exists(f):
            rows.append({"PATNO": p})
            continue
        a = pd.read_csv(f, header=None, sep=r"\s+").values
        rec = {"PATNO": p}
        for r in REG:
            rec[r + "_thick"] = float(a[masks[r], 0].mean())
            rec[r + "_gwr"] = float(a[masks[r], 1].mean())
        rows.append(rec)
    feat = pd.DataFrame(rows)
    feat.to_csv(cache, index=False, encoding="utf-8-sig")
    print("[已保存] %s" % cache)

df = df.merge(feat, on="PATNO", how="left")
pdf = df[df.group == "PD"].copy()
print("各网络平均 GWR 可用例数: %d / %d\n" % (pdf["sns_gwr"].notna().sum(), len(pdf)))


def star(p):
    return "***" if p < .001 else "**" if p < .01 else "*" if p < .05 else ""


# ============ 1) brain-PAD 与 厚度 / GWR 的相关 ============
print("=" * 78)
print("B2-1  各网络 brain-PAD 与该网络平均皮层厚度、平均 GWR 的相关 (PD 组)")
print("=" * 78)
rows = []
for r in REG:
    d = pdf.dropna(subset=[r + "_PAD", r + "_thick", r + "_gwr"])
    rt, pt = stats.pearsonr(d[r + "_PAD"], d[r + "_thick"])
    rg, pg = stats.pearsonr(d[r + "_PAD"], d[r + "_gwr"])
    rtg, ptg = stats.pearsonr(d[r + "_thick"], d[r + "_gwr"])
    rows.append([EN[r], rt, pt, rg, pg, rtg, ptg, len(d)])
c1 = pd.DataFrame(rows, columns=["Network", "r_PAD_thickness", "P_thick",
                                 "r_PAD_GWR", "P_gwr", "r_thickness_GWR", "P_tg", "n"])
c1.to_csv(BASE + "/rev_B2_correlations.csv", index=False, encoding="utf-8-sig")
print(c1.round(4).to_string(index=False))

# ============ 2) VIF 共线性 ============
print("\n" + "=" * 78)
print("B2-2  共线性检验 (VIF): 模型 outcome ~ brain-PAD + thickness + GWR + 协变量")
print("=" * 78)
rows = []
for r in REG:
    d = pdf.dropna(subset=[r + "_z", r + "_thick", r + "_gwr", "age_at_visit",
                           "EDUCYRS"]).copy()
    X = d[[r + "_z", r + "_thick", r + "_gwr", "age_at_visit", "EDUCYRS"]].copy()
    X["SEX"] = (d["SEX"] == d["SEX"].iloc[0]).astype(float)
    X.insert(0, "const", 1.0)
    v = [variance_inflation_factor(X.values, i) for i in range(1, X.shape[1])]
    rows.append([EN[r]] + [round(x, 3) for x in v])
vif = pd.DataFrame(rows, columns=["Network", "VIF_brainPAD", "VIF_thickness",
                                  "VIF_GWR", "VIF_age", "VIF_educ", "VIF_sex"])
vif.to_csv(BASE + "/rev_B2_VIF.csv", index=False, encoding="utf-8-sig")
print(vif.to_string(index=False))
print("\n  判读: VIF < 5 视为无实质共线性, < 10 为可接受。最大 VIF = %.2f"
      % vif[[c for c in vif.columns if c.startswith("VIF")]].values.max())


# ============ 3) 控制 GWR / 同时控制厚度与 GWR 后重跑主关联 ============
def fam(data, outcome, terms, label):
    rows = []
    for r in REG + ["global"]:
        tt = [t.replace("<R>", r) for t in terms]
        need = [t for t in tt if t in data.columns] + [outcome, "age_at_visit",
                                                       "SEX", "EDUCYRS", "batch"]
        d = data.dropna(subset=need).copy()
        if r == "global" and any("_gwr" in t or "_thick" in t for t in tt):
            rows.append([EN[r], np.nan, np.nan, np.nan, len(d)])
            continue
        f = "%s ~ %s + age_at_visit + C(SEX) + EDUCYRS" % (outcome, " + ".join(tt))
        m = smf.mixedlm(f, d, groups=d["batch"]).fit()
        k = r + "_z"
        rows.append([EN[r], m.params[k], m.pvalues[k], m.params.get(r + "_gwr", np.nan), len(d)])
    out = pd.DataFrame(rows, columns=["Network", "beta_brainPAD", "P", "beta_GWR", "n"])
    ok = out["P"].notna()
    out.loc[ok, "fdr"] = multipletests(out.loc[ok, "P"], method="fdr_bh")[1]
    out["sig"] = out["fdr"].map(lambda p: star(p) if pd.notna(p) else "")
    out.to_csv(BASE + "/rev_B2_%s.csv" % label, index=False, encoding="utf-8-sig")
    return out


for outcome, tag in [("updrs3_score", "UPDRS-III"), ("moca", "MoCA")]:
    print("\n" + "=" * 78)
    print("B2-3  %s : brain-PAD 在控制 GWR / 控制厚度+GWR 后是否存活" % tag)
    print("=" * 78)
    base = fam(pdf, outcome, ["<R>_z"], "%s_base" % outcome)
    gwr = fam(pdf, outcome, ["<R>_z", "<R>_gwr"], "%s_adjGWR" % outcome)
    both = fam(pdf, outcome, ["<R>_z", "<R>_thick", "<R>_gwr"], "%s_adjBoth" % outcome)
    cmp = base[["Network", "beta_brainPAD", "fdr"]].rename(
        columns={"beta_brainPAD": "b_base", "fdr": "fdr_base"})
    cmp["b_adjGWR"] = gwr["beta_brainPAD"]
    cmp["fdr_adjGWR"] = gwr["fdr"]
    cmp["beta_GWR_itself"] = gwr["beta_GWR"]
    cmp["b_adjBoth"] = both["beta_brainPAD"]
    cmp["fdr_adjBoth"] = both["fdr"]
    cmp.to_csv(BASE + "/rev_B2_%s_summary.csv" % outcome, index=False, encoding="utf-8-sig")
    print(cmp.round(4).to_string(index=False))
    print("  FDR<0.05  仅PAD: %s" % cmp[cmp.fdr_base < .05]["Network"].tolist())
    print("            +GWR: %s" % cmp[cmp.fdr_adjGWR < .05]["Network"].tolist())
    print("     +厚度+GWR  : %s" % cmp[cmp.fdr_adjBoth < .05]["Network"].tolist())
