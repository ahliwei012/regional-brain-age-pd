# -*- coding: utf-8 -*-
"""SAP核心分析: delta校正 -> PD vs HC区域PAD -> PAD↔UPDRS-III/MoCA -> RF特征重要性。"""
import pandas as pd, numpy as np
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
import warnings; warnings.filterwarnings("ignore")

BASE = r"F:/PD_brainage_data"
REGIONS = ["sns","fpn","drs","vnt","dft","slt","lng","adt","vsl","lmb"]
NAME = {"sns":"感觉运动","fpn":"额顶","drs":"背侧注意","vnt":"腹侧注意","dft":"默认模式",
        "slt":"突显","lng":"语言","adt":"听觉","vsl":"视觉","lmb":"边缘"}

ba = pd.read_csv(f"{BASE}/brain_ages/all_brain_ages.csv")
ba["PATNO"] = ba["subject"].str.extract(r"(\d+)").astype(int)
cl = pd.read_csv(f"{BASE}/master_clinical.csv")
df = ba.merge(cl, on="PATNO")
print(f"合并后: {len(df)} 例  (PD={sum(df.group=='PD')}, HC={sum(df.group=='HC')})")

# ---- delta校正: 在HC上拟合 brain~age, 全体去偏 ----
hc = df[df.group == "HC"]
for col in REGIONS + ["global"]:
    a, b = np.polyfit(hc["age_at_visit"], hc[col], 1)
    df[col+"_PAD"] = df[col] - (a*df["age_at_visit"] + b)
print("delta校正完成 (HC上拟合age趋势, 去偏)")

# ================= 1. PD vs HC 区域PAD =================
print("\n================ 1. PD vs HC 区域脑龄PAD (校正 age/sex/educ) ================")
rows = []
for col in REGIONS + ["global"]:
    d = df.dropna(subset=[col+"_PAD","age_at_visit","SEX","EDUCYRS"]).copy()
    d["grp"] = (d.group == "PD").astype(int)
    m = smf.ols(f"Q('{col}_PAD') ~ grp + age_at_visit + C(SEX) + EDUCYRS", d).fit()
    beta, p = m.params["grp"], m.pvalues["grp"]
    pd_m = d[d.grp==1][col+"_PAD"]; hc_m = d[d.grp==0][col+"_PAD"]
    sd_pool = np.sqrt((pd_m.var()*(len(pd_m)-1)+hc_m.var()*(len(hc_m)-1))/(len(pd_m)+len(hc_m)-2))
    cohend = (pd_m.mean()-hc_m.mean())/sd_pool
    rows.append([col, NAME.get(col,col), pd_m.mean(), hc_m.mean(), beta, cohend, p])
r1 = pd.DataFrame(rows, columns=["region","名称","PD_PAD","HC_PAD","beta","Cohen_d","p"])
r1["FDR_p"] = multipletests(r1["p"], method="fdr_bh")[1]
r1["sig"] = r1["FDR_p"].apply(lambda x: "***" if x<0.001 else "**" if x<0.01 else "*" if x<0.05 else "")
print(r1.round(4).to_string(index=False))

# ================= 2. 区域PAD ↔ UPDRS-III (PD内) =================
print("\n================ 2. 区域PAD ↔ UPDRS-III (PD内, 校正 age/sex/educ) ================")
pdf = df[df.group=="PD"].copy()
rows = []
for col in REGIONS + ["global"]:
    d = pdf.dropna(subset=[col+"_PAD","updrs3_score","age_at_visit","SEX","EDUCYRS"])
    m = smf.ols(f"updrs3_score ~ Q('{col}_PAD') + age_at_visit + C(SEX) + EDUCYRS", d).fit()
    rows.append([col, NAME.get(col,col), m.params[f"Q('{col}_PAD')"], m.pvalues[f"Q('{col}_PAD')"], len(d)])
r2 = pd.DataFrame(rows, columns=["region","名称","beta(对UPDRS-III)","p","n"])
r2["FDR_p"] = multipletests(r2["p"], method="fdr_bh")[1]
r2["sig"] = r2["FDR_p"].apply(lambda x: "***" if x<0.001 else "**" if x<0.01 else "*" if x<0.05 else "")
print(r2.round(4).to_string(index=False))

# ================= 3. 区域PAD ↔ MoCA (PD内) =================
print("\n================ 3. 区域PAD ↔ MoCA (PD内, 校正 age/sex/educ) ================")
rows = []
for col in REGIONS + ["global"]:
    d = pdf.dropna(subset=[col+"_PAD","moca","age_at_visit","SEX","EDUCYRS"])
    m = smf.ols(f"moca ~ Q('{col}_PAD') + age_at_visit + C(SEX) + EDUCYRS", d).fit()
    rows.append([col, NAME.get(col,col), m.params[f"Q('{col}_PAD')"], m.pvalues[f"Q('{col}_PAD')"], len(d)])
r3 = pd.DataFrame(rows, columns=["region","名称","beta(对MoCA)","p","n"])
r3["FDR_p"] = multipletests(r3["p"], method="fdr_bh")[1]
r3["sig"] = r3["FDR_p"].apply(lambda x: "***" if x<0.001 else "**" if x<0.01 else "*" if x<0.05 else "")
print(r3.round(4).to_string(index=False))

# ================= 4. RF特征重要性 (预测运动好/差) =================
print("\n================ 4. 随机森林特征重要性 (预测运动好/差, 1000次bootstrap) ================")
feat_cols = [c+"_PAD" for c in REGIONS]
d = pdf.dropna(subset=feat_cols+["updrs3_score"]).copy()
thr = d["updrs3_score"].median()
d["poor"] = (d["updrs3_score"] > thr).astype(int)
X = d[feat_cols].values; y = d["poor"].values
print(f"PD n={len(d)}, 运动差(UPDRS-III>{thr:.0f})={y.sum()}, 好={len(y)-y.sum()}")
imp = np.zeros((1000, len(feat_cols)))
for i in range(1000):
    Xtr,_,ytr,_ = train_test_split(X, y, test_size=0.2, random_state=i, stratify=y)
    rf = RandomForestClassifier(n_estimators=200, random_state=i, n_jobs=-1).fit(Xtr, ytr)
    imp[i] = rf.feature_importances_
mi = imp.mean(0); si = imp.std(0)
r4 = pd.DataFrame({"region":[NAME[c] for c in REGIONS], "importance":mi, "sd":si}).sort_values("importance", ascending=False)
print(r4.round(4).to_string(index=False))

# 保存
r1.to_csv(f"{BASE}/results_pd_vs_hc.csv", index=False, encoding="utf-8-sig")
r2.to_csv(f"{BASE}/results_pad_updrs.csv", index=False, encoding="utf-8-sig")
r3.to_csv(f"{BASE}/results_pad_moca.csv", index=False, encoding="utf-8-sig")
r4.to_csv(f"{BASE}/results_rf_importance.csv", index=False, encoding="utf-8-sig")
df.to_csv(f"{BASE}/analysis_merged.csv", index=False, encoding="utf-8-sig")
print("\n结果表已保存到 F:/PD_brainage_data/results_*.csv")
