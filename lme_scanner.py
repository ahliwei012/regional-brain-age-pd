# -*- coding: utf-8 -*-
"""稳健性检验: (A)LME加扫描仪随机效应(贴合原文); (B)PD内ComBat保留UPDRS后再测。"""
import pandas as pd, numpy as np
from neuroCombat import neuroCombat
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import warnings; warnings.filterwarnings("ignore")

BASE=r"F:/PD_brainage_data"
REG=["sns","fpn","drs","vnt","dft","slt","lng","adt","vsl","lmb"]
EN={"sns":"Sensorimotor","fpn":"Frontoparietal","drs":"Dorsal att","vnt":"Ventral att","dft":"Default",
    "slt":"Salience","lng":"Language","adt":"Auditory","vsl":"Visual","lmb":"Limbic"}

df=pd.read_csv(f"{BASE}/analysis_merged_combat.csv")  # 已含 _PAD(原) 和 scanner/batch
pdf=df[df.group=="PD"].copy()

# ===== (A) LME: updrs3 ~ PAD + age + sex + educ + (1|scanner)  (原文做法: site随机效应) =====
print("===== (A) LME 加扫描仪随机效应 (未谐和PAD; 贴合原文site随机效应) =====")
rows=[]
for c in REG:
    d=pdf.dropna(subset=[c+"_PAD","updrs3_score","age_at_visit","SEX","EDUCYRS","scanner"]).copy()
    d=d.rename(columns={c+"_PAD":"pad"})
    try:
        m=smf.mixedlm("updrs3_score ~ pad + age_at_visit + C(SEX) + EDUCYRS", d, groups=d["scanner"]).fit()
        rows.append([c, m.params["pad"], m.pvalues["pad"]])
    except Exception as e:
        rows.append([c, np.nan, np.nan])
r=pd.DataFrame(rows,columns=["r","beta","p"]); r["fdr"]=multipletests(r["p"].fillna(1),method="fdr_bh")[1]
for i,c in enumerate(REG):
    print(f"  {EN[c]:14s} β={r.beta[i]:7.3f}  p={r.p[i]:.3f}  FDR={r.fdr[i]:.3f} {'*' if r.fdr[i]<0.05 else ''}")
print(f"  >>> LME+site随机效应 FDR<0.05: {(r.fdr<0.05).sum()}个")

# ===== (B) PD内ComBat保留UPDRS(+age/sex), 再测 =====
print("\n===== (B) PD内ComBat(保留 age/sex/UPDRS) 后再测 =====")
ALL=REG+["global","AD"]
ba=pd.read_csv(f"{BASE}/brain_ages/all_brain_ages.csv"); ba["PATNO"]=ba["subject"].str.extract(r"(\d+)").astype(int)
p2=ba.merge(df[["PATNO","group","age_at_visit","SEX","EDUCYRS","updrs3_score","Manuf","scanner","batch"]],on="PATNO")
p2=p2[(p2.group=="PD")].dropna(subset=["batch","age_at_visit","SEX","updrs3_score"]).reset_index(drop=True)
# 批次≥2
vc=p2["batch"].value_counts(); p2=p2[p2["batch"].isin(vc[vc>=3].index)].reset_index(drop=True)
cov=p2[["batch","age_at_visit","SEX","updrs3_score"]].copy()
out=neuroCombat(dat=p2[ALL].T.values, covars=cov, batch_col="batch",
                categorical_cols=["SEX"], continuous_cols=["age_at_visit","updrs3_score"])
H=out["data"].T
for j,c in enumerate(ALL): p2[c+"_h"]=H[:,j]
# delta校正用全体PD(无HC, 用PD自身age趋势去偏仅为去年龄, 主要看关联)
rows=[]
for c in REG:
    a,b=np.polyfit(p2["age_at_visit"],p2[c+"_h"],1); p2["padh"]=p2[c+"_h"]-(a*p2["age_at_visit"]+b)
    d=p2.dropna(subset=["padh","updrs3_score","age_at_visit","SEX","EDUCYRS"])
    m=smf.ols("updrs3_score ~ padh + age_at_visit + C(SEX) + EDUCYRS",d).fit()
    rows.append([c,m.params["padh"],m.pvalues["padh"]])
r2=pd.DataFrame(rows,columns=["r","beta","p"]); r2["fdr"]=multipletests(r2["p"],method="fdr_bh")[1]
for i,c in enumerate(REG):
    print(f"  {EN[c]:14s} β={r2.beta[i]:7.3f}  p={r2.p[i]:.3f}  FDR={r2.fdr[i]:.3f} {'*' if r2.fdr[i]<0.05 else ''}")
print(f"  >>> PD内ComBat(保留UPDRS)后 FDR<0.05: {(r2.fdr<0.05).sum()}个")
