# -*- coding: utf-8 -*-
"""ComBat谐和化区域脑龄(批次=扫描仪) -> 重算PAD -> 重跑分析, 对比谐和化前后。"""
import pandas as pd, numpy as np
from neuroCombat import neuroCombat
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import warnings; warnings.filterwarnings("ignore")

BASE=r"F:/PD_brainage_data"
REG=["sns","fpn","drs","vnt","dft","slt","lng","adt","vsl","lmb"]; ALL=REG+["global","AD"]
EN={"sns":"Sensorimotor","fpn":"Frontoparietal","drs":"Dorsal att","vnt":"Ventral att","dft":"Default",
    "slt":"Salience","lng":"Language","adt":"Auditory","vsl":"Visual","lmb":"Limbic","global":"Global"}

ba=pd.read_csv(f"{BASE}/brain_ages/all_brain_ages.csv"); ba["PATNO"]=ba["subject"].str.extract(r"(\d+)").astype(int)
cl=pd.read_csv(f"{BASE}/master_clinical.csv")
sc=pd.read_csv(f"{BASE}/scanner_info.csv")[["PATNO","Manuf","scanner"]]
df=ba.merge(cl,on="PATNO").merge(sc,on="PATNO").dropna(subset=["Manuf","age_at_visit","SEX"]).reset_index(drop=True)

# 批次: 厂商×型号, <10例并入"厂商_other"
df["batch"]=df["scanner"]; vc=df["batch"].value_counts()
df.loc[df["batch"].isin(vc[vc<10].index),"batch"]=df["Manuf"].astype(str)+"_other"
print(f"样本={len(df)}  批次数={df['batch'].nunique()}")
print("批次分布:",df["batch"].value_counts().to_dict())

# ---- ComBat 谐和12个脑龄 (保留 age/sex/group) ----
covars=df[["batch","age_at_visit","SEX"]].copy(); covars["grp01"]=(df["group"]=="PD").astype(int)
out=neuroCombat(dat=df[ALL].T.values, covars=covars, batch_col="batch",
                categorical_cols=["SEX","grp01"], continuous_cols=["age_at_visit"])
harm=out["data"].T
for j,c in enumerate(ALL): df[c+"_h"]=harm[:,j]
print("ComBat谐和完成")

# ---- PAD(谐和后) + delta校正(HC) ----
hc=df[df.group=="HC"]
for c in REG+["global"]:
    a,b=np.polyfit(hc["age_at_visit"],hc[c+"_h"],1); df[c+"_PADh"]=df[c+"_h"]-(a*df["age_at_visit"]+b)
    a2,b2=np.polyfit(hc["age_at_visit"],hc[c],1);     df[c+"_PAD"]=df[c]-(a2*df["age_at_visit"]+b2)

pdf=df[df.group=="PD"]
def assoc(outcome,suffix):
    rows=[]
    for c in REG+["global"]:
        d=pdf.dropna(subset=[c+suffix,outcome,"age_at_visit","SEX","EDUCYRS"])
        m=smf.ols(f"{outcome} ~ Q('{c+suffix}') + age_at_visit + C(SEX) + EDUCYRS",d).fit()
        rows.append([c,m.params[f"Q('{c+suffix}')"],m.pvalues[f"Q('{c+suffix}')"]])
    r=pd.DataFrame(rows,columns=["r","beta","p"]); r["fdr"]=multipletests(r["p"],method="fdr_bh")[1]; return r

print("\n========= UPDRS-III 关联: 谐和前 vs 谐和后 (FDR<0.05标*) =========")
u0=assoc("updrs3_score","_PAD"); u1=assoc("updrs3_score","_PADh")
print(f"{'区域':14s}{'前β':>9s}{'前FDR':>9s}{'后β':>9s}{'后FDR':>9s}")
for i,c in enumerate(REG+["global"]):
    s0="*" if u0.fdr[i]<0.05 else ""; s1="*" if u1.fdr[i]<0.05 else ""
    print(f"{EN[c]:14s}{u0.beta[i]:9.3f}{u0.fdr[i]:8.3f}{s0:1s}{u1.beta[i]:9.3f}{u1.fdr[i]:8.3f}{s1:1s}")
print(f"谐和前FDR<0.05: {(u0.fdr<0.05).sum()}个  谐和后: {(u1.fdr<0.05).sum()}个")

# PD vs HC (谐和后)
print("\n========= PD vs HC 区域PAD: 谐和后 (校正age/sex/educ) =========")
rows=[]
for c in REG+["global"]:
    d=df.dropna(subset=[c+"_PADh","age_at_visit","SEX","EDUCYRS"]).copy(); d["g"]=(d.group=="PD").astype(int)
    m=smf.ols(f"Q('{c}_PADh') ~ g + age_at_visit + C(SEX) + EDUCYRS",d).fit()
    rows.append([c,m.params["g"],m.pvalues["g"]])
r=pd.DataFrame(rows,columns=["r","beta","p"]); r["fdr"]=multipletests(r["p"],method="fdr_bh")[1]
for i,c in enumerate(REG+["global"]):
    print(f"  {EN[c]:14s} β={r.beta[i]:7.3f}  FDR={r.fdr[i]:.3f} {'*' if r.fdr[i]<0.05 else ''}")

df.to_csv(f"{BASE}/analysis_merged_combat.csv",index=False,encoding="utf-8-sig")
print("\n已保存 analysis_merged_combat.csv (含谐和后PAD列 *_PADh)")
