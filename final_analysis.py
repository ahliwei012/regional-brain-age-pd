# -*- coding: utf-8 -*-
"""最终分析(550 PD + 261 HC): LME+站点随机效应(主) + ComBat(敏感性)。"""
import pandas as pd, numpy as np
from neuroCombat import neuroCombat
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import warnings; warnings.filterwarnings("ignore")

BASE=r"F:/PD_brainage_data"
REG=["sns","fpn","drs","vnt","dft","slt","lng","adt","vsl","lmb"]; ALL=REG+["global","AD"]
EN={"sns":"Sensorimotor","fpn":"Frontoparietal","drs":"Dorsal attention","vnt":"Ventral attention","dft":"Default mode",
    "slt":"Salience","lng":"Language","adt":"Auditory","vsl":"Visual","lmb":"Limbic","global":"Global"}

ba=pd.read_csv(f"{BASE}/brain_ages/all_brain_ages_final.csv"); ba["PATNO"]=ba["subject"].str.extract(r"(\d+)").astype(int)
cl=pd.read_csv(f"{BASE}/master_clinical.csv")
sc=pd.read_csv(f"{BASE}/scanner_info.csv")[["PATNO","Manuf","scanner"]]
df=ba.merge(cl,on="PATNO").merge(sc,on="PATNO").dropna(subset=["Manuf","age_at_visit","SEX"]).reset_index(drop=True)
# 批次: 厂商×型号, <10并入 厂商_other
df["batch"]=df["scanner"]; vc=df["batch"].value_counts()
df.loc[df["batch"].isin(vc[vc<10].index),"batch"]=df["Manuf"].astype(str)+"_other"
print(f"最终样本: {len(df)} (PD={sum(df.group=='PD')}, HC={sum(df.group=='HC')}), 批次={df['batch'].nunique()}")

# delta校正(HC上拟合, 未谐和PAD)
hc=df[df.group=="HC"]
for c in REG+["global"]:
    a,b=np.polyfit(hc["age_at_visit"],hc[c],1); df[c+"_PAD"]=df[c]-(a*df["age_at_visit"]+b)
    df[c+"_z"]=(df[c+"_PAD"]-df[c+"_PAD"].mean())/df[c+"_PAD"].std()
pdf=df[df.group=="PD"].copy()

def lme(data,outcome,pred):
    d=data.dropna(subset=[pred,outcome,"age_at_visit","SEX","EDUCYRS","batch"]).copy().rename(columns={pred:"x"})
    m=smf.mixedlm(f"{outcome} ~ x + age_at_visit + C(SEX) + EDUCYRS", d, groups=d["batch"]).fit()
    ci=m.conf_int().loc["x"]
    return m.params["x"], ci[0], ci[1], m.pvalues["x"], len(d)

def run_family(data,outcome,tag):
    rows=[]
    for c in REG+["global"]:
        b,lo,hi,p,n=lme(data,outcome,c+"_z")
        rows.append([c,EN[c],b,lo,hi,p,n])
    r=pd.DataFrame(rows,columns=["r","name","beta","lo","hi","p","n"])
    r["fdr"]=multipletests(r["p"],method="fdr_bh")[1]
    r["sig"]=r["fdr"].apply(lambda x:"***" if x<0.001 else "**" if x<0.01 else "*" if x<0.05 else "")
    r.to_csv(f"{BASE}/results_final_{tag}.csv",index=False,encoding="utf-8-sig")
    print(f"\n=== {tag} (LME+站点随机效应, 标准化β) ===")
    print(r[["name","beta","fdr","sig","n"]].round(4).to_string(index=False))
    print(f"  FDR<0.05: {(r.fdr<0.05).sum()}个 -> {r[r.fdr<0.05]['name'].tolist()}")
    return r

# PD vs HC: 组别做预测(仍用LME+site)
dfz=df.copy(); dfz["grp"]=(dfz.group=="PD").astype(int)
rows=[]
for c in REG+["global"]:
    d=dfz.dropna(subset=[c+"_PAD","age_at_visit","SEX","EDUCYRS","batch"]).copy().rename(columns={c+"_PAD":"y"})
    m=smf.mixedlm("y ~ grp + age_at_visit + C(SEX) + EDUCYRS", d, groups=d["batch"]).fit()
    ci=m.conf_int().loc["grp"]; rows.append([c,EN[c],m.params["grp"],ci[0],ci[1],m.pvalues["grp"],len(d)])
pvh=pd.DataFrame(rows,columns=["r","name","beta","lo","hi","p","n"]); pvh["fdr"]=multipletests(pvh["p"],method="fdr_bh")[1]
pvh["sig"]=pvh["fdr"].apply(lambda x:"***" if x<0.001 else "**" if x<0.01 else "*" if x<0.05 else "")
pvh.to_csv(f"{BASE}/results_final_pdvshc.csv",index=False,encoding="utf-8-sig")
print("\n=== PD vs HC 区域PAD (LME+站点) ===")
print(pvh[["name","beta","fdr","sig"]].round(4).to_string(index=False))
print(f"  FDR<0.05: {(pvh.fdr<0.05).sum()}个 -> {pvh[pvh.fdr<0.05]['name'].tolist()}")

U=run_family(pdf,"updrs3_score","updrs")
M=run_family(pdf,"moca","moca")

# ComBat敏感性(UPDRS)
covars=df[["batch","age_at_visit","SEX"]].copy(); covars["g"]=(df.group=="PD").astype(int)
H=neuroCombat(dat=df[ALL].T.values,covars=covars,batch_col="batch",categorical_cols=["SEX","g"],continuous_cols=["age_at_visit"])["data"].T
for j,c in enumerate(ALL): df[c+"_h"]=H[:,j]
hc2=df[df.group=="HC"]; pdf2=df[df.group=="PD"]
rows=[]
for c in REG+["global"]:
    a,b=np.polyfit(hc2["age_at_visit"],hc2[c+"_h"],1); pad=pdf2[c+"_h"]-(a*pdf2["age_at_visit"]+b)
    d=pdf2.assign(padh=(pad-pad.mean())/pad.std()).dropna(subset=["padh","updrs3_score","age_at_visit","SEX","EDUCYRS"])
    m=smf.ols("updrs3_score ~ padh + age_at_visit + C(SEX) + EDUCYRS",d).fit()
    rows.append([EN[c],m.params["padh"],m.pvalues["padh"]])
cb=pd.DataFrame(rows,columns=["name","beta","p"]); cb["fdr"]=multipletests(cb["p"],method="fdr_bh")[1]
print("\n=== [敏感性] ComBat谐和后 UPDRS关联 ===")
print(f"  FDR<0.05: {(cb.fdr<0.05).sum()}个 -> {cb[cb.fdr<0.05]['name'].tolist()}")

df.to_csv(f"{BASE}/analysis_merged_final.csv",index=False,encoding="utf-8-sig")
print("\n最终结果表已保存 results_final_*.csv + analysis_merged_final.csv")
