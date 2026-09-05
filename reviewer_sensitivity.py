# -*- coding: utf-8 -*-
"""Sensitivity analyses:
   (1) 临床关联额外校正 病程(duration_yrs)+LEDD;
   (2) 跨全部检验的多重比较(pooled FDR / Bonferroni);
   (3) younger悖论: PD内 brain-PAD符号与分布;
   (4) 明确 β 单位(每SD)并给出每年等价。"""
import pandas as pd, numpy as np
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import warnings; warnings.filterwarnings("ignore")
BASE=r"F:/PD_brainage_data"
REG=["sns","fpn","drs","vnt","dft","slt","lng","adt","vsl","lmb"]
EN={"sns":"Sensorimotor","fpn":"Frontoparietal","drs":"Dorsal attention","vnt":"Ventral attention","dft":"Default mode",
    "slt":"Salience","lng":"Language","adt":"Auditory","vsl":"Visual","lmb":"Limbic"}
df=pd.read_csv(f"{BASE}/analysis_merged_final.csv")
pdf=df[df.group=="PD"].copy()

def lme(data,outcome,pred,extra=""):
    d=data.dropna(subset=[pred,outcome,"age_at_visit","SEX","EDUCYRS","batch"]+
                  ([c for c in extra.replace('+',' ').split() if c in data.columns]) ).copy().rename(columns={pred:"x"})
    f=f"{outcome} ~ x + age_at_visit + C(SEX) + EDUCYRS{extra}"
    m=smf.mixedlm(f, d, groups=d["batch"]).fit()
    ci=m.conf_int().loc["x"]
    return m.params["x"], ci[0], ci[1], m.pvalues["x"], len(d)

print("### (1) 加入 病程+LEDD 校正后的敏感性 (标准化brain-PAD为预测)")
for tag,out in [("UPDRS-III","updrs3_score"),("MoCA","moca")]:
    rows=[]
    for c in REG:
        b,lo,hi,p,n=lme(pdf,out,c+"_z",extra=" + duration_yrs + LEDD")
        rows.append([EN[c],b,lo,hi,p,n])
    r=pd.DataFrame(rows,columns=["name","beta","lo","hi","p","n"]); r["fdr"]=multipletests(r["p"],method="fdr_bh")[1]
    sig=r[r.fdr<0.05]["name"].tolist()
    print(f"  {tag}: n={r['n'].iloc[0]}, FDR<0.05 -> {sig}")
    print(r.round(3)[["name","beta","p","fdr"]].to_string(index=False))
    r.to_csv(f"{BASE}/sens_{out}_durLEDD.csv",index=False,encoding="utf-8-sig")

print("\n### (2) 跨全部检验 pooled 多重比较 (10网络 x 3结局 = 30 检验)")
allp=[]
for tag in ["updrs","moca","pdvshc"]:
    t=pd.read_csv(f"{BASE}/results_final_{tag}.csv"); t=t[t.r.isin(REG)]
    for _,row in t.iterrows(): allp.append([tag,row["name"],row["p"]])
A=pd.DataFrame(allp,columns=["family","name","p"])
A["fdr_pooled"]=multipletests(A["p"],method="fdr_bh")[1]
A["bonf"]=multipletests(A["p"],method="bonferroni")[1]
print("  存活 pooled-FDR<0.05:")
print(A[A.fdr_pooled<0.05].round(4).to_string(index=False))
print("  存活 Bonferroni<0.05:", A[A.bonf<0.05][["family","name"]].values.tolist())
A.to_csv(f"{BASE}/sens_pooled_multipletest.csv",index=False,encoding="utf-8-sig")

print("\n### (3) younger悖论: PD组各网络 brain-PAD 的均值(年)")
for c in REG:
    v=pdf[c+"_PAD"].dropna()
    print(f"  {EN[c]:16s} mean {v.mean():+.2f}y  median {v.median():+.2f}y  %>0: {100*(v>0).mean():.0f}%")

print("\n### (4) brain-PAD 每SD对应多少年 (PD组内SD)")
for c in ["vsl","vnt","slt","lng","fpn"]:
    sd=pdf[c+"_PAD"].std()
    print(f"  {EN[c]:16s} 1 SD = {sd:.2f} years")

print("\n### (5) UPDRS OFF vs ON 一致性 (核心4网络, 用ON评分重跑)")
if "updrs3_score_on" in pdf.columns:
    for c in ["vsl","vnt","slt","lng"]:
        b,lo,hi,p,n=lme(pdf,"updrs3_score_on",c+"_z")
        print(f"  {EN[c]:16s} beta(ON)={b:+.2f} p={p:.3f} n={n}")
