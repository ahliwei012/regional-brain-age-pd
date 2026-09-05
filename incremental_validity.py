# -*- coding: utf-8 -*-
"""Incremental validity and paradox checks:
   (1) 增量效度: 控制同网络原始平均皮层厚度后, 区域 brain-PAD 是否仍预测临床?
   (2) 悖论判别: '更年轻'网络(sns/lmb)的原始厚度是否 PD<HC? (判 模型伪迹 vs 生物学)
   附: 置换检验空基线。"""
import os, glob, numpy as np, pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import warnings; warnings.filterwarnings("ignore")
np.random.seed(2026)
BASE=r"F:/PD_brainage_data"; REPO=r"F:/Deep_learning_derived_ MRI regional brain age/regional_Brain_age-main"
EN={"sns":"Sensorimotor","fpn":"Frontoparietal","drs":"Dorsal attention","vnt":"Ventral attention","dft":"Default mode",
    "slt":"Salience","lng":"Language","adt":"Auditory","vsl":"Visual","lmb":"Limbic"}
# 网络->AAL区域 (取自模型源码)
NET={"sns":[1,2,19,20,57,58,69,70],"fpn":[5,6,7,8,9,10,65,66],"drs":[3,4,59,60],
     "vnt":[11,12,13,14,63,64],"dft":[21,22,25,26,27,28,35,36,65,66,67,68,85,86],
     "slt":[29,30,31,32],"lng":[11,12,13,17,63],"adt":[79,80,81,82],
     "vsl":[43,44,45,46,47,48,49,50,51,52,53,54,55,56,89,90],
     "lmb":[15,16,23,24,33,34,39,40,83,84,87,88]}
REG=list(NET.keys())

aal=pd.read_table(f"{REPO}/surface_information/aal_atlas_20k.txt",header=None).values.ravel().astype(int)
assert aal.shape[0]==20484, aal.shape
masks={r:np.isin(aal,NET[r]) for r in REG}
for r in REG: assert masks[r].sum()>0

df=pd.read_csv(f"{BASE}/analysis_merged_final.csv")
pats=df["PATNO"].astype(int).tolist()
# 逐被试算各网络原始平均厚度(col0)
rows={r+"_thick":[] for r in REG}; keep=[]
for p in pats:
    f=f"{BASE}/features_v2/{p}_features_20k.txt"
    if not os.path.exists(f):
        for r in REG: rows[r+"_thick"].append(np.nan)
        continue
    a=pd.read_csv(f,header=None,sep=r"\s+").values
    th=a[:,0]
    for r in REG: rows[r+"_thick"].append(float(th[masks[r]].mean()))
for k,v in rows.items(): df[k]=v
print("各网络原始平均厚度已计算。样例(全样本均值):")
for r in REG: print(f"  {EN[r]:16s} {df[r+'_thick'].mean():.3f} mm")

pdf=df[df.group=="PD"].copy()

def lme(data,outcome,terms):
    d=data.dropna(subset=[outcome,"age_at_visit","SEX","EDUCYRS","batch"]+[t for t in terms if t in data.columns]).copy()
    f=f"{outcome} ~ {' + '.join(terms)} + age_at_visit + C(SEX) + EDUCYRS"
    m=smf.mixedlm(f,d,groups=d["batch"]).fit()
    return m,len(d)

print("\n########## G1 增量效度: 控制原始厚度后 brain-PAD 是否仍显著 ##########")
SIG={"updrs3_score":["vsl","vnt","slt","lng"],"moca":["fpn","vnt"]}
g1=[]
for out,nets in SIG.items():
    print(f"\n=== {out} ===")
    for r in nets:
        pad=r+"_z"; thick=r+"_thick"
        # 模型A: 仅brain-PAD; 模型B: brain-PAD + 原始厚度
        mA,nA=lme(pdf,out,[pad])
        mB,nB=lme(pdf,out,[pad,thick])
        bA,pA=mA.params[pad],mA.pvalues[pad]
        bB,pB=mB.params[pad],mB.pvalues[pad]
        bt,pt=mB.params[thick],mB.pvalues[thick]
        drop=100*(1-abs(bB)/abs(bA)) if bA!=0 else np.nan
        print(f"  {EN[r]:16s} PAD单独 b={bA:+.3f}(p={pA:.3f}) | 加厚度后 PAD b={bB:+.3f}(p={pB:.3f}) 厚度b={bt:+.3f}(p={pt:.3f}) | PAD效应衰减{drop:.0f}%")
        g1.append([out,EN[r],bA,pA,bB,pB,bt,pt,drop,nB])
G1=pd.DataFrame(g1,columns=["outcome","network","b_PADonly","p_PADonly","b_PADadj","p_PADadj","b_thick","p_thick","PAD_atten_pct","n"])
G1.to_csv(f"{BASE}/G1_incremental_validity.csv",index=False,encoding="utf-8-sig")

print("\n########## 置换检验空基线 (打乱brain-PAD, 1000次, 看|b|>=观测的概率) ##########")
for out,r in [("updrs3_score","vsl"),("updrs3_score","vnt"),("moca","fpn")]:
    d=pdf.dropna(subset=[out,r+"_z","age_at_visit","SEX","EDUCYRS","batch"]).copy()
    m,_=lme(d,out,[r+"_z"]); obs=abs(m.params[r+"_z"])
    null=[]
    for i in range(1000):
        dd=d.copy(); dd[r+"_z"]=np.random.permutation(dd[r+"_z"].values)
        mm=smf.mixedlm(f"{out} ~ {r}_z + age_at_visit + C(SEX) + EDUCYRS",dd,groups=dd["batch"]).fit()
        null.append(abs(mm.params[r+"_z"]))
    pperm=(np.sum(np.array(null)>=obs)+1)/1001
    print(f"  {out} / {EN[r]:16s} 观测|b|={obs:.3f}, 置换p={pperm:.4f}")

print("\n########## G2 悖论判别: 原始平均厚度 PD vs HC ##########")
dfz=df.copy(); dfz["grp"]=(dfz.group=="PD").astype(int)
g2=[]
for r in REG:
    d=dfz.dropna(subset=[r+"_thick","age_at_visit","SEX","EDUCYRS","batch"]).copy()
    m=smf.mixedlm(f"{r}_thick ~ grp + age_at_visit + C(SEX) + EDUCYRS",d,groups=d["batch"]).fit()
    b,p=m.params["grp"],m.pvalues["grp"]
    padb=pd.read_csv(f"{BASE}/results_final_pdvshc.csv").set_index("r").loc[r,"beta"]
    g2.append([EN[r],b,p,padb])
G2=pd.DataFrame(g2,columns=["network","thick_PDminusHC_mm","p","brainPAD_PDminusHC_yr"])
G2["thick_fdr"]=multipletests(G2["p"],method="fdr_bh")[1]
print(G2.round(4).to_string(index=False))
print("\n判读: 若某网络 brain-PAD显示PD更年轻(负) 而 原始厚度PD并不更薄(thick差≥0或n.s.) -> 与厚度一致,非纯伪迹;")
print("      若 brain-PAD更年轻 但 厚度显著更薄(负) -> brain-PAD与厚度矛盾 -> 提示模型低估(伪迹)。")
G2.to_csv(f"{BASE}/G2_thickness_pdvshc.csv",index=False,encoding="utf-8-sig")
df.to_csv(f"{BASE}/analysis_with_thickness.csv",index=False,encoding="utf-8-sig")
print("\nDONE. 输出: G1_incremental_validity.csv / G2_thickness_pdvshc.csv / analysis_with_thickness.csv")
