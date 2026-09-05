# -*- coding: utf-8 -*-
"""严格版: 先把 age/sex/edu/scanner 残差化, 再对区域脑龄谱聚类。
   若去混杂后仍有临床可区分的亚型 -> 真新颖; 否则为诚实阴性。"""
import numpy as np, pandas as pd
import statsmodels.formula.api as smf
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler
from scipy import stats
from statsmodels.stats.multitest import multipletests
import warnings; warnings.filterwarnings("ignore")
np.random.seed(2026)
BASE=r"F:/PD_brainage_data"
REG=["sns","fpn","drs","vnt","dft","slt","lng","adt","vsl","lmb"]
EN={"sns":"Sensorimotor","fpn":"Frontoparietal","drs":"Dorsal attn","vnt":"Ventral attn","dft":"Default",
    "slt":"Salience","lng":"Language","adt":"Auditory","vsl":"Visual","lmb":"Limbic"}
df=pd.read_csv(f"{BASE}/analysis_with_thickness.csv")
pdf=df[df.group=="PD"].dropna(subset=[r+"_PAD" for r in REG]+["age_at_visit","SEX","EDUCYRS","batch"]).reset_index(drop=True)
print("PD n=",len(pdf))

# 1) 残差化: 去掉 age/sex/edu/scanner
res=np.zeros((len(pdf),len(REG)))
for j,r in enumerate(REG):
    m=smf.ols(f"{r}_PAD ~ age_at_visit + C(SEX) + EDUCYRS + C(batch)",data=pdf).fit()
    res[:,j]=m.resid
R=pd.DataFrame(res,columns=REG)
# 残差与年龄相关性检查
print("残差化后 各网络残差与年龄的|r|最大:",max(abs(np.corrcoef(R[c],pdf.age_at_visit)[0,1]) for c in REG).round(3))

Z=StandardScaler().fit_transform(R.values)
print("\n--- 去混杂后聚类结构 (silhouette) ---")
best=(None,-1)
for k in range(2,6):
    lab=KMeans(n_clusters=k,n_init=100,random_state=2026).fit_predict(Z)
    s=silhouette_score(Z,lab); print(f"  k={k}: {s:.3f}")
    if s>best[1]: best=(k,s)
k,sil=best
print(f"最佳 k={k}, silhouette={sil:.3f}  (>0.5佳, 0.25-0.5弱, <0.25基本无结构)")
lab=KMeans(n_clusters=k,n_init=100,random_state=2026).fit_predict(Z)
pdf['cl']=lab; ks=sorted(pdf.cl.unique())
print("每簇n:",[int((lab==kk).sum()) for kk in ks])

rows=[]
for var,nm in [("age_at_visit","Age"),("EDUCYRS","Education"),("duration_yrs","Duration"),
               ("NHY","H&Y"),("updrs3_score","UPDRS-III"),("moca","MoCA"),("LEDD","LEDD")]:
    g=[pdf[pdf.cl==kk][var].dropna() for kk in ks]; g=[x for x in g if len(x)>2]
    if len(g)<2: continue
    F,p=stats.f_oneway(*g)
    rows.append([nm]+[f"{pdf[pdf.cl==kk][var].mean():.2f}" for kk in ks]+[p])
T=pd.DataFrame(rows,columns=["Variable"]+[f"C{kk+1}" for kk in ks]+["p"])
T["fdr"]=multipletests(T["p"],method="fdr_bh")[1]
print("\n===== 去混杂后 亚型临床差异 =====")
print(T.round(4).to_string(index=False))
ct=pd.crosstab(pdf.cl,pdf.batch); _,pb,_,_=stats.chi2_contingency(ct)
cs=pd.crosstab(pdf.cl,pdf.SEX); _,psx,_,_=stats.chi2_contingency(cs)
print(f"  站点混杂 chi2 p={pb:.3f} | 性别 p={psx:.3f}")
prof=pd.DataFrame(res,columns=[EN[r] for r in REG]).assign(cl=lab).groupby('cl').mean().T
prof.columns=[f"C{kk+1}" for kk in ks]
print("\n各亚型 残差化区域脑龄谱:"); print(prof.round(2).to_string())
print("\n判读: 若 silhouette<0.25 且临床FDR全不显著 -> 区域脑龄谱不能定义PD亚型(阴性)")
