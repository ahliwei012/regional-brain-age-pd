# -*- coding: utf-8 -*-
"""PD特异新分析: 用10网络区域脑龄谱做影像驱动亚型划分。
   A) 原始PAD聚类(严重度轴)  B) 去个体均值的"谱形"聚类(空间模式轴)
   然后检验亚型间临床差异 + 扫描仪混杂检查。"""
import numpy as np, pandas as pd
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
pdf=df[df.group=="PD"].copy()
padc=[r+"_PAD" for r in REG]
pdf=pdf.dropna(subset=padc).reset_index(drop=True)
print(f"PD 用于聚类: n={len(pdf)}")

X_raw=pdf[padc].values
X_shape=X_raw-X_raw.mean(axis=1,keepdims=True)   # 去个体总体水平 -> 只剩空间谱形

def pick_k(X,name,kmax=5):
    Z=StandardScaler().fit_transform(X)
    print(f"\n--- {name}: 选k (silhouette) ---")
    best=(None,-1)
    for k in range(2,kmax+1):
        lab=KMeans(n_clusters=k,n_init=50,random_state=2026).fit_predict(Z)
        s=silhouette_score(Z,lab)
        print(f"  k={k}: silhouette={s:.3f}")
        if s>best[1]: best=(k,s)
    return best[0],Z

def clin_compare(pdf,lab,tag):
    pdf=pdf.copy(); pdf['cl']=lab
    ks=sorted(pdf.cl.unique())
    print(f"\n===== {tag}: 亚型临床特征 (n per cluster: {[int((lab==k).sum()) for k in ks]}) =====")
    rows=[]
    for var,nm in [("age_at_visit","Age"),("EDUCYRS","Education"),("duration_yrs","Duration"),
                   ("NHY","H&Y"),("updrs3_score","UPDRS-III"),("moca","MoCA"),("LEDD","LEDD"),
                   ("global_PAD","Global brain-PAD")]:
        if var not in pdf.columns: continue
        groups=[pdf[pdf.cl==k][var].dropna() for k in ks]
        groups=[g for g in groups if len(g)>2]
        if len(groups)<2: continue
        F,p=stats.f_oneway(*groups)
        means=[f"{pdf[pdf.cl==k][var].mean():.2f}" for k in ks]
        rows.append([nm]+means+[p])
    R=pd.DataFrame(rows,columns=["Variable"]+[f"C{k+1}" for k in ks]+["p"])
    R["fdr"]=multipletests(R["p"],method="fdr_bh")[1]
    print(R.round(4).to_string(index=False))
    # 性别/扫描仪混杂
    ct=pd.crosstab(pdf.cl,pdf.SEX); _,psex,_,_=stats.chi2_contingency(ct)
    ctb=pd.crosstab(pdf.cl,pdf.batch); _,pbatch,_,_=stats.chi2_contingency(ctb)
    print(f"  性别 chi2 p={psex:.3f} | 扫描仪批次 chi2 p={pbatch:.4f}  <-- 若显著说明亚型受站点混杂")
    # 各网络脑龄谱
    prof=pdf.groupby('cl')[padc].mean().T
    prof.index=[EN[r] for r in REG]
    prof.columns=[f"C{k+1}" for k in ks]
    print("\n  各亚型区域brain-PAD均值(年):")
    print(prof.round(2).to_string())
    return R,prof,pdf

# ---------- A) 原始PAD ----------
kA,ZA=pick_k(X_raw,"A 原始brain-PAD")
labA=KMeans(n_clusters=kA,n_init=50,random_state=2026).fit_predict(ZA)
RA,profA,pdfA=clin_compare(pdf,labA,f"A 原始brain-PAD 聚类 (k={kA})")

# ---------- B) 谱形(去个体均值) ----------
kB,ZB=pick_k(X_shape,"B 谱形(去个体均值)")
labB=KMeans(n_clusters=kB,n_init=50,random_state=2026).fit_predict(ZB)
RB,profB,pdfB=clin_compare(pdf,labB,f"B 谱形 聚类 (k={kB})")

pdfA[['PATNO','cl']].rename(columns={'cl':'subtype_raw'}).to_csv(f"{BASE}/pd_subtype_raw.csv",index=False)
pdfB[['PATNO','cl']].rename(columns={'cl':'subtype_shape'}).to_csv(f"{BASE}/pd_subtype_shape.csv",index=False)
profA.to_csv(f"{BASE}/pd_subtype_profile_raw.csv"); profB.to_csv(f"{BASE}/pd_subtype_profile_shape.csv")
print("\nDONE")
