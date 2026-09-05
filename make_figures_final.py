# -*- coding: utf-8 -*-
"""最终图(550 PD + 261 HC, LME+站点): 森林图/热图/PD-vs-HC/散点。"""
import pandas as pd, numpy as np, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
import warnings; warnings.filterwarnings("ignore")

BASE=r"F:/PD_brainage_data"; FIG=f"{BASE}/figures_final"; os.makedirs(FIG,exist_ok=True)
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.spines.top":False,
                     "axes.spines.right":False,"figure.dpi":300,"savefig.dpi":300})
def star(p): return "***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else ""
U=pd.read_csv(f"{BASE}/results_final_updrs.csv"); M=pd.read_csv(f"{BASE}/results_final_moca.csv")
PVH=pd.read_csv(f"{BASE}/results_final_pdvshc.csv"); df=pd.read_csv(f"{BASE}/analysis_merged_final.csv")

# Fig1 森林图 UPDRS
fig,ax=plt.subplots(figsize=(6,4.2)); s=U.sort_values("beta").reset_index(drop=True)
for i in range(len(s)):
    c="#c0392b" if s["fdr"][i]<0.05 else "#95a5a6"
    ax.plot([s["lo"][i],s["hi"][i]],[i,i],color=c,lw=2); ax.scatter(s["beta"][i],i,color=c,s=45,zorder=2)
    if star(s["fdr"][i]): ax.text(s["hi"][i]+0.03,i,star(s["fdr"][i]),va="center",color="#c0392b",fontsize=11)
ax.axvline(0,color="k",lw=0.8,ls="--"); ax.set_yticks(range(len(s))); ax.set_yticklabels(s["name"])
ax.set_xlabel("Standardized β (regional brain-PAD → UPDRS-III)")
ax.set_title("Regional cortical brain age vs motor severity in PD\n(LME + scanner random effect; n=550 PD, 261 HC; red=FDR<0.05)",fontsize=9.5)
plt.tight_layout(); plt.savefig(f"{FIG}/Fig1_forest_updrs.png"); plt.close()

# Fig2 热图 区域×{UPDRS,MoCA}
mat=np.column_stack([U.set_index("r").loc[[c for c in U["r"]]]["beta"].values, M.set_index("r").loc[[c for c in U["r"]]]["beta"].values])
names=U["name"].values; uf=U["fdr"].values; mf=M.set_index("r").loc[U["r"]]["fdr"].values
fig,ax=plt.subplots(figsize=(4.3,5))
norm=TwoSlopeNorm(vcenter=0,vmin=-abs(mat).max(),vmax=abs(mat).max())
im=ax.imshow(mat,cmap="RdBu_r",norm=norm,aspect="auto")
ax.set_xticks([0,1]); ax.set_xticklabels(["UPDRS-III\n(motor)","MoCA\n(cognition)"])
ax.set_yticks(range(len(names))); ax.set_yticklabels(names)
for i in range(len(names)):
    ax.text(0,i,star(uf[i]),ha="center",va="center",fontsize=12); ax.text(1,i,star(mf[i]),ha="center",va="center",fontsize=12)
ax.set_title("Regional brain-PAD ↔ outcomes\n(std β; LME+site; * FDR<0.05)",fontsize=9.5)
plt.colorbar(im,ax=ax,fraction=0.046,pad=0.04,label="Standardized β")
plt.tight_layout(); plt.savefig(f"{FIG}/Fig2_heatmap.png"); plt.close()

# Fig3 PD vs HC 森林图
fig,ax=plt.subplots(figsize=(6,4.2)); s=PVH.sort_values("beta").reset_index(drop=True)
for i in range(len(s)):
    c="#2471a3" if s["fdr"][i]<0.05 else "#95a5a6"
    ax.plot([s["lo"][i],s["hi"][i]],[i,i],color=c,lw=2); ax.scatter(s["beta"][i],i,color=c,s=45,zorder=2)
    if star(s["fdr"][i]): ax.text(s["hi"][i]+0.05,i,star(s["fdr"][i]),va="center",color="#2471a3",fontsize=11)
ax.axvline(0,color="k",lw=0.8,ls="--"); ax.set_yticks(range(len(s))); ax.set_yticklabels(s["name"])
ax.set_xlabel("Regional brain-PAD difference (PD − HC, years)")
ax.set_title("PD vs HC regional brain age\n(LME+site; blue=FDR<0.05; negative=younger in PD)",fontsize=9.5)
plt.tight_layout(); plt.savefig(f"{FIG}/Fig3_pdvshc.png"); plt.close()

# Fig4 RF重要性
REG=["sns","fpn","drs","vnt","dft","slt","lng","adt","vsl","lmb"]
EN={"sns":"Sensorimotor","fpn":"Frontoparietal","drs":"Dorsal att","vnt":"Ventral att","dft":"Default",
    "slt":"Salience","lng":"Language","adt":"Auditory","vsl":"Visual","lmb":"Limbic"}
pdf=df[df.group=="PD"]; fc=[c+"_PAD" for c in REG]
d=pdf.dropna(subset=fc+["updrs3_score"]); thr=d["updrs3_score"].median()
X=d[fc].values; y=(d["updrs3_score"]>thr).astype(int).values
imp=np.zeros((500,len(fc)))
for i in range(500):
    Xtr,_,ytr,_=train_test_split(X,y,test_size=0.2,random_state=i,stratify=y)
    imp[i]=RandomForestClassifier(n_estimators=200,random_state=i,n_jobs=-1).fit(Xtr,ytr).feature_importances_
rf=pd.DataFrame({"name":[EN[c] for c in REG],"imp":imp.mean(0),"sd":imp.std(0)}).sort_values("imp")
fig,ax=plt.subplots(figsize=(6,4.2))
ax.barh(rf["name"],rf["imp"],xerr=rf["sd"],color="#aed6f1",edgecolor="#2980b9",error_kw=dict(ecolor="#34495e",capsize=2))
ax.set_xlabel("Mean feature importance (Gini, 500 bootstraps)")
ax.set_title("Regional brain-PAD importance for poor motor outcome (PD)",fontsize=9.5)
plt.tight_layout(); plt.savefig(f"{FIG}/Fig4_rf.png"); plt.close()

print("最终图已生成:", os.listdir(FIG))
