# -*- coding: utf-8 -*-
"""用最终数据重做 S2(热图) 与 S3(随机森林)。"""
import os; os.environ["MPLCONFIGDIR"]=r"G:/Temp/mpl"
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
import warnings; warnings.filterwarnings("ignore")
plt.rcParams['font.family']='sans-serif'; plt.rcParams['font.sans-serif']=['Arial','DejaVu Sans']
plt.rcParams['svg.fonttype']='none'; plt.rcParams['pdf.fonttype']=42
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':0.9})
BASE=r"F:/PD_brainage_data"; OUT=f"{BASE}/figures_nature"
REG=["sns","fpn","drs","vnt","dft","slt","lng","adt","vsl","lmb"]
EN={"sns":"Sensorimotor","fpn":"Frontoparietal","drs":"Dorsal attention","vnt":"Ventral attention","dft":"Default mode",
    "slt":"Salience","lng":"Language","adt":"Auditory","vsl":"Visual","lmb":"Limbic"}
def star(p): return '***' if p<0.001 else '**' if p<0.01 else '*' if p<0.05 else ''
def save(fig,name):
    for f in ['png','pdf','svg']: fig.savefig(f"{OUT}/{name}.{f}",dpi=400 if f=='png' else None,bbox_inches='tight')
    plt.close(fig); print(name,"done")

U=pd.read_csv(f"{BASE}/results_final_updrs.csv").set_index("r"); M=pd.read_csv(f"{BASE}/results_final_moca.csv").set_index("r")
# S2 热图
fig,ax=plt.subplots(figsize=(4.3,5))
mat=np.column_stack([U.loc[REG,"beta"].values, M.loc[REG,"beta"].values]); vmax=np.abs(mat).max()
im=ax.imshow(mat,cmap='RdBu_r',norm=TwoSlopeNorm(0,-vmax,vmax),aspect='auto')
ax.set_xticks([0,1]); ax.set_xticklabels(['UPDRS-III','MoCA']); ax.set_yticks(range(len(REG))); ax.set_yticklabels([EN[c] for c in REG],fontsize=8)
ax.tick_params(length=0)
for i,c in enumerate(REG):
    if star(U.loc[c,"fdr"]): ax.text(0,i,star(U.loc[c,"fdr"]),ha='center',va='center',fontsize=11)
    if star(M.loc[c,"fdr"]): ax.text(1,i,star(M.loc[c,"fdr"]),ha='center',va='center',fontsize=11)
for sp in ax.spines.values(): sp.set_visible(False)
ax.set_xticks(np.arange(-.5,2,1),minor=True); ax.set_yticks(np.arange(-.5,len(REG),1),minor=True)
ax.grid(which='minor',color='white',lw=1.3); ax.tick_params(which='minor',length=0)
cb=fig.colorbar(im,ax=ax,fraction=0.05,pad=0.06); cb.set_label('Standardized β',fontsize=8)
ax.set_title("Motor vs cognitive dissociation\n(* FDR<0.05)",fontsize=9.5,loc='left',fontweight='bold')
save(fig,"Figure_S2_motor_cognitive_heatmap")

# S3 随机森林
df=pd.read_csv(f"{BASE}/analysis_merged_final.csv"); pdf=df[df.group=="PD"]
fc=[c+"_PAD" for c in REG]; d=pdf.dropna(subset=fc+["updrs3_score"]); thr=d["updrs3_score"].median()
X=d[fc].values; y=(d["updrs3_score"]>thr).astype(int).values
imp=np.zeros((1000,len(fc)))
for i in range(1000):
    Xtr,_,ytr,_=train_test_split(X,y,test_size=0.2,random_state=i,stratify=y)
    imp[i]=RandomForestClassifier(n_estimators=200,random_state=i,n_jobs=-1).fit(Xtr,ytr).feature_importances_
rf=pd.DataFrame({"name":[EN[c] for c in REG],"imp":imp.mean(0),"sd":imp.std(0)}).sort_values("imp")
fig,ax=plt.subplots(figsize=(6.2,4.0))
ax.barh(rf["name"],rf["imp"],xerr=rf["sd"],color="#AED6F1",edgecolor="#2980B9",error_kw=dict(ecolor="#34495E",capsize=2))
ax.set_xlabel("Mean feature importance (Gini, 1000 bootstraps)")
ax.set_title(f"Regional brain-PAD importance for poor motor outcome (PD, n={len(d)})",fontsize=9.5,loc='left',fontweight='bold')
save(fig,"Figure_S3_random_forest_importance")
print("DONE")
