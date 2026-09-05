# -*- coding: utf-8 -*-
"""补齐剩余图: Fig1(设计+验证), FigS1(复刻验证), FigS4(ComBat敏感性), Fig3(PD-vs-HC森林)。"""
import os; os.environ["MPLCONFIGDIR"]=r"G:/Temp/mpl"
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from neuroCombat import neuroCombat
import warnings; warnings.filterwarnings("ignore")

plt.rcParams['font.family']='sans-serif'; plt.rcParams['font.sans-serif']=['Arial','DejaVu Sans']
plt.rcParams['svg.fonttype']='none'; plt.rcParams['pdf.fonttype']=42
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':0.9})
P={'blue':'#0F4D92','red':'#B64342','grey':'#9A9A9A','teal':'#42949E','dark':'#272727','blue_soft':'#3775BA'}
BASE=r"F:/PD_brainage_data"; OUT=f"{BASE}/figures_nature"
REG=["sns","fpn","drs","vnt","dft","slt","lng","adt","vsl","lmb"]
EN={"sns":"Sensorimotor","fpn":"Frontoparietal","drs":"Dorsal attention","vnt":"Ventral attention","dft":"Default mode",
    "slt":"Salience","lng":"Language","adt":"Auditory","vsl":"Visual","lmb":"Limbic","global":"Global"}
def star(p): return '***' if p<0.001 else '**' if p<0.01 else '*' if p<0.05 else ''
def save(fig,name):
    for f in ['png','pdf','svg']: fig.savefig(f"{OUT}/{name}.{f}",dpi=400 if f=='png' else None,bbox_inches='tight')
    plt.close(fig); print(name,"done")

ba=pd.read_csv(f"{BASE}/brain_ages/all_brain_ages_final.csv"); ba["PATNO"]=ba["subject"].str.extract(r"(\d+)").astype(int)
cl=pd.read_csv(f"{BASE}/master_clinical.csv"); sc=pd.read_csv(f"{BASE}/scanner_info.csv")[["PATNO","Manuf","scanner"]]
df=ba.merge(cl,on="PATNO").merge(sc,on="PATNO").dropna(subset=["Manuf","age_at_visit","SEX"]).reset_index(drop=True)
hc=df[df.group=="HC"]

# ================= Fig 1 =================
fig=plt.figure(figsize=(7.0,6.0)); gs=fig.add_gridspec(2,1,height_ratios=[0.85,1.3],hspace=0.30)
axa=fig.add_subplot(gs[0]); axa.axis('off'); axa.set_xlim(0,10); axa.set_ylim(0,3)
steps=["T1-weighted\nMRI","CIVET cortical\nfeatures","Graph conv.\nnetwork","10 network\nbrain ages","Regional\nbrain-PAD"]
cols=[P['grey'],P['teal'],P['blue'],P['blue'],P['red']]; xs=np.linspace(1.1,8.9,5)
for x,s,c in zip(xs,steps,cols):
    axa.add_patch(FancyBboxPatch((x-0.78,1.05),1.56,1.15,boxstyle="round,pad=0.06",fc='white',ec=c,lw=1.8))
    axa.text(x,1.62,s,ha='center',va='center',fontsize=7.3,color=P['dark'])
for x0,x1 in zip(xs[:-1],xs[1:]):
    axa.annotate('',xy=(x1-0.82,1.62),xytext=(x0+0.82,1.62),arrowprops=dict(arrowstyle='-|>',color=P['dark'],lw=1.4))
axa.text(5,2.75,'n = 542 PD + 261 HC   (PPMI, multi-site; 1.5 T and 3 T)',ha='center',fontsize=8.5,style='italic',color=P['dark'])
axa.text(0.0,1.0,'A',transform=axa.transAxes,fontsize=16,fontweight='bold',va='top')
axb=fig.add_subplot(gs[1]); x=hc["age_at_visit"].values; y=hc["global"].values
a,b=np.polyfit(x,y,1); y_c=y-(a*x+b)+x; r=np.corrcoef(x,y)[0,1]; mae=np.abs(y_c-x).mean()
axb.scatter(x,y_c,s=15,alpha=0.5,color=P['blue_soft'],edgecolor='none')
lim=[x.min()-3,x.max()+3]; axb.plot(lim,lim,color=P['dark'],ls='--',lw=1); axb.set_xlim(lim); axb.set_ylim(lim)
axb.set_xlabel("Chronological age (years)"); axb.set_ylabel("Predicted brain age (years)")
axb.text(0.04,0.95,f"HC  n = {len(hc)}\nr = {r:.2f}\nMAE = {mae:.2f} years",transform=axb.transAxes,va='top',fontsize=9)
axb.text(-0.12,1.0,'B',transform=axb.transAxes,fontsize=16,fontweight='bold',va='top')
axb.set_title("Brain-age model validation in healthy controls",fontsize=9.5,loc='left',fontweight='bold')
save(fig,"Figure_1_design_validation")

# ================= Fig S1: 特征复刻验证 =================
th=[]; gw=[]
for i in df["PATNO"]:
    f=f"{BASE}/features_v2/{i}_features_20k.txt"
    if os.path.exists(f):
        arr=np.loadtxt(f); th.append(arr[:,0].mean()); gw.append(arr[:,1].mean())
fig,ax=plt.subplots(1,2,figsize=(7.0,3.0))
ax[0].hist(th,bins=30,color=P['teal'],edgecolor='white',lw=0.4); ax[0].axvline(3.12,color=P['red'],ls='--',lw=1.6)
ax[0].set_xlabel("Mean cortical thickness (mm)"); ax[0].set_ylabel("Subjects")
ax[0].text(0.97,0.92,"model target\n3.12 mm",transform=ax[0].transAxes,ha='right',va='top',color=P['red'],fontsize=8)
ax[0].text(-0.14,1.02,'A',transform=ax[0].transAxes,fontsize=15,fontweight='bold')
ax[1].hist(gw,bins=30,color=P['blue_soft'],edgecolor='white',lw=0.4); ax[1].axvline(-0.257,color=P['red'],ls='--',lw=1.6)
ax[1].set_xlabel("Mean GM/WM feature  (log ratio)"); ax[1].set_ylabel("Subjects")
ax[1].text(0.05,0.92,"reference\n−0.257",transform=ax[1].transAxes,ha='left',va='top',color=P['red'],fontsize=8)
ax[1].text(-0.14,1.02,'B',transform=ax[1].transAxes,fontsize=15,fontweight='bold')
fig.suptitle("Cortical feature reconstruction matches expected distributions",fontsize=9.5,fontweight='bold',y=1.02)
save(fig,"Figure_S1_reconstruction_validation")

# ================= Fig 3: PD vs HC 森林图 =================
PVH=pd.read_csv(f"{BASE}/results_final_pdvshc.csv").sort_values('beta').reset_index(drop=True)
fig,axc=plt.subplots(figsize=(6.2,4.0))
for i in range(len(PVH)):
    sig=PVH['fdr'][i]<0.05; col=P['blue'] if sig else P['grey']
    axc.plot([PVH['lo'][i],PVH['hi'][i]],[i,i],color=col,lw=1.7,solid_capstyle='round')
    axc.scatter(PVH['beta'][i],i,s=36 if sig else 22,color=col,zorder=3,edgecolor='white',linewidth=0.6)
    if star(PVH['fdr'][i]): axc.text(PVH['lo'][i]-0.08,i,star(PVH['fdr'][i]),va='center',ha='right',color=P['blue'],fontsize=10,fontweight='bold')
axc.axvline(0,color=P['grey'],ls=(0,(4,3)),lw=0.9)
axc.set_yticks(range(len(PVH))); axc.set_yticklabels(PVH['name']); axc.tick_params(axis='y',length=0)
axc.set_xlabel("Regional brain-PAD difference  (PD − HC, years)")
axc.set_title("PD vs HC regional brain age  (LME + scanner random effect)\nnegative = younger-appearing in PD; blue = FDR<0.05",fontsize=9,loc='left',fontweight='bold')
save(fig,"Figure_3_PD_vs_HC_forest")

# ================= Fig S4: ComBat 敏感性 =================
df2=df.copy(); df2["batch"]=df2["scanner"]; vc=df2["batch"].value_counts()
df2.loc[df2["batch"].isin(vc[vc<10].index),"batch"]=df2["Manuf"].astype(str)+"_other"
ALLc=REG+["global","AD"]; cov=df2[["batch","age_at_visit","SEX"]].copy(); cov["g"]=(df2.group=="PD").astype(int)
H=neuroCombat(dat=df2[ALLc].T.values,covars=cov,batch_col="batch",categorical_cols=["SEX","g"],continuous_cols=["age_at_visit"])["data"].T
for j,c in enumerate(ALLc): df2[c+"_h"]=H[:,j]
hc2=df2[df2.group=="HC"]; pdf2=df2[df2.group=="PD"]
U=pd.read_csv(f"{BASE}/results_final_updrs.csv").set_index("r")
rows=[]
for c in REG+["global"]:
    a,b=np.polyfit(hc2["age_at_visit"],hc2[c+"_h"],1); pad=pdf2[c+"_h"]-(a*pdf2["age_at_visit"]+b)
    d=pdf2.assign(z=(pad-pad.mean())/pad.std()).dropna(subset=["z","updrs3_score","age_at_visit","SEX","EDUCYRS"])
    m=smf.ols("updrs3_score ~ z + age_at_visit + C(SEX) + EDUCYRS",d).fit()
    rows.append([c,EN[c],U.loc[c,"beta"],U.loc[c,"fdr"],m.params["z"],m.pvalues["z"]])
S=pd.DataFrame(rows,columns=["r","name","b_lme","fdr_lme","b_cb","p_cb"]); S["fdr_cb"]=multipletests(S["p_cb"],method="fdr_bh")[1]
S=S.sort_values("b_lme").reset_index(drop=True)
fig,ax=plt.subplots(figsize=(6.4,4.2))
for i in range(len(S)):
    ax.plot([S['b_lme'][i],S['b_cb'][i]],[i,i],color=P['grey'],lw=1.2,zorder=1)
    ax.scatter(S['b_lme'][i],i,s=40,color=P['red'] if S['fdr_lme'][i]<0.05 else P['grey'],zorder=3,edgecolor='white',lw=0.6)
    ax.scatter(S['b_cb'][i],i,s=40,marker='D',color=P['blue'] if S['fdr_cb'][i]<0.05 else '#C8C8C8',zorder=3,edgecolor='white',lw=0.6)
ax.axvline(0,color=P['grey'],ls=(0,(4,3)),lw=0.9)
ax.set_yticks(range(len(S))); ax.set_yticklabels(S['name']); ax.tick_params(axis='y',length=0)
ax.set_xlabel("Standardized β  (regional brain-PAD → UPDRS-III)")
ax.scatter([],[],color=P['red'],label='LME + scanner random effect (FDR<0.05)')
ax.scatter([],[],marker='D',color=P['blue'],label='ComBat harmonized')
ax.legend(loc='lower right',fontsize=7.5)
ax.set_title("Sensitivity to scanner harmonization",fontsize=9.5,loc='left',fontweight='bold')
save(fig,"Figure_S4_combat_sensitivity")
print("ALL DONE")
