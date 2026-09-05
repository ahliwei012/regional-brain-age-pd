# -*- coding: utf-8 -*-
"""Supplementary Figures:
   S5 = 增量效度森林图(控制原始皮层厚度前/后 brain-PAD 的β对比)
   S6 = 纵向 time×brain-PAD 交互效应(基线脑龄预测进展速率, 阴性)"""
import os; os.environ["MPLCONFIGDIR"]=r"G:/Temp/mpl"
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
import warnings; warnings.filterwarnings("ignore")
plt.rcParams['font.family']='sans-serif'; plt.rcParams['font.sans-serif']=['Arial','DejaVu Sans']
plt.rcParams['svg.fonttype']='none'; plt.rcParams['pdf.fonttype']=42
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':0.9})
BASE=r"F:/PD_brainage_data"; OUT=f"{BASE}/manuscript_figures"
BLUE='#0F4D92'; ORANGE='#D55E00'; GREY='#9A9A9A'
def star(p): return '***' if p<0.001 else '**' if p<0.01 else '*' if p<0.05 else ''
def save(fig,name):
    for f in ['png','pdf','svg']: fig.savefig(f"{OUT}/{name}.{f}",dpi=400 if f=='png' else None,bbox_inches='tight')
    plt.close(fig); print(name,"done")

df=pd.read_csv(f"{BASE}/analysis_with_thickness.csv"); pdf=df[df.group=="PD"].copy()
def lme_ci(outcome,terms):
    d=pdf.dropna(subset=[outcome,"age_at_visit","SEX","EDUCYRS","batch"]+[t for t in terms if t in pdf.columns]).copy()
    m=smf.mixedlm(f"{outcome} ~ {' + '.join(terms)} + age_at_visit + C(SEX) + EDUCYRS",d,groups=d["batch"]).fit()
    ci=m.conf_int(); key=terms[0]
    return m.params[key],ci.loc[key,0],ci.loc[key,1],m.pvalues[key]

# ---------- S5 增量效度 ----------
ITEMS=[("updrs3_score","vsl","Visual"),("updrs3_score","vnt","Ventral attention"),
       ("updrs3_score","slt","Salience"),("updrs3_score","lng","Language"),
       ("moca","fpn","Frontoparietal"),("moca","vnt","Ventral attention")]
fig,axes=plt.subplots(1,2,figsize=(9.2,3.6),gridspec_kw={'width_ratios':[2,1],'wspace':0.45})
for ax,(out,lab,rng,ttl) in zip(axes,[("updrs3_score","UPDRS-III  (points / SD)",range(0,4),"A  Motor severity"),
                                        ("moca","MoCA  (points / SD)",range(4,6),"B  Global cognition")]):
    rows=list(rng); ypos=np.arange(len(rows))[::-1]
    for yi,idx in zip(ypos,rows):
        o,r,name=ITEMS[idx]
        bo,lo,ho,po=lme_ci(o,[r+"_z"])            # PAD only
        ba,la,ha,pa=lme_ci(o,[r+"_z",r+"_thick"]) # PAD | thickness
        _,_,_,pt=lme_ci(o,[r+"_thick",r+"_z"])    # thickness term p (first term)
        ax.plot([lo,ho],[yi+0.16,yi+0.16],color=BLUE,lw=1.8,solid_capstyle='round')
        ax.scatter(bo,yi+0.16,s=42,color=BLUE,zorder=3,edgecolor='white',lw=0.6)
        ax.plot([la,ha],[yi-0.16,yi-0.16],color=ORANGE,lw=1.8,solid_capstyle='round')
        ax.scatter(ba,yi-0.16,s=42,color=ORANGE,marker='D',zorder=3,edgecolor='white',lw=0.6)
        ax.text(1.02,yi,f"thickness {'n.s.' if pt>=0.05 else star(pt)}",transform=ax.get_yaxis_transform(),
                fontsize=7.5,va='center',ha='left',color='black')
    ax.axvline(0,color=GREY,ls=(0,(4,3)),lw=0.9)
    ax.set_yticks(ypos); ax.set_yticklabels([ITEMS[i][2] for i in rows],fontsize=8.5)
    ax.set_ylim(-0.6,len(rows)-0.4); ax.set_xlabel(out if False else "")
    ax.set_xlabel([("UPDRS-III  β  (points / SD)") if out=="updrs3_score" else "MoCA  β  (points / SD)"][0],fontsize=8.5)
    ax.set_title(ttl,fontsize=9.5,loc='left',fontweight='bold',pad=8,color='black')
    ax.tick_params(colors='black'); ax.xaxis.label.set_color('black'); ax.yaxis.label.set_color('black')
from matplotlib.lines import Line2D
leg=[Line2D([0],[0],color=BLUE,marker='o',lw=1.8,ls='-',label='brain-PAD only'),
     Line2D([0],[0],color=ORANGE,marker='D',lw=1.8,ls='-',label='brain-PAD | cortical thickness')]
fig.legend(handles=leg,loc='upper right',bbox_to_anchor=(0.99,1.02),ncol=2,fontsize=8,frameon=False,handletextpad=0.3,columnspacing=1.4,labelcolor='black')
fig.suptitle("Regional brain age carries clinical information beyond cortical thickness",
             fontsize=10.5,fontweight='bold',x=0.02,ha='left',y=1.03,color='black')
save(fig,"Figure_4_Incremental_Validity")

# ---------- S6 纵向 time×PAD ----------
U=pd.read_csv(f"{BASE}/long_updrs_progression.csv"); M=pd.read_csv(f"{BASE}/long_moca_progression.csv")
fig,axes=plt.subplots(1,2,figsize=(9.2,4.0),sharey=True,gridspec_kw={'wspace':0.08})
for ax,(D,ttl,xl) in zip(axes,[(U,"A  Motor progression (UPDRS-III)","β  time × brain-PAD  (points per SD per year)"),
                               (M,"B  Cognitive progression (MoCA)","β  time × brain-PAD  (points per SD per year)")]):
    D=D.iloc[::-1].reset_index(drop=True); y=np.arange(len(D))
    for i in range(len(D)):
        nom=D['p'][i]<0.05
        col=ORANGE if nom else GREY
        ax.plot([D['lo'][i],D['hi'][i]],[i,i],color=col,lw=1.7,solid_capstyle='round')
        ax.scatter(D['beta_timePAD'][i],i,s=34 if nom else 24,color=col,zorder=3,edgecolor='white',lw=0.6)
    ax.axvline(0,color='#444',ls=(0,(4,3)),lw=0.9)
    ax.set_yticks(y); ax.set_yticklabels(D['network'],fontsize=8.5)
    ax.set_xlabel(xl,fontsize=8); ax.set_title(ttl,fontsize=9.5,loc='left',fontweight='bold',pad=8)
axes[0].tick_params(axis='y',length=0)
fig.suptitle("Baseline regional brain age does not predict clinical progression  (all FDR n.s.; n≈195, ~3 y)",
             fontsize=10,fontweight='bold',x=0.02,ha='left',y=1.02)
# 图例
from matplotlib.lines import Line2D
axes[0].legend([Line2D([0],[0],color=ORANGE,marker='o',lw=1.7,ls='-'),
                Line2D([0],[0],color=GREY,marker='o',lw=1.7,ls='-')],
               ['nominal P<0.05 (n.s. after FDR)','P≥0.05'],loc='upper right',fontsize=7,frameon=False)
save(fig,"Figure_S5_Longitudinal_Progression")
print("ALL DONE")
