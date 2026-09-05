# -*- coding: utf-8 -*-
"""Figure 5: 皮层网络脑龄 vs 多巴胺能失神经(DaTSCAN SBR)
   a) 解耦: 各网络 brain-PAD ~ 壳核SBR (全部n.s.)
   b) 增量效度: 控制SBR前后 brain-PAD 的β
   c) 方差分解: 脑龄独有 / DAT独有 / 共享"""
import os; os.environ["MPLCONFIGDIR"]=r"G:/Temp/mpl"
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import warnings; warnings.filterwarnings("ignore")
plt.rcParams['font.family']='sans-serif'; plt.rcParams['font.sans-serif']=['Arial','DejaVu Sans']
plt.rcParams['svg.fonttype']='none'; plt.rcParams['pdf.fonttype']=42
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':0.9})
BASE=r"F:/PD_brainage_data"; OUT=f"{BASE}/manuscript_figures"
BLUE='#0F4D92'; ORANGE='#D55E00'; GREY='#9A9A9A'; GREEN='#1B7837'; PURP='#762A83'
REG=["sns","fpn","drs","vnt","dft","slt","lng","adt","vsl","lmb"]
EN={"sns":"Sensorimotor","fpn":"Frontoparietal","drs":"Dorsal attention","vnt":"Ventral attention",
    "dft":"Default mode","slt":"Salience","lng":"Language","adt":"Auditory","vsl":"Visual","lmb":"Limbic"}

sb=pd.read_csv(f'{BASE}/DaTSCAN/Xing_Core_Lab_-_Quant_SBR_24Jul2026.csv')
order={'SC':0,'BL':1,'V04':2,'V06':3,'V10':4}
sb['ord']=sb.EVENT_ID.map(order).fillna(9)
sb=sb.sort_values('ord').drop_duplicates('PATNO')[['PATNO','PUTAMEN_REF_CWM']].rename(columns={'PUTAMEN_REF_CWM':'putamen'})
df=pd.read_csv(f'{BASE}/analysis_with_thickness.csv').merge(sb,on='PATNO',how='inner')
pdf=df[df.group=='PD'].copy()

# ---- Panel a: 解耦 ----
rows=[]
for r in REG:
    d=pdf.dropna(subset=[f'{r}_PAD','putamen','age_at_visit','SEX','EDUCYRS','batch'])
    m=smf.mixedlm(f"{r}_PAD ~ putamen + age_at_visit + C(SEX) + EDUCYRS",d,groups=d['batch']).fit()
    ci=m.conf_int().loc['putamen']
    rows.append([EN[r],m.params['putamen'],ci[0],ci[1],m.pvalues['putamen']])
A=pd.DataFrame(rows,columns=['name','beta','lo','hi','p']); A['fdr']=multipletests(A['p'],method='fdr_bh')[1]
A.to_csv(f'{BASE}/DAT_coupling.csv',index=False,encoding='utf-8-sig')

B=pd.read_csv(f'{BASE}/DAT_incremental_validity.csv')
V=pd.read_csv(f'{BASE}/DAT_variance_partition.csv')

fig=plt.figure(figsize=(11.4,3.9))
gs=fig.add_gridspec(1,3,width_ratios=[1.05,1.0,0.95],wspace=0.52)

# a) 解耦森林图
ax=fig.add_subplot(gs[0]); Ar=A.iloc[::-1].reset_index(drop=True); y=np.arange(len(Ar))
for i in range(len(Ar)):
    ax.plot([Ar['lo'][i],Ar['hi'][i]],[i,i],color=GREY,lw=1.8,solid_capstyle='round')
    ax.scatter(Ar['beta'][i],i,s=28,color=GREY,zorder=3,edgecolor='white',lw=0.6)
ax.axvline(0,color='black',ls=(0,(4,3)),lw=0.9)
ax.set_yticks(y); ax.set_yticklabels(Ar['name'],fontsize=8.2); ax.tick_params(axis='y',length=0)
ax.set_xlabel('β  regional brain-PAD per unit putamen SBR',fontsize=8.2,color='black')
ax.set_title('A  Cortical brain age is decoupled\n    from dopaminergic denervation',fontsize=9.5,loc='left',fontweight='bold',color='black')
ax.text(0.98,0.03,'all FDR > 0.47',transform=ax.transAxes,ha='right',fontsize=8,color='black',style='italic')

# b) 增量效度
ax=fig.add_subplot(gs[1])
OCN={'updrs3_score':'UPDRS-III','moca':'MoCA','UPDRS-III':'UPDRS-III','MoCA':'MoCA'}
B['ocn']=B['outcome'].map(OCN)
lab=[f"{r['network']}\n({r['ocn']})" for _,r in B.iterrows()]
yp=np.arange(len(B))[::-1]
for i,(yi,(_,r)) in enumerate(zip(yp,B.iterrows())):
    ax.scatter(abs(r['b_PADonly']),yi+0.16,s=40,color=BLUE,zorder=3,edgecolor='white',lw=0.6)
    ax.scatter(abs(r['b_PADadjDAT']),yi-0.16,s=40,color=ORANGE,marker='D',zorder=3,edgecolor='white',lw=0.6)
    ax.plot([abs(r['b_PADonly']),abs(r['b_PADadjDAT'])],[yi+0.16,yi-0.16],color='#BBB',lw=0.8,zorder=1)
    st='*' if r['p_PADadjDAT']<0.05 else ''
    ax.text(1.02,yi,st,transform=ax.get_yaxis_transform(),va='center',fontsize=11,color='black')
ax.set_yticks(yp); ax.set_yticklabels(lab,fontsize=7.6); ax.tick_params(axis='y',length=0)
ax.set_xlabel('|β|  clinical score per SD of brain-PAD',fontsize=8.2,color='black')
ax.set_title('B  Association survives adjustment\n    for dopaminergic denervation',fontsize=9.5,loc='left',fontweight='bold',color='black')
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([0],[0],marker='o',color=BLUE,ls='',label='brain-PAD only'),
                   Line2D([0],[0],marker='D',color=ORANGE,ls='',label='brain-PAD | DAT')],
          loc='lower right',fontsize=7,frameon=False)

# c) 方差分解
ax=fig.add_subplot(gs[2])
grp=V.groupby('outcome')[['unique_PAD','unique_SBR','shared']].mean()*100
grp=grp.reindex(['UPDRS-III','MoCA'])
xs=np.arange(len(grp)); w=0.55
b1=ax.bar(xs,grp['unique_PAD'],w,color=BLUE,label='Brain age (unique)')
b2=ax.bar(xs,grp['unique_SBR'],w,bottom=grp['unique_PAD'],color=ORANGE,label='DAT-SPECT (unique)')
b3=ax.bar(xs,grp['shared'].clip(lower=0),w,bottom=grp['unique_PAD']+grp['unique_SBR'],color=GREY,label='Shared')
ax.set_xticks(xs); ax.set_xticklabels(grp.index,fontsize=9)
ax.set_ylabel('Unique variance explained (ΔR², %)',fontsize=8.2,color='black')
ax.set_title('C  Variance partition',fontsize=9.5,loc='left',fontweight='bold',color='black')
ax.legend(fontsize=7,frameon=False,loc='upper right')
for xi,(pa,sb_) in enumerate(zip(grp['unique_PAD'],grp['unique_SBR'])):
    ax.text(xi,pa/2,f'{pa:.2f}',ha='center',va='center',fontsize=7.5,color='white',fontweight='bold')
    if sb_>0.15: ax.text(xi,pa+sb_/2,f'{sb_:.2f}',ha='center',va='center',fontsize=7.5,color='white',fontweight='bold')
    else: ax.text(xi,pa+sb_+0.08,f'{sb_:.2f}',ha='center',va='bottom',fontsize=7.5,color=ORANGE)
ax.tick_params(colors='black')

fig.suptitle("Regional cortical brain age carries clinical information independent of dopaminergic denervation",
             fontsize=10.5,fontweight='bold',x=0.02,ha='left',y=1.04,color='black')
for f in ['png','pdf','svg']:
    fig.savefig(f"{OUT}/Figure_5_DAT_independence.{f}",dpi=400 if f=='png' else None,bbox_inches='tight')
plt.close(fig)
print("Figure_5_DAT_independence done")
print(A.round(4).to_string(index=False))
