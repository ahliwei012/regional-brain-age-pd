# -*- coding: utf-8 -*-
"""多面板主图 (nature-figure规范): PD区域脑龄与临床结局。"""
import os; os.environ["MPLCONFIGDIR"]=r"G:/Temp/mpl"
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import matplotlib.gridspec as gridspec
import statsmodels.api as sm

# ── 强制可编辑SVG文本 + 排版 ──
plt.rcParams['font.family']='sans-serif'
plt.rcParams['font.sans-serif']=['Arial','DejaVu Sans','Liberation Sans']
plt.rcParams['svg.fonttype']='none'; plt.rcParams['pdf.fonttype']=42
plt.rcParams.update({'font.size':8,'axes.spines.right':False,'axes.spines.top':False,
                     'axes.linewidth':0.9,'legend.frameon':False,'xtick.major.width':0.9,
                     'ytick.major.width':0.9,'xtick.major.size':3,'ytick.major.size':3})
P={'blue':'#0F4D92','blue_soft':'#3775BA','red':'#B64342','red_soft':'#E9A6A1',
   'grey':'#9A9A9A','grey_l':'#CFCECE','dark':'#272727','teal':'#42949E'}
def stars(p): return '***' if p<0.001 else '**' if p<0.01 else '*' if p<0.05 else ''
def panel(ax,l,x=-0.20,y=1.04): ax.text(x,y,l,transform=ax.transAxes,fontsize=11,fontweight='bold',va='bottom',ha='left')

BASE=r"F:/PD_brainage_data"; FIG=f"{BASE}/figures_nature"; os.makedirs(FIG,exist_ok=True)
U=pd.read_csv(f"{BASE}/results_final_updrs.csv"); M=pd.read_csv(f"{BASE}/results_final_moca.csv")
PVH=pd.read_csv(f"{BASE}/results_final_pdvshc.csv"); df=pd.read_csv(f"{BASE}/analysis_merged_final.csv")

fig=plt.figure(figsize=(7.2,6.4))
gs=gridspec.GridSpec(2,2,figure=fig,hspace=0.5,wspace=0.55,width_ratios=[1.12,1],height_ratios=[1,1],
                     left=0.13,right=0.95,top=0.93,bottom=0.10)

# ── a: HERO 森林图 UPDRS ──
axA=fig.add_subplot(gs[0,0]); s=U.sort_values('beta').reset_index(drop=True)
for i in range(len(s)):
    sig=s['fdr'][i]<0.05; col=P['red'] if sig else P['grey']
    axA.plot([s['lo'][i],s['hi'][i]],[i,i],color=col,lw=1.6,solid_capstyle='round',zorder=1)
    axA.scatter(s['beta'][i],i,s=34 if sig else 22,color=col,zorder=2,edgecolor='white',linewidth=0.6)
    if stars(s['fdr'][i]): axA.text(s['hi'][i]+0.06,i,stars(s['fdr'][i]),va='center',color=P['red'],fontsize=9,fontweight='bold')
axA.axvline(0,color=P['grey'],ls=(0,(4,3)),lw=0.9)
axA.set_yticks(range(len(s))); axA.set_yticklabels(s['name'],fontsize=8)
axA.set_xlabel('Standardized β  (older regional brain age → higher UPDRS-III)',fontsize=8)
axA.tick_params(axis='y',length=0); axA.set_xlim(-0.9,2.3)
panel(axA,'a')
axA.set_title('Regional brain age vs motor severity',fontsize=9,pad=6,loc='left',fontweight='bold')

# ── b: 热图 区域×{UPDRS,MoCA} ──
axB=fig.add_subplot(gs[0,1]); order=list(U['r'])
ub=U.set_index('r').loc[order]; mb=M.set_index('r').loc[order]
mat=np.column_stack([ub['beta'].values,mb['beta'].values]); vmax=np.abs(mat).max()
im=axB.imshow(mat,cmap='RdBu_r',norm=TwoSlopeNorm(0,-vmax,vmax),aspect='auto')
axB.set_xticks([0,1]); axB.set_xticklabels(['UPDRS-III','MoCA'],fontsize=8)
axB.set_yticks(range(len(order))); axB.set_yticklabels(ub['name'],fontsize=7.5)
axB.tick_params(length=0)
for i in range(len(order)):
    if stars(ub['fdr'].iloc[i]): axB.text(0,i,stars(ub['fdr'].iloc[i]),ha='center',va='center',fontsize=9,color=P['dark'])
    if stars(mb['fdr'].iloc[i]): axB.text(1,i,stars(mb['fdr'].iloc[i]),ha='center',va='center',fontsize=9,color=P['dark'])
for sp in axB.spines.values(): sp.set_visible(False)
axB.set_xticks(np.arange(-.5,2,1),minor=True); axB.set_yticks(np.arange(-.5,len(order),1),minor=True)
axB.grid(which='minor',color='white',lw=1.2); axB.tick_params(which='minor',length=0)
cb=fig.colorbar(im,ax=axB,fraction=0.05,pad=0.06); cb.set_label('Standardized β',fontsize=7.5); cb.ax.tick_params(labelsize=7)
panel(axB,'b',x=-0.32); axB.set_title('Motor vs cognitive dissociation',fontsize=9,pad=6,loc='left',fontweight='bold')

# ── c: PD vs HC 森林图 ──
axC=fig.add_subplot(gs[1,0]); s=PVH.sort_values('beta').reset_index(drop=True)
for i in range(len(s)):
    sig=s['fdr'][i]<0.05; col=P['blue'] if sig else P['grey']
    axC.plot([s['lo'][i],s['hi'][i]],[i,i],color=col,lw=1.6,solid_capstyle='round',zorder=1)
    axC.scatter(s['beta'][i],i,s=34 if sig else 22,color=col,zorder=2,edgecolor='white',linewidth=0.6)
    if stars(s['fdr'][i]): axC.text(s['lo'][i]-0.08,i,stars(s['fdr'][i]),va='center',ha='right',color=P['blue'],fontsize=9,fontweight='bold')
axC.axvline(0,color=P['grey'],ls=(0,(4,3)),lw=0.9)
axC.set_yticks(range(len(s))); axC.set_yticklabels(s['name'],fontsize=8); axC.tick_params(axis='y',length=0)
axC.set_xlabel('Brain-PAD difference  (PD − HC, years)',fontsize=8)
panel(axC,'c'); axC.set_title('PD vs HC regional brain age',fontsize=9,pad=6,loc='left',fontweight='bold')

# ── d: 散点 视觉PAD vs UPDRS + 95%CI带 ──
axD=fig.add_subplot(gs[1,1]); pdf=df[df.group=='PD'].dropna(subset=['vsl_PAD','updrs3_score'])
x=pdf['vsl_PAD'].values; y=pdf['updrs3_score'].values
axD.scatter(x,y,s=9,alpha=0.35,color=P['blue_soft'],edgecolor='none',rasterized=True)
xs=np.linspace(np.percentile(x,1),np.percentile(x,99),80)
Xd=sm.add_constant(x); res=sm.OLS(y,Xd).fit(); pr=res.get_prediction(sm.add_constant(xs)).summary_frame(alpha=0.05)
axD.plot(xs,pr['mean'],color=P['red'],lw=1.8)
axD.fill_between(xs,pr['mean_ci_lower'],pr['mean_ci_upper'],color=P['red_soft'],alpha=0.4,linewidth=0)
r=np.corrcoef(x,y)[0,1]
axD.text(0.04,0.94,f'r = {r:.2f},  n = {len(pdf)}',transform=axD.transAxes,fontsize=8,va='top')
axD.set_xlabel('Visual network brain-PAD (years)',fontsize=8); axD.set_ylabel('UPDRS-III',fontsize=8)
panel(axD,'d',x=-0.24); axD.set_title('Exemplar network association',fontsize=9,pad=6,loc='left',fontweight='bold')

for f in ['svg','pdf','png']:
    fig.savefig(f"{FIG}/Figure_PD_brainage.{f}",dpi=600 if f=='png' else None,bbox_inches='tight')
fig.savefig(f"{FIG}/Figure_PD_brainage.tiff",dpi=600,bbox_inches='tight')
plt.close(fig)
print("已生成:", sorted(os.listdir(FIG)))
