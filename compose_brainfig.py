# -*- coding: utf-8 -*-
"""三面板主图: a=UPDRS脑图(红), b=MoCA脑图(蓝), c=森林图(区域脑龄→UPDRS-III)。"""
import os; os.environ["MPLCONFIGDIR"]=r"G:/Temp/mpl"
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from PIL import Image
plt.rcParams['font.family']='sans-serif'; plt.rcParams['font.sans-serif']=['Arial','DejaVu Sans']
plt.rcParams['svg.fonttype']='none'; plt.rcParams['pdf.fonttype']=42
plt.rcParams.update({'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':0.9})
FIG=r"F:/PD_brainage_data/figures_nature"
def star(p): return '***' if p<0.001 else '**' if p<0.01 else '*' if p<0.05 else ''

def load_crop(p,pad=12):
    im=np.array(Image.open(p).convert("RGB")); mask=(im<245).any(2); ys,xs=np.where(mask)
    return im[max(0,ys.min()-pad):min(im.shape[0],ys.max()+pad), max(0,xs.min()-pad):min(im.shape[1],xs.max()+pad)]

A=load_crop(f"{FIG}/updrs_beta_sig.png"); B=load_crop(f"{FIG}/moca_beta_sig.png")
arA,arB=A.shape[0]/A.shape[1], B.shape[0]/B.shape[1]
U=pd.read_csv(r"F:/PD_brainage_data/results_final_updrs.csv").sort_values('beta').reset_index(drop=True)

W=6.8
fig=plt.figure(figsize=(W, W*(arA+arB)+3.0))
gs=gridspec.GridSpec(3,1,figure=fig,hspace=0.16,height_ratios=[arA,arB,0.60])

for gi,(img,lab,sub) in enumerate([(A,'A','Older regional brain age  →  worse motor (MDS-UPDRS-III)'),
                                    (B,'B','Older regional brain age  →  worse cognition (MoCA)')]):
    ax=fig.add_subplot(gs[gi]); ax.imshow(img); ax.axis('off')
    ax.text(-0.02,1.015,lab,transform=ax.transAxes,fontsize=18,fontweight='bold',va='bottom',ha='left')
    ax.text(0.5,1.02,sub,transform=ax.transAxes,fontsize=9.5,fontweight='bold',ha='center',va='bottom',color='#272727')

# c 森林图
axc=fig.add_subplot(gs[2]); s=U
for i in range(len(s)):
    sig=s['fdr'][i]<0.05; col='#c0392b' if sig else '#9A9A9A'
    axc.plot([s['lo'][i],s['hi'][i]],[i,i],color=col,lw=1.7,solid_capstyle='round')
    axc.scatter(s['beta'][i],i,s=36 if sig else 22,color=col,zorder=3,edgecolor='white',linewidth=0.6)
    if star(s['fdr'][i]): axc.text(s['hi'][i]+0.05,i,star(s['fdr'][i]),va='center',color='#c0392b',fontsize=10,fontweight='bold')
axc.axvline(0,color='#9A9A9A',ls=(0,(4,3)),lw=0.9)
axc.set_yticks(range(len(s))); axc.set_yticklabels(s['name'],fontsize=8.5); axc.tick_params(axis='y',length=0)
axc.set_xlabel('Standardized β  (regional brain-PAD → UPDRS-III)',fontsize=9)
axc.set_xlim(-0.9,2.35)
axc.text(-0.02,1.03,'C',transform=axc.transAxes,fontsize=18,fontweight='bold',va='bottom',ha='left')
axc.set_title('Effect sizes across networks (LME + scanner random effect)',fontsize=9.5,loc='center',pad=6,fontweight='bold')

for f in ['png','pdf','svg']:
    fig.savefig(f"{FIG}/Figure_main_composite.{f}",dpi=400 if f=='png' else None,bbox_inches='tight')
plt.close(fig); print("已生成 Figure_main_composite.png/pdf/svg")
