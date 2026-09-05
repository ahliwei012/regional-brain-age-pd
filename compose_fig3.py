# -*- coding: utf-8 -*-
"""Figure 3: a=PD vs HC脑图, b=森林图。"""
import os; os.environ["MPLCONFIGDIR"]=r"G:/Temp/mpl"
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.gridspec as gridspec
from PIL import Image
plt.rcParams['font.family']='sans-serif'; plt.rcParams['font.sans-serif']=['Arial','DejaVu Sans']
plt.rcParams['svg.fonttype']='none'; plt.rcParams['pdf.fonttype']=42
plt.rcParams.update({'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':0.9})
P={'blue':'#0F4D92','grey':'#9A9A9A'}; MF=r"F:/PD_brainage_data/manuscript_figures"
def star(p): return '***' if p<0.001 else '**' if p<0.01 else '*' if p<0.05 else ''
def load_crop(p,pad=12):
    im=np.array(Image.open(p).convert("RGB")); m=(im<245).any(2); ys,xs=np.where(m)
    return im[max(0,ys.min()-pad):min(im.shape[0],ys.max()+pad),max(0,xs.min()-pad):min(im.shape[1],xs.max()+pad)]

A=load_crop(f"{MF}/pdvshc_beta_sig.png"); ar=A.shape[0]/A.shape[1]
PVH=pd.read_csv(r"F:/PD_brainage_data/results_final_pdvshc.csv").sort_values('beta').reset_index(drop=True)
W=6.8; fig=plt.figure(figsize=(W, W*ar+3.4))
gs=gridspec.GridSpec(2,1,figure=fig,height_ratios=[ar,0.62],hspace=0.20)

axa=fig.add_subplot(gs[0]); axa.imshow(A); axa.axis('off')
axa.text(-0.02,1.075,'A',transform=axa.transAxes,fontsize=18,fontweight='bold',va='bottom',ha='left')
axa.text(0.5,1.02,'PD vs HC: younger-appearing brain age (sensorimotor, limbic)',
         transform=axa.transAxes,fontsize=9.5,fontweight='bold',ha='center',va='bottom',color='#272727')

axb=fig.add_subplot(gs[1]); s=PVH
for i in range(len(s)):
    sig=s['fdr'][i]<0.05; col=P['blue'] if sig else P['grey']
    axb.plot([s['lo'][i],s['hi'][i]],[i,i],color=col,lw=1.7,solid_capstyle='round')
    axb.scatter(s['beta'][i],i,s=36 if sig else 22,color=col,zorder=3,edgecolor='white',linewidth=0.6)
    if star(s['fdr'][i]): axb.text(s['lo'][i]-0.06,i,star(s['fdr'][i]),va='center',ha='right',color=P['blue'],fontsize=10,fontweight='bold')
axb.axvline(0,color=P['grey'],ls=(0,(4,3)),lw=0.9)
axb.set_yticks(range(len(s))); axb.set_yticklabels(s['name'],fontsize=8.5); axb.tick_params(axis='y',length=0)
axb.set_xlabel('Regional brain-PAD difference  (PD − HC, years)',fontsize=9)
axb.text(-0.02,1.03,'B',transform=axb.transAxes,fontsize=18,fontweight='bold',va='bottom',ha='left')
axb.set_title('Effect sizes across networks  (LME + scanner random effect; blue = FDR<0.05)',fontsize=9,loc='center',pad=6,fontweight='bold')

for f in ['png','pdf','svg']:
    fig.savefig(f"{MF}/Figure_3_PD_vs_HC.{f}",dpi=400 if f=='png' else None,bbox_inches='tight')
plt.close(fig); print("Figure_3_PD_vs_HC 更新为 脑图+森林图 双面板")
