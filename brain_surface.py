# -*- coding: utf-8 -*-
"""真实皮层表面渲染: 把区域脑龄→UPDRS-III效应(标准化β)投射到大脑表面(nilearn)。"""
import os; os.environ["MPLCONFIGDIR"]=r"G:/Temp/mpl"
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from nilearn import plotting

plt.rcParams['font.family']='sans-serif'; plt.rcParams['font.sans-serif']=['Arial','DejaVu Sans']
plt.rcParams['svg.fonttype']='none'; plt.rcParams['pdf.fonttype']=42

BASE=r"F:/PD_brainage_data"; SURF=r"F:/Deep_learning_derived_ MRI regional brain age/regional_Brain_age-main/surface_information"
CIV=f"{BASE}/civet_out/100005"; FIG=f"{BASE}/figures_nature"

# 网络 -> AAL区号 (来自模型源码)
NET={"sns":[1,2,19,20,57,58,69,70],"fpn":[5,6,7,8,9,10,65,66],"drs":[3,4,59,60],
     "vnt":[11,12,13,14,63,64],"dft":[21,22,25,26,27,28,35,36,65,66,67,68,85,86],
     "slt":[29,30,31,32],"lng":[11,12,13,17,63],"adt":[79,80,81,82],
     "vsl":[43,44,45,46,47,48,49,50,51,52,53,54,55,56,89,90],
     "lmb":[15,16,23,24,33,34,39,40,83,84,87,88]}

def read_obj(path):
    t=open(path).read().split(); assert t[0]=='P'; n=int(t[6]); i=7
    coords=np.array(t[i:i+3*n],float).reshape(n,3); i+=3*n; i+=3*n     # coords, skip normals
    nit=int(t[i]); i+=1; cf=int(t[i]); i+=1
    i+= 4 if cf==0 else (4*nit if cf==1 else 4*n)                       # skip colour
    i+=nit                                                              # skip end_indices
    faces=np.array(t[i:i+3*nit],int).reshape(nit,3)
    return coords, faces

# UPDRS 标准化β per network
U=pd.read_csv(f"{BASE}/results_final_updrs.csv").set_index("r")
beta={k:U.loc[k,"beta"] for k in NET}; fdr={k:U.loc[k,"fdr"] for k in NET}

aal=np.loadtxt(f"{SURF}/aal_atlas.txt").astype(int)   # 81924 = L40962 + R40962
NF=40962
def vertex_beta(aal_h):
    v=np.full(len(aal_h), np.nan)
    for k,regs in NET.items():
        m=np.isin(aal_h, regs)
        # 取该顶点所属网络中β最大者(突出显著网络)
        cur=v.copy(); nv=np.where(np.isnan(cur)|(beta[k]>cur), beta[k], cur)
        v[m]=np.where(np.isnan(v[m]) | (beta[k]>np.nan_to_num(v[m],nan=-9)), beta[k], v[m])
    return v

# 只高亮 FDR<0.05 的显著网络, 其余顶点=NaN(灰)
SIG=[k for k in NET if fdr[k]<0.05]
print("显著网络:", SIG)
def vbeta(aal_h):
    best=np.full(len(aal_h), -np.inf)
    for k in SIG:
        m=np.isin(aal_h,NET[k]); best[m]=np.maximum(best[m],beta[k])
    best[best==-np.inf]=np.nan; return best

meshes={"left":read_obj(f"{CIV}/PPMI_100005_mid_surface_rsl_left_81920.obj"),
        "right":read_obj(f"{CIV}/PPMI_100005_mid_surface_rsl_right_81920.obj")}
maps={"left":vbeta(aal[:NF]),"right":vbeta(aal[NF:2*NF])}
print("L顶点:",len(meshes['left'][0]),"面:",len(meshes['left'][1]),"有效映射:",np.sum(~np.isnan(maps['left'])))

# 红色系colormap(高β=红=脑龄老→运动差), NaN显灰
vmax=np.nanmax([np.nanmax(maps['left']),np.nanmax(maps['right'])])
cmap=LinearSegmentedColormap.from_list("mot",["#F2E6E0","#E9A6A1","#B64342","#6E1F1E"])

fig=plt.figure(figsize=(6.6,4.2))
views=[("left","lateral"),("left","medial"),("right","lateral"),("right","medial")]
titles=["Left lateral","Left medial","Right lateral","Right medial"]
for idx,((hemi,view),tit) in enumerate(zip(views,titles)):
    ax=fig.add_subplot(2,2,idx+1,projection='3d')
    c,f=meshes[hemi]
    bg=np.full(len(c),0.62)                       # 均匀浅灰背景, 消除内侧壁网格
    plotting.plot_surf_stat_map([c,f], maps[hemi], bg_map=bg, hemi=hemi, view=view, colorbar=False,
                                cmap=cmap, vmax=vmax, threshold=0.001, bg_on_data=False,
                                darkness=1.0, axes=ax, engine="matplotlib")
    ax.set_title(tit,fontsize=8,pad=-2)
    try: ax.dist=7.2
    except Exception: pass
fig.subplots_adjust(left=0.02,right=0.98,top=0.90,bottom=0.14,wspace=-0.05,hspace=-0.05)
fig.suptitle("Networks whose brain age tracks motor severity in PD (UPDRS-III)",fontsize=10,fontweight='bold',y=0.99)
fig.text(0.5,0.53,"Significant networks (FDR < 0.05)\nVisual · Ventral attention · Salience · Language",
         ha='center',va='center',fontsize=9,color='#4D4D4D',linespacing=1.5)
# 单独colorbar
cax=fig.add_axes([0.4,0.06,0.22,0.025])
import matplotlib as mpl
cb=mpl.colorbar.ColorbarBase(cax,cmap=cmap,norm=mpl.colors.Normalize(0,vmax),orientation='horizontal')
cb.set_label("Standardized β (older brain age → worse motor)",fontsize=7); cb.ax.tick_params(labelsize=6.5)
for f in ['png','svg','pdf']: fig.savefig(f"{FIG}/Figure_brain_surface.{f}",dpi=600 if f=='png' else None,bbox_inches='tight')
plt.close(fig); print("已生成 Figure_brain_surface.*")
