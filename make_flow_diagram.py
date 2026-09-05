# -*- coding: utf-8 -*-
"""参与者流程图 (STROBE item 13): 提交CIVET -> 失败 -> 成功 -> 缺临床变量排除 -> 最终分析样本。
   数字来自真实日志(civet_pd_all.log/civet_hc_batch.log/civet_newhc.log)与最终分析表交叉核对。"""
import os; os.environ["MPLCONFIGDIR"]=r"G:/Temp/mpl"
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
plt.rcParams['font.family']='sans-serif'; plt.rcParams['font.sans-serif']=['Arial','DejaVu Sans']
plt.rcParams['svg.fonttype']='none'; plt.rcParams['pdf.fonttype']=42
OUT=r"F:/PD_brainage_data/manuscript_figures"
INK='#111111'; BLUE='#0F4D92'; RED='#C1362F'; GREY='#7A7A7A'; BOXFC='white'

fig=plt.figure(figsize=(8.6,4.3)); ax=fig.add_axes([0,0,1,1])
ax.set_xlim(0,100); ax.set_ylim(42,103); ax.axis('off')

def box(x,y,w,h,text,ec=INK,fs=9.2):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.3,rounding_size=1.4",
                 lw=1.6,edgecolor=ec,facecolor=BOXFC))
    ax.text(x+w/2,y+h/2,text,ha='center',va='center',fontsize=fs,color=INK,linespacing=1.35)

def arrow(x0,y0,x1,y1,color=INK):
    ax.annotate("",xy=(x1,y1),xytext=(x0,y0),arrowprops=dict(arrowstyle='-|>',color=color,lw=1.5))

def sidebox(x,y,w,h,text,color=RED):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.25,rounding_size=1.2",
                 lw=1.3,edgecolor=color,facecolor='#FBEEEE' if color==RED else '#EEF3FA'))
    ax.text(x+w/2,y+h/2,text,ha='center',va='center',fontsize=8.2,color=color,linespacing=1.3)

# 顶部: 提交处理
box(28,88,44,9,"PPMI participants with a baseline T1-weighted scan\nsubmitted for CIVET cortical processing\n(560 PD; 282 controls; n = 842)",fs=8.8)
arrow(50,88,50,80)

# CIVET失败(侧支)
sidebox(70,79.5,26,8,"Excluded: CIVET processing\nfailed to yield cortical thickness\n(10 PD; 17 controls; n = 27)")
arrow(58,83.5,70,83.5,color=RED)

box(28,71,44,9,"Successfully processed with valid\ncortical thickness and GM/WM features\n(550 PD; 265 controls; n = 815)",fs=8.8)
arrow(50,71,50,63)

# 缺临床变量(侧支)
sidebox(70,62.5,26,8,"Excluded: missing MDS-UPDRS-III,\nMoCA or other required clinical\nvariables (8 PD; 4 controls; n = 12)")
arrow(58,66.5,70,66.5,color=RED)

box(28,54,44,9,"Final analytic sample\n542 patients with Parkinson's disease\n261 healthy controls (n = 803)",ec=BLUE,fs=9.3)

ax.text(50,97,"Participant flow",fontsize=12.5,fontweight='bold',ha='center',color=INK)
ax.text(50,48,"Numbers reconciled across CIVET processing logs, the feature-extraction output\n"
             "directory (815 files), and the final analytic dataset (803 rows).",
        ha='center',fontsize=7.8,style='italic',color='#555')

for f in ['png','pdf']:
    fig.savefig(f"{OUT}/Figure_S6_Participant_Flow.{f}",dpi=400 if f=='png' else None,bbox_inches='tight',facecolor='white')
plt.close(fig)
print("Figure_S6_Participant_Flow done")
